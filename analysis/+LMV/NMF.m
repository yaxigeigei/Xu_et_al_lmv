classdef NMF
    
    methods(Static)
        function sBatch = BatchClustering(ce, kList, regionGroups, evalMethod)
            % Run bootstrap NMF clustering with different number of components and regional groupings
            % 
            %   sBatch = BatchClustering(ce, kList)
            %   sBatch = BatchClustering(ce, kList, regionGroups)
            %   sBatch = BatchClustering(ce, kList, regionGroups, evalMethod)
            % 
            % Inputs
            %   ce              A NP.CodingExplorer object.
            %   kList           A vector for the numbers of components (or factors) to compute for.
            %   regionGroups    A cell array of string vectors. Each vector includes the regions of interest 
            %                   for a clustering. Default is a single group including all regions.
            % Outputs
            %   sBatch          A struct array of output. Each element corresponds to a region group.
            %     ce            A duplicate of input ce but only with units from the regions of interest. The 
            %                   clusTb in ce has cluster info added by LMV.NMF.ComputeOptimalComp.
            %     sBoot         A struct array of bootstrap clustering reuslts returned by LMV.NMF.BootClustering.
            %     sOpti         
            % 
            
            if nargin < 4 || isempty(evalMethod)
                evalMethod = 'gap';
            end
            
            % Make sure regionGroups is a cell array of string vectors
            if ~exist('regionGroups', 'var') || isempty(regionGroups)
                regionGroups = "all";
            end
            if iscellstr(regionGroups)
                regionGroups = string(regionGroups);
            end
            if ~iscell(regionGroups)
                regionGroups = num2cell(regionGroups);
            end
            
            % Iterate through region groups
            cTb = ce.clusTb;
            sBatch = cell(size(regionGroups));
            for i = 1 : numel(regionGroups)
                % Create a region specific ce
                fprintf("\nCompute bootstrap NMF clustering for %s\n", strjoin(regionGroups{i}, ', '));
                if regionGroups{i} == "all"
                    rm = [];
                else
                    rm = find(~ismember(cTb.region, regionGroups{i}));
                end
                ceSub = ce.Duplicate;
                tsNames = {'resp', 'sd', 'sem', 'spikeRate'};
                for j = 1 : numel(tsNames)
                    ceSub.SetColumn(tsNames{j}, rm+1, []);
                end
                etNames = {'spikeTime'};
                for j = 1 : numel(etNames)
                    ceSub.SetColumn(etNames{j}, rm, []);
                end
                ceSub.clusTb(rm,:) = [];
                
                % Get unit PETHs
                [~, mm] = ceSub.GetArray('resp', 'Normalization', 'minmaxsoft');
                
                s = struct;
                switch evalMethod
                    case 'stability'
                        % Iterate through the numbers of components
                        sClus = cell(size(kList));
                        for j = 1 : numel(kList)
                            fprintf("# of components: %i\n", kList(j));
                            sClus{j} = NP.BootStability(mm, 'nComp', kList(j), 'nBoot', 1e2, 'fraction', .5);
                        end
                        sClus = cat(1, sClus{:});
                        
                        % Compute optimal number of components
                        [ceSub, sEval] = LMV.NMF.FindOptimalK(ceSub, sClus);
                        
                        % Pack outputs
                        s.region = regionGroups{i};
                        s.ce = ceSub;
                        s.sClus = sClus;
                        s.sEval = sEval;
                        
                    case {'CalinskiHarabasz', 'DaviesBouldin', 'silhouette', 'gap'}
                        % Iterate through the numbers of components
                        [sClus, sEval] = LMV.NMF.BootCriterion(mm, kList, evalMethod, 'nBoot', 1e2, 'fraction', .9);
                        
                        % Add cluster info to unit table
                        sOpti = sClus(kList==sEval.optiK);
                        ceSub.clusTb.nmfcId = sOpti.clustId;
                        ceSub.clusTb.nmfcW = max(sOpti.W, [], 1)';
                        
                        % Pack outputs
                        s.region = regionGroups{i};
                        s.ce = ceSub;
                        s.sClus = sClus;
                        s.sEval = sEval;
                        
                    otherwise
                        error("'%s' is not a valid method", evalMethod);
                end
                sBatch{i} = s;
            end
            sBatch = cat(1, sBatch{:});
        end
        
        function [sClus, sEval] = BootCriterion(X, kList, criterion, varargin)
            % Bootstrap clustering evaluated by Calinski-Harabasz criterion 
            % 
            %   [sClus, sEval] = BootCriterion(X, kList, criterion)
            %   [sClus, sEval] = BootCriterion(..., 'nBoot', 1e3)
            %   [sClus, sEval] = BootCriterion(..., 'fraction', .9)
            % 
            % Inputs
            %   X               t-by-u array, where t is the number of samples, u is the number of units.
            %   kList           A vector of the numbers of clusters to compute for.
            %   criterion       
            %   'nBoot'         The number of bootstrap iterations.
            %   'fraction'      The fraction of units to run bootstrap clustering. The rest are held out 
            %                   for testing.
            % Output
            %   Struct sClus with the following fields:
            %     nComp         The number of clusters computed, copied from 'nComp'.
            %     H             t-by-k array of template cluster centroid patterns (i.e. nonnegative 
            %                   left factor of X) computed from all units.
            %     W             k-by-u array of component weights (i.e. nonnegative right factor of X)
            %                   computed from all units.
            %     clustId       u-by-1 array of assgined cluster IDs.
            %     clustW        u-by-1 array of weights on the assigned clusters.
            % 
            %   Struct sEval with the following fields:
            %     kList         Same as the input kList.
            %     fullScore     k-by-1 array.
            %     fullError     k-by-1 array.
            %     fullOptiK     
            %     bootScore     k-by-nBoot array.
            %     bootError     k-by-nBoot array.
            %     bootOptiK     1-by-nBoot vector.
            %     optiK         
            % 
            
            p = inputParser();
            p.addParameter('nBoot', 1e2, @isscalar);
            p.addParameter('fraction', .9, @isscalar);
            p.parse(varargin{:});
            nBoot = p.Results.nBoot;
            frac = p.Results.fraction;
            
            if lower(criterion) == "gap"
                args = {'Distance', 'sqEuclidean', 'B', 20, 'ReferenceDistribution', 'uniform', 'SearchMethod', 'firstMaxSE'};
            else
                args = {};
            end
            
            % Cluster using full data
            rng(61);
            sClus = arrayfun(@(k) LMV.NMF.Clustering(X, k), kList);
            
            % Evaluate full data clustering
            evaObj = evalclusters(X', @LMV.NMF.ClusteringForEvalclusters, criterion, 'KList', kList, args{:});
            fullScore = evaObj.CriterionValues;
            fullError = NaN(size(fullScore));
            if isprop(evaObj, 'SE')
                fullError = evaObj.SE;
            end
            fullOptiK = evaObj.OptimalK;
            
            % Bootstrap evaluation
            [~, nUnit] = size(X);
            nSample = round(nUnit*frac);
            bootScore = NaN(numel(kList), nBoot);
            bootError = NaN(numel(kList), nBoot);
            bootOptiK = NaN(1, nBoot);
            
            sc = parallel.pool.Constant(RandStream('Threefry'));
            M = 10;
            parfor (n = 1 : nBoot, M)
                % Print progress
                str = sprintf('Iteration %i/%i\n', n, nBoot);
                fprintf(str);
                
                % Random sampling without replacement
                stream = sc.Value;
                stream.Substream = n;
                ind = randsample(stream, nUnit, nSample, false);
                ind = MMath.Ind2Logical(ind, nUnit);
                
                % Find optimal K
                evaObj = evalclusters(X(:,ind)', @LMV.NMF.ClusteringForEvalclusters, criterion, 'KList', kList, args{:});
                bootOptiK(n) = evaObj.OptimalK;
                bootScore(:,n) = evaObj.CriterionValues;
                if isprop(evaObj, 'SE')
                    bootError(:,n) = evaObj.SE;
                end
                
                % Clear printed
                if M == 0 && n < nBoot
                    fprintf(repmat('\b', 1, numel(str)));
                end
            end
            
            sEval.kList = kList;
            sEval.fullScore = fullScore;
            sEval.fullError = fullError;
            sEval.fullOptiK = fullOptiK;
            sEval.bootScore = bootScore;
            sEval.bootError = bootError;
            sEval.bootOptiK = bootOptiK;
            sEval.optiK = median(bootOptiK);
        end
        
        function C = ClusteringForEvalclusters(X, k)
            sClus = LMV.NMF.Clustering(X', k);
            C = sClus.clustId;
        end
        
        function s = Clustering(X, k)
            
            % Perform NMF, note that H is normalized to unit vectors
            [W, H] = nnmf(X', k);
            W = W'; % make it column vectors of weights, each for a neuron
            H = H'; % make it column vectors of patterns, each for a cluster
            
            % Sort components by peak time
            [~, I] = LMV.NMF.SortPatterns(H, 'peak');
            W = W(I,:);
            H = H(:,I);
            
            % Find cluster membership by the maximum coefficient
            [clustW, clustIdx] = max(W, [], 1);
            
            s.nComp = k;
            s.H = H;
            s.W = W;
            s.clustId = clustIdx';
            s.clustW = clustW';
        end
        
        function [pp, ind] = SortPatterns(pp, methodStr)
            % Sort patterns (or any timeseries) by certain criteria
            % 
            %   [pp, ind] = SortPatterns(pp, methodStr)
            % 
            if nargin < 2
                methodStr = '';
            end
            switch methodStr
                case 'sum'
                    [~, ind] = sort(sum(pp), 'descend');
                case 'csc'
                    [~, ind] = sort(sum(cumsum(pp)), 'descend');
                case 'cosine'
                    r_vect = pdist(pp', 'cosine');
                    r = squareform(r_vect);
                    r_sum = sum(r, 'omitnan');
                    [~, ind] = sort(r_sum, 'ascend');
                case 'dist'
                    r_vect = pdist(pp', 'hamming');
                    r = squareform(r_vect);
                    r_sum = sum(r, 'omitnan');
                    [~, ind] = sort(r_sum, 'descend');
                case 'peak'
                    [~, maxInd] = max(pp);
                    [~, ind] = sort(maxInd, 'ascend');
                case 'com'
                    cs = cumsum(pp);
                    cs = cs > max(cs)./2;
                    [~, ind] = sort(sum(cs), 'descend');
                otherwise
                    ind = 1 : size(pp,2);
            end
            pp = pp(:,ind);
        end
        
        function PlotBases(sBatch, k2Plot)
            % Plot components or bases from the H matrix for different component numbers
            
            if nargin < 2
                k2Plot = 3 : 10;
            end
            
            % Extract data
            sClus = sBatch.sClus;
            kList = [sClus.nComp];
            ce = sBatch.ce;
            t = ce.GetArray('resp');
            
            % Plotting
            tl = tiledlayout(1, numel(k2Plot));
            tl.Padding = 'compact';
            
            for i = 1 : numel(k2Plot)
                ax = nexttile(tl);
                
                k = kList == k2Plot(i);
                H = sClus(k).H * 5; % with scaling
                h = 1 : size(H,2);
                MPlot.PlotTraceLadder(t, H(:,:,1), h, 'Color', 'k', 'Parent', ax);
                
                colNames = {'cue1', 'stim', 'cue3', 'prod'};
                yRange = [1 max(k2Plot)+1];
                NP.TaskBaseClass.PlotEventWindows(ax, ce, 1, [], colNames, 'YRange', yRange);
                
                ax.YLim = yRange;
                ax.XLabel.String = "Time (s)";
                ax.YTick = 1 : max(k2Plot);
                ax.Title.String = "# clusters = " + k2Plot(i);
                ax.Box = 'off';
                
                if k2Plot(i) ~= sBatch.sEval.optiK
                    ax.Title.FontWeight = 'normal';
                end
            end
        end
        
        function PlotClusterMean(ce)
            % 
            %   PlotClusterMean(ce)
            % 
            
            uTb = ce.clusTb;
            cidList = unique(uTb.nmfcId);
%             cc = MPlot.Rainbow(numel(cidList));
            cc = zeros(numel(cidList), 3);
            
            ax = gca;
            hold(ax, 'on');
            
            colNames = {'cue1', 'stim', 'cue3', 'prod'};
            yRange = [-numel(cidList) 0];
            NP.TaskBaseClass.PlotEventWindows(ax, ce, 1, [], colNames, 'YRange', yRange);
            
            for i = 1 : numel(cidList)
                isClus = uTb.nmfcId == cidList(i);
                [t, y] = ce.GetArray('resp', [], find(isClus)+1, 'Normalization', 'minmaxsoft');
                [m, e] = MMath.MeanStats(y, 2);
                m = m-i;
                MPlot.ErrorShade(t, m, e, 'Color', cc(i,:), 'Parent', ax);
                plot(ax, t, m, 'Color', cc(i,:), 'LineWidth', 1.5);
            end
            
            ax.XLim = t([1 end]);
            ax.YLim = yRange;
            ax.XLabel.String = "Time (s)";
        end
        
        % Stability method
        function s = BootStability(X, varargin)
            % Bootstrap clustering to estimate the reliability of cluster membership
            % 
            %   s = BootClustering(X)
            %   s = BootClustering(..., 'nComp', 10)
            %   s = BootClustering(..., 'nBoot', 1e3)
            %   s = BootClustering(..., 'fraction', .9)
            % 
            % Inputs
            %   X               t-by-u array, where t is the number of samples, u is the number of units.
            %   'nComp'         The number of clusters to compute for.
            %   'nBoot'         The number of bootstrap iterations.
            %   'fraction'      The fraction of units to run bootstrap clustering. The rest are held out 
            %                   for testing.
            % Output
            %   Struct s with the following fields:
            %   s.nComp         The number of clusters computed, copied from 'nComp'.
            %   s.H             t-by-nComp array of template cluster centroid patterns (i.e. nonnegative 
            %                   left factor of X) computed from all units.
            %   s.W             nComp-by-u array of component weights (i.e. nonnegative right factor of X)
            %                   computed from all units.
            %   s.clustId       u-by-1 array of assgined cluster IDs.
            %   s.clustW        u-by-1 array of weights on the assigned clusters.
            %   s.bootH         t-by-nComp-by-nBoot array of bootstrap cluster centroid patterns.
            %   s.bootId        u-by-nBoot array of cluster membership IDs. Only the elements of testing 
            %                   units have values, others are NaN.
            %   s.idProb        u-by-nComp array
            %   s.maxProb       u-by-1 array
            %   s.bestId        u-by-1 array
            %   s.bestMeanW     u-by-1 array
            % 
            
            p = inputParser();
            p.addParameter('nComp', 10, @isscalar);
            p.addParameter('nBoot', 1e3, @isscalar);
            p.addParameter('fraction', .9, @isscalar);
            p.parse(varargin{:});
            nComp = p.Results.nComp;
            nBoot = p.Results.nBoot;
            frac = p.Results.fraction;
            
            [nBin, nUnit] = size(X);
            nSample = round(nUnit*frac);
            
            % Compute template
            rng(nComp);
            s = LMV.NMF.Clustering(X, nComp);
            H0 = s.H;
            
            % Bootstrap iterations
            rng(61);
            rsInd = cell(nBoot, 1);
            bootH = NaN(nBin, nComp, nBoot);
            for n = 1 : nBoot
                % Print progress
                str = sprintf('Iteration %i/%i\n', n, nBoot);
                fprintf(str);
                
                % Random sampling without replacement
                rng(nComp+n, 'twister');
                ind = randsample(nUnit, nSample, false);
                ind = MMath.Ind2Logical(ind, nUnit);
                
                % Factorization
                [~, H] = nnmf(X(:,ind)', nComp);
                
                rsInd{n} = ind;
                bootH(:,:,n) = H';
                
                % Clear printed
                if n < nBoot
                    fprintf(repmat('\b', 1, numel(str)));
                end
            end
            
            % 
            bootId = NaN(nUnit, nBoot);
            bootMaxW = NaN(nUnit, nBoot);
            for n = 1 : nBoot
                % Compute pairwise correlation matrix bootstrap and template patterns
                H = bootH(:,:,n);
                r = corr(H, H0);
                
                % Reorder bootstrap patterns based on best matches
                [~, rRank] = sort(r(:), 'descend'); % rank matches from best to worst
                I = zeros(nComp, 1);
                for m = 1 : numel(rRank)
                    % Find out which template (i) is matching which boot (j) pattern
                    [i, j] = ind2sub([nComp nComp], rRank(m));
                    
                    % Assign temp position index if not already
                    isAssigned = I(j);
                    isTempUsed = ismember(i, I); % each template can only be used once
                    if ~isAssigned && ~isTempUsed
                        I(j) = i;
                    end
                    
                    % End loop once all boot patterns are assigned with temp indices
                    if all(I)
                        break
                    end
                end
                H = H(:,I);
                
                % Compute cluster membership of the held-out units
                ind = rsInd{n};
                Wout = H \ X(:,~ind);
                [maxWout, I] = max(Wout);
                id = NaN(nUnit, 1);
                maxW = NaN(nUnit, 1);
                id(~ind) = I;
                maxW(~ind) = maxWout;
                
                bootH(:,:,n) = H;
                bootId(:,n) = id;
                bootMaxW(:,n) = maxW;
            end
            
            % Find most likely cluster membership
            for i = nUnit : -1 : 1
                % Compute probability distribution of cluster ID for each unit
                idProb(i,:) = histcounts(bootId(i,:), (1:nComp+1)-0.5);
                idProb(i,:) = idProb(i,:) ./ sum(idProb(i,:));
                
                % Find maximum likelihood and the corresponding cluster ID
                [maxProb(i,1), bestId(i,1)] = max(idProb(i,:));
                
                % Compute mean cluster weight of this cluster
                bestMeanW(i,1) = mean(bootMaxW(i, bootId(i,:) == bestId(i)));
                
                % Compute correlation between a PETH and its best factor
                [rho(i,1), pval(i,1)] = corr(s.H(:,bestId(i)), X(:,i));
            end
            
            s.bootH = bootH;
            s.bootId = bootId;
            s.idProb = idProb;
            s.maxProb = maxProb;
            s.bestId = bestId;
            s.bestMeanW = bestMeanW;
            s.rho = rho;
            s.pval = pval;
        end
        
        function [ce, sOpti] = FindOptimalK(ce, sBoot)
            % Compute the optimal number of components
            % 
            %   [ce, sOpti] = FindOptimalNComp(ce, sBoot)
            % 
            
            if isfield(sBoot(1), 'maxProb')
                % Compute probabilities above chance
                nCompList = [sBoot.nComp];
                P = cat(2, sBoot.maxProb);
                [Pmean, ~, ~, Pci] = MMath.MeanStats(P);
                dPmean = Pmean - (1./nCompList);
                dPci = Pci - (1./nCompList);
                
                % Choose the # of clusters
                [~, opIdx] = max(dPmean);
                opNum = nCompList(opIdx);
                
                % Add clustering info to the clusTb in ce
                ce.clusTb.nmfcId = sBoot(opIdx).bestId;
                ce.clusTb.nmfcProb = sBoot(opIdx).maxProb;
                ce.clusTb.nmfcW = sBoot(opIdx).bestMeanW;
                
                % Pack outputs
                sOpti.nCompList = nCompList;
                sOpti.Pmean = Pmean;
                sOpti.Pci = Pci;
                sOpti.dPmean = dPmean;
                sOpti.dPci = dPci;
                sOpti.optiIdx = opIdx;
                sOpti.optiNum = opNum;
            else
                % Compute average loss
                trainLoss = cat(1, sBoot.trainLoss);
                [mTrain, ~, ~, ciTrain] = MMath.MeanStats(trainLoss, 2);
                testLoss = cat(1, sBoot.testLoss);
                [mTest, ~, ~, ciTest] = MMath.MeanStats(testLoss, 2);
                
                % Find the # of clusters with lowest testing loss
                [~, opIdx] = min(mTrain);
                nCompList = [sBoot.nComp];
                opNum = nCompList(opIdx);
                
                % Add clustering info to the clusTb in ce
                ce.clusTb.nmfcId = sBoot(opIdx).clustId;
                ce.clusTb.nmfcW = mean(sBoot(opIdx).unitTestLoss, 2, 'omitnan');
                
                % Pack outputs
                sOpti.nCompList = nCompList;
                sOpti.meanTrain = mTrain;
                sOpti.ciTrain = ciTrain;
                sOpti.meanTest = mTest;
                sOpti.ciTest = ciTest;
                sOpti.optiIdx = opIdx;
                sOpti.optiNum = opNum;
            end
        end
        
        function PlotStability(sOpti, nComp2Plot)
            % The probability of a unit being grouped in the same cluster across bootstrap iterations
            
            % Optionally select a subset of components
            nCompList = sOpti.nCompList;
            if nargin < 2
                nComp2Plot = nCompList;
            end
            m = ismember(nCompList, nComp2Plot);
            
            % Mean±CI probability across different cluster #
            dPmean = sOpti.dPmean(m);
            dPci = sOpti.dPci(:,m);
            errorbar(nComp2Plot, dPmean, dPmean-dPci(1,:), dPci(2,:)-dPmean, 'Color', 'k'); hold on
            ax = gca;
            ax.XLim = [-1 1] + nComp2Plot([1 end]);
            ax.XTick = nComp2Plot;
            xlabel('# of clusters');
            ylabel('P(stable) - P(chance)');
            % title('Consistency of cluster membership');
            MPlot.Axes(gca);
        end
        
        function PlotStabilityCDF(sBatch)
            % CDFs of the cluster probability with different cluster #
            
            nCompList = sBatch.sOpti.nCompList;
            P = cat(2, sBatch.sBoot.maxProb);
            dPmean = sBatch.sOpti.dPmean;
            
            f = MPlot.Figure(46532); clf
            for i = 1 : numel(nCompList)
                ax = subplot(3,5,i);
                binEdges = 0:.02:1;
                N = histcounts(P(:,i), binEdges, 'Normalization', 'cdf');
                stairs(binEdges(1:end-1), N); hold on
                plot(1./nCompList([i i]), [0 1], '--', 'Color', [0 0 0]);
                xlim([0 1]);
                ylim([0 1]);
                xlabel("P(stable)");
                ylabel("Frac. of units");
                title(nCompList(i) + " comp, mean \DeltaP = " + round(dPmean(i),2));
                grid on
            end
            MPlot.Axes(gca);
            MPlot.Paperize(f, 'ColumnsWide', 2.5, 'ColumnsHigh', 1.2);
            MPlot.SavePNG(f, fullfile(anaDir, 'nmf_p_stable_cdf'));
        end
        
    end
    
end

