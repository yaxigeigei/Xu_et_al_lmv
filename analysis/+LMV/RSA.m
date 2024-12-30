classdef RSA
    
    properties(Constant)
        senMdl = "rsa_sen_id"; % sentence identity
        % senMdl = "rsa_sen_phc"; % sentence phoneme composition (histogram)
    end
    
    methods(Static)
        % Sentence representation similarity
        function mdls = Sentences(rTb, fTb, stimIds, phaseName)
            % Compute similarity between neuronal responses and sentence identities
            % 
            %   mdls = Sentences(rTb, fTb, stimIds, phaseName)
            % 
            
            % Get predictors
            X = rTb{:,:};
            
            % Get class labels
            Y = fTb{:,stimIds};
            assert(all(sum(Y,2)<=1), "Found more than one class label in the same observation(s)");
            Y = Y .* (1 : size(Y,2));
            y = sum(Y, 2);
            y = categorical(y, 1:numel(stimIds), stimIds); % i.e. indices of stim
            
            % Get observations from the task phase of interest
            M = fTb{:,phaseName};
            m = any(M,2) & ~isundefined(y);
            X = X(m,:);
            y = y(m);
            
            % Fit models
            fprintf("%i predictors\n", size(X,2));
            fprintf("%i raw samples, %i in task phase\n", numel(m), sum(m));
            fprintf("\nCompute RSA.\n");
            switch LMV.RSA.senMdl
                case "rsa_sen_id"
                    mdls = LMV.RSA.BatchRSA(X, y, 'DistFunY', @(x) squareform(1-(x==x')));
                case "rsa_sen_phc"
                    V = NP.Phone.PhonemeHistcounts(stimIds, fullfile(NP.Data.GetProjectRoot, "code", "tasks", "lmv", "TIMIT"));
                    Y = V(y,:);
                    mdls = LMV.RSA.BatchRSA(X, Y, 'DistFunY', @(x) pdist(x, "cityblock"));
                otherwise
                    error("'%s' is not a valid model name", LMV.RSA.senMdl);
            end
        end
        
        function mdls = BatchRSA(X, Y, varargin)
            % Run RSA between Y and each column of X
            % 
            %   mdls = BatchRSA(X, Y)
            %   mdls = BatchRSA(X, Y)
            %   mdls = BatchRSA(X, Y, ..., 'DistFunX', @(x) pdist(x, "cityblock"))
            %   mdls = BatchRSA(X, Y, ..., 'DistFunY', @pdist)
            %   mdls = BatchRSA(X, Y, ..., 'NPerm', 100)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('DistFunX', @(x) pdist(x, "cityblock"), @(x) isa(x, "function_handle"));
            p.addParameter('DistFunY', @pdist, @(x) isa(x, "function_handle"));
            p.addParameter('NPerm', 100, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Verbose', true, @islogical);
            p.parse(varargin{:});
            fx = p.Results.DistFunX;
            fy = p.Results.DistFunY;
            nPerm = p.Results.NPerm;
            isVerbose = p.Results.Verbose;
            funcArgs = [fieldnames(p.Unmatched) struct2cell(p.Unmatched)]';
            
            % Compute feature RDM
            Dy = fy(Y);
            [C, ~, ic] = unique(Y, 'rows');
            
            nMdls = size(X, 2);
            mdls = cell(nMdls, 1);
            for i = 1 : nMdls
                if isVerbose
                    str = sprintf('Fitting models %i/%i\n', i, nMdls);
                    fprintf(str);
                end
                
                % Prepare input for this iteration
                x = X(:,i);
                
                % Remove outliers in each group
                for k = 1 : numel(C)
                    m = ic == k;
                    xClean = x(m);
                    xClean(isoutlier(xClean)) = NaN;
                    x(m) = xClean;
                end
                
                % RSA
                % Mx = abs(x-x');
                % Mx(logical(eye(size(Mx)))) = 0;
                % Dx = squareform(Mx)';
                Dx = fx(x);
                r = corr(Dx(:), Dy(:), 'Rows', 'complete');
                
                s = struct;
                s.x = x;
                s.y = ic;
                s.r = r;
                mdls{i} = s;
                
                % Fit null models with permutation
                if ~nPerm
                    continue
                end
                
                rNull = zeros(1, numel(x));
                for n = 1 : nPerm
                    I = randperm(numel(x));
                    % Mx = abs(x(I)-x(I)');
                    % Mx(logical(eye(size(Mx)))) = 0;
                    % DxNull = squareform(Mx)';
                    DxNull = fx(x(I));
                    rNull(n) = corr(DxNull(:), Dy(:), 'Rows', 'complete');
                end
                
                s.nPerm = nPerm;
                s.rNull = rNull;
                s.pval = MMath.EstimatePval(r, rNull, 'Tail', 'right');
                mdls{i} = s;
            end
        end
        
        % Speech feature representation
        function ops = GetOptions(ops)
            % 
            if nargin < 1
                ops = struct;
            end
            ops.stimIdList = LMV.Param.stimIdList14;
        end
        
        function s = Run(ce, ops)
            % 
            s = NP.RSA.PrepareInput(ce, ops);
            s = NP.RSA.ComputeRDMs(s);
            s = NP.RSA.ComputeSim(s);
        end
        
        function s = PrepareEventInputs(ce, ops)
            % Prepare input patterns for computing RDMs
            % 
            %   s = PrepareEventInputs(ce, ops)
            % 
            % Inputs
            %   ce                  
            %   ops.unitInd
            %   ops.repFeats
            % Outputs
            %   s.ce
            %   s.t
            %   s.stimIdx
            %   s.repNames
            %   s.input
            % 
            
            % Select units
            ce = ce.Duplicate;
            unitMask = true(size(ce.respNames));
            unitMask(ops.unitInd) = false;
            ce.RemoveUnits(unitMask);
            
            % Get arrays
            [t, R] = ce.GetArray('resp', 'Normalization', 'zscore');
            repFeats = struct2cell(ops.repFeats)';
            [~, FF] = cellfun(@(x) ce.GetArray('feat', [], x, 'Normalization', 'zscore'), repFeats, 'Uni', false);
            XX = [{R}, FF];
            
            % Mask samples
            lb = t.GetParentLabel;
            m = ismember(lb, ops.eventMask);
            t = t(m);
            lb = lb(m);
            XX = cellfun(@(x) x(m,:), XX, 'Uni', false);
            
            % Compute group mean by event labels
            [G, lb] = findgroups(lb);
            for i = 1 : numel(XX)
                X = XX{i};
                X(isnan(X)) = 0;
                X = splitapply(@(x) mean(x,1), X, G);
                XX{i} = X;
            end
            
            % Pack results
            s = ops;
            s.ce = ce;
            s.t = double(t);
            s.eventNames = lb;
            s.repNames = ["Neural", fieldnames(ops.repFeats)'];
            s.input = XX;
        end
        
        function s = PrepareTimeseriesInput(ce, ops)
            % Prepare input patterns for computing RDMs
            % 
            %   s = PrepareTimeseriesInput(ce, ops)
            % 
            % Inputs
            %   ce                  
            %   ops.unitInd
            %   ops.dtNeural
            %   ops.featFam
            %   ops.maskName
            %   ops.downsample
            % Outputs
            %   s.ce
            %   s.t
            %   s.stimIdx
            %   s.repNames
            %   s.input
            % 
            
            % Select units
            ce = ce.Duplicate;
            unitMask = true(size(ce.respNames));
            unitMask(ops.unitInd) = false;
            ce.RemoveUnits(unitMask);
            
            % Get arrays
            [t, stimIdx] = ce.GetArray('feat', [], 'stimNumId');
            
            [~, R] = ce.GetArray('resp', 'TimeShifts', ops.dtNeural, 'Normalization', 'zscore');
            
            featFam = struct2cell(ops.featFam)';
            [~, FF] = cellfun(@(x) ce.GetArray('feat', [], x, 'Normalization', 'zscore'), featFam, 'Uni', false);
            
            XX = [{R}, FF];
            
            % Mask samples
            [~, M] = ce.GetArray('feat', [], ops.maskName);
            m = all(logical(M), 2);
            t = t(m);
            stimIdx = stimIdx(m);
            XX = cellfun(@(x) x(m,:), XX, 'Uni', false);
            
            % Downsampling
            if ops.downsample > 1
                t = downsample(t, ops.downsample);
                stimIdx = downsample(stimIdx, ops.downsample);
                XX = cellfun(@(x) downsample(x, ops.downsample), XX, 'Uni', false);
            end
            
            % Pack results
            s = ops;
            s.ce = ce;
%             s.recMeta = ce.userData.recMeta;
%             s.clusTb = ce.clusTb(unitInd,:);
            s.t = t;
            s.stimIdx = stimIdx;
            s.repNames = ["Neural", fieldnames(ops.featFam)'];
            s.input = XX;
        end
        
        function s = ComputeRDMs(s)
            % Compute RDMs
            % 
            %   s = ComputeRDMs(s)
            % 
            % Inputs
            %   s.input
            % Outputs
            %   s.RDM
            %   s.RDV
            
            XX = s.input;
            
            [RDM, RDV] = NP.RSA.IComputeRDMs(XX);
            
            s.RDM = RDM;
            s.RDV = RDV;
        end
        
        function s = ComputeSim(s)
            % Compute similarity between RDMs
            % 
            %   s = ComputeSim(s)
            % 
            % Inputs
            %   s.input
            % Outputs
            %   s.RDM
            %   s.RDV
            
            % Compute pairwise similarity
            [PSM, PSV] = NP.RSA.IComputePSM(s.RDV);
            
            s.PSM = PSM;
            s.PSV = PSV;
            
            % Bootstrap null similarity and pval
            nBoot = s.nBoot;
            if nBoot == 0
                s.pval = zeros(size(PSM));
                return
            end
            nSp = numel(s.t);
            nRDM = size(PSM, 1);
            bPSV = zeros(numel(PSV), nBoot);
            rng(61);
            for n = 1 : nBoot
                nShift = randi([1 nSp], [1 nRDM]);
                shRDM = cellfun(@(x,k) circshift(circshift(x,k,1),k,2), s.RDM, num2cell(nShift), 'Uni', false);
                [~, bPSV(:,n)] = NP.RSA.IComputePSM(shRDM);
            end
            PV = MMath.EstimatePval(PSV, bPSV, 'Method', 'Gumbel', 'Tail', 'left');
            PM = squareform(PV);
            
            s.pval = PM;
        end
        
        function s = LightenResult(s)
            % 
            fld2rm = {'ce', 'input', 't', 'stimIdx', 'RDM', 'RDV'};
            fld2rm = intersect(fld2rm, fieldnames(s));
            s = rmfield(s, fld2rm);
        end
        
        function s = RunSelection(s)
            % Optimize feature and unit selection via monte-carlo sampling and credit assignment
            % 
            %   s = RunSelection(s)
            % 
            % Inputs
            %   s.selection.nIter
            %   s.selection.nFrac
            %   s.selection.fFrac
            % Outputs
            %   s.selection.isUnit
            %   s.selection.isFeat
            %   s.selection.sim
            % 
            
            if ~isfield(s, 'selection')
                return
            end
            sel = s.selection;
            
            XX = s.input;
            for i = 1 : numel(XX)
                XX{i}(isnan(XX{i})) = 0;
            end
            Xn = XX{1};
            Xf = XX{2};
                
            nu = size(Xn,2);
            nf = size(Xf,2);
            nuRS = round(nu * sel.uFrac);
            nfRS = round(nf * sel.fFrac);
            
            % Preallocate outputs
            nIter = sel.nIter;
            isU = false(nu, nIter);
            isF = false(nf, nIter);
            sim = NaN(1,nIter);
            
            % Iterate RSA
            parfor i = 1 : nIter
                mU = false(nu,1);
                mU(randperm(nu, nuRS)) = true;
                isU(:,i) = mU;
                
                mF = false(nf,1);
                mF(randperm(nf, nfRS)) = true;
                isF(:,i) = mF;
                
                Du = pdist(Xn(:,mU), 'seuclidean');
                Df = pdist(Xf(:,mF), 'seuclidean');
                sim(i) = corr(Du', Df');
            end
            
            % Compute mean similarity as credits
            sU = sim .* isU;
            sU = sum(sU,2) ./ sum(isU,2);
            
            sF = sim .* isF;
            sF = sum(sF,2) ./ sum(isF,2);
            
            % Pack outputs
%             sel.isUnit = isU;
%             sel.isFeat = isF;
%             sel.sim = sim;
            sel.simUnit = sU;
            sel.simFeat = sF;
            s.selection = sel;
        end
        
        function [MM, VV] = IComputeRDMs(XX)
            % 
            
            if ~iscell(XX)
                XX = {XX};
            end
            
            for i = numel(XX) : -1 : 1
                X = XX{i};
                X(isnan(X)) = 0;
                
%                 D = pdist(X, 'mahalanobis');
                D = pdist(X, 'correlation');
%                 D = pdist(X, 'seuclidean');
                M = squareform(D);
                
                MM{i} = M;
                VV{i} = D;
            end
        end
        
        function r = IComputeSim(A, B)
            % Compute similarity score between A and each RDM in B
            % 
            %   r = IComputeSim(A, B)
            % 
            
            % Standardize distances to vector form and concatenate them into a single matrix
            if iscell(A)
                A = A{1};
            end
            if ismatrix(A)
                A = squareform(A, 'tovector');
            end
            
            if ndims(B) == 3
                B = mat2cell(B, [], [], ones(size(B,3),1)); % Seperate matrices into a cell array
            end
            if iscell(B)
                for i = 1 : numel(B)
                    if ~isvector(B{i})
                        B{i} = squareform(B{i}, 'tovector');
                    end
                    B{i} = B{i}(:);
                end
                B = cat(2, B{:});
            end
            
            % Compute pairwise similarity
            r = corr(A, B, 'rows','complete');
        end
        
        function [M, V] = IComputePSM(DD)
            % 
            % 
            %   [M, V] = IComputePSM(DD)
            % 
            
            % Seperate matrices into a cell array
            if ndims(DD) == 3
                DD = mat2cell(DD, [], [], ones(size(DD,3),1));
            end
            
            % Standardize distances to vector form and concatenate them into a single matrix
            if iscell(DD)
                for i = 1 : numel(DD)
                    if ~isvector(DD{i})
                        DD{i} = squareform(DD{i}, 'tovector');
                    end
                    DD{i} = DD{i}(:);
                end
                DD = cat(2, DD{:});
            end
            
            % Compute pairwise similarity
            V = 1 - pdist(DD', 'correlation');
            M = squareform(V);
        end
        
        function PlotRDM(RDM, G)
            % 
            
            if nargin > 1 && ~isempty(G)
                % Sort input based on ordinal information in G
                if iscategorical(G) && isordinal(G)
                    [G, I] = sort(G);
                    RDM = RDM(I,I);
                    gNames = unique(G);
                else
                    G = categorical(G);
                    gNames = unique(G, 'stable');
                end
            end
            
            % Plot matrix
            imagesc(RDM);
            colormap(flip(brewermap([], 'RdBu')));
            
            ax = gca;
            ax.Title.Interpreter = 'none';
            MPlot.Axes(ax);
            axis xy equal tight off
            
            if nargin < 2 || isempty(G)
                return
            end
            
            n = size(RDM,1);
            if numel(G) > numel(gNames)
                % Plot ribbons
                ribPos = ([-0.02 0] - 0.005) * n;
                txtPos = -0.01*n;
                
                [xc, yc] = MPlot.GroupRibbon(G, ribPos);
                text(xc, yc+txtPos, gNames, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8);
                
                [xc, yc] = MPlot.GroupRibbon(ribPos, G);
                text(xc+txtPos, yc, gNames, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 8);
                
            elseif numel(G) == numel(gNames)
                % Label ticks
                ax.XTick = 1 : n;
                ax.YTick = 1 : n;
                ax.XTickLabel = gNames;
                ax.YTickLabel = gNames;
                ax.XTickLabelRotation = 0;
                axis on
            end
        end
        
        function PlotPSM(s)
            % 
            
            C = s.PSM;
            P = s.pval;
            repNames = s.repNames;
            n = numel(repNames);
            C(logical(eye(n))) = NaN;
            C(P > 0.05) = NaN;
            
            h = heatmap(repNames, repNames, C);
            h.CellLabelFormat = '%0.2f';
        end
        
        function PlotSimTimecourse(sArray, recIds)
            % 
            
            SS4 = cellfun(@(x) permute(x.PSM(2:end,1), [2 3 4 1]), sArray, 'Uni', false);
            S = cell2mat(SS4);
            
            [~, nPhase, ~] = size(sArray);
            repNames = sArray{1}.repNames(2:end);
            nFam = numel(repNames);
            
            % Set up figure
            tl = tiledlayout(nPhase, nFam);
            tl.Padding = 'compact';
            for i = 1 : nPhase
                % Get task phase and time offsets
                mn = sArray{1,i,1}.maskName;
                dt = cellfun(@(x) x.dtNeural, sArray(1,i,:));
                dt = squeeze(dt);
                
                % Plot through representations
                for j = 1 : nFam
                    ax = nexttile;
                    ss = squeeze(S(:,i,:,j))';
                    hh = plot(dt, ss, 'Color', [0 0 0 0.2]); hold on
                    plot(dt, mean(ss,2), 'LineWidth', 2);
                    
                    for k = 1 : numel(hh)
                        hh(k).DataTipTemplate.Interpreter = 'none';
                        hh(k).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("recId", repelem(recIds(k), size(ss,1)));
                    end
                    
                    ax.XLim = dt([1 end]);
                    ax.YLim = [0 .2];
                    ax.XLabel.String = 'Neural \Deltat (s)';
                    ax.YLabel.String = 'Similarity';
                    ax.Title.String = mn + ": " + repNames(j);
                    MPlot.Axes(ax);
                end
            end
        end
        
        function PlotSimOptimal(sArray, recTb, regions)
            % Plot sets of bars showing similarity where each set is a representation and each bar in a set is a brain region
            % 
            %   PlotSimOptimal(sArray, recTb)
            % 
            % Inputs
            %   sArray          rec-by-dt cell array.
            %   recTb           ...
            % 
            
            repNames = sArray{1}.repNames(2:end);
            maskName = sArray{1}.maskName;
            
            % Get similarity array
            sArray = squeeze(sArray);
            s = cellfun(@(x) permute(x.PSM(2:end,1), [2 3 1]), sArray, 'Uni', false);
            s = cell2mat(s); % rec x dt x fam
            s = permute(s, [1 3 2]); % rec x fam x dt
            
            % Find optimal similarity values
            if nargin < 3
                regions = unique(recTb.region);
            end
            regCell = cell(size(regions));
            for i = 1 : numel(regions)
                % Get recordings within a region
                m = recTb.region==regions(i);
                sReg = s(m,:,:); % rec x fam x dt
                
                % Find the maximal time point from the mean for each family
                sRegMean = squeeze(mean(sReg, 1)); % fam x dt
                [~, I] = max(sRegMean, [], 2);
                if maskName == "stim"
                    I = ones(size(I))*16;
                else
                    I = ones(size(I))*11;
                end
                
                % Take the values at optimal timepoints for each family
                sRegOpti = zeros(size(sReg,1), size(sReg,2)); % rec x fam
                for j = 1 : numel(I)
                    sRegOpti(:,j) = sReg(:,j,I(j));
                end
                regCell{i} = sRegOpti;
                
                % Compute mean stats
                regMean = MMath.MeanStats(sRegOpti);
                
                % 
                ax = nexttile;
                h = bar(ax, regMean, 'FaceColor', 'none');
                hold(ax, 'on');
                plot(ax, h.XData', sRegOpti', 'x-', 'Color', [0 0 0 .15]);
                ax.XTickLabel = repNames;
                ax.YLim = [0 .15];
                ax.YLabel.String = 'Similarity';
                ax.Title.String = regions(i) + " " + maskName;
                MPlot.Axes(ax);
            end
            
%             % Compute mean stats
%             regMean = cellfun(@(x) MMath.MeanStats(x), regCell, 'Uni', false);
%             regMean = cell2mat(regMean); % reg x fam
%             
%             ax = nexttile;
%             h = bar(ax, regMean, 'EdgeColor', 'flat', 'FaceColor', 'none');
%             ax.XTickLabel = regions;
%             ax.YLim = [0 .1];
%             ax.YLabel.String = 'Similarity';
%             MPlot.Axes(ax);
            
%             [~, nPhase, ~] = size(sArray);
%             repNames = sArray{1}.repNames(2:end);
%             nFam = numel(repNames);
%             
%             % Plot through representations
%             for j = 1 : nFam
%                 ax = nexttile;
%                 ss = squeeze(s(:,i,:,j))';
%                 hh = plot(dt, ss, 'Color', [0 0 0 0.2]); hold on
%                 plot(dt, mean(ss,2), 'LineWidth', 2);
%                 
%                 for k = 1 : numel(hh)
%                     hh(k).DataTipTemplate.Interpreter = 'none';
%                     hh(k).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("recId", repelem(recIds(k), size(ss,1)));
%                 end
%                 
%                 ax.XLim = dt([1 end]);
%                 ax.YLim = [0 .2];
%                 ax.XLabel.String = 'Neural \Deltat (s)';
%                 ax.YLabel.String = 'Similarity';
%                 ax.Title.String = mn + ": " + repNames(j);
%                 MPlot.Axes(ax);
%             end
        end
        
        function PlotInputFeatures(sRez)
            % 
            
            s = sRez.stimIdx;
            XX = sRez.input(2:end);
            repNames = sRez.repNames(2:end);
            nRep = numel(XX);
            
            sp = cell(size(XX));
            yL = cell(size(XX));
            for i = 1 : nRep
                X = XX{i};
                
                % Get labels
                if ismember(repNames(i), fieldnames(sRez.featFam))
                    yL{i} = sRez.featFam.(repNames(i));
                else
                    yL{i} = cell(1, size(X,2));
                end
                if repNames(i) == "Phone"
                    isRm = ~any(X>0, 1);
                    X(:,isRm) = [];
                    yL{i} = MLing.ARPA2IPA(yL{i}(~isRm));
                end
                
                % Add sentence indicator row
                if i == 1
                    X = [s, X];
                    yL{i} = [{''} yL{i}];
                end
                
                % Normalize between 0 and 1
                X = MMath.Normalize(X, 'minmax');
                
                % Adjust spacing
                sp{i} = ones(size(X,2), 1);
                if repNames(i) ~= "Mel"
                    X = repelem(X, 1, 8);
                    sp{i} = sp{i} * 8;
                else
                    sp{i} = sum(sp{i})/2;
                end
                
                XX{i} = X;
            end
            Xcat = cat(2, XX{:});
            ycat = cumsum(cat(1, sp{:})) - 4;
            yLcat = cat(2, yL{:});
            
            % Plot features heatmap
            t = (1:size(Xcat,1)) * sRez.rsBinSize*sRez.downsample;
            CData = 1 - Xcat';
            CData(isnan(CData)) = 1;
            imagesc(t, [], CData); hold on
            colormap gray
            
            % Plot phonetic labels
            tt = sRez.ce.GetTable('taskTime');
            tge = tt.(sRez.maskName);
            if iscell(tge)
                tge = cat(1, tge{:});
            end
            iSen0 = [1; find(diff(s))+1];
            tge = tge - double(tge) + (t(iSen0))';
%             tge = tge - double(tge) + cumsum([sRez.rsBinSize*sRez.downsample/2; tge(1:end-1).GetDuration]);
            tge.PlotTiers(-40, 10);
            
            ax = gca;
            ax.YLim(1) = -30;
            ax.YTick = ycat;
            ax.YTickLabel = yLcat;
            ax.Title.String = NP.SE.GetID(sRez.ce);
            ax.Title.Interpreter = 'none';
            MPlot.Axes(ax);
        end
        
        function PlotCredits(sRez)
            % 
            
            s = sRez.selection;
            sU = s.simUnit;
            sF = s.simFeat;
            mU = isoutlier(sU);
            mF = isoutlier(sF);
            
            ax = nexttile;
            stem(ax, sU); hold on
            stem(ax, find(mU), sU(mU), 'filled');
            ax.XLim = [0 numel(sU)+1];
            ax.XTick = 1 : numel(sU);
            ax.XTickLabel = "u"+sRez.ce.clusTb.clusId;
            ax.YLabel.String = "Similarity (r)";
            MPlot.Axes(ax);
            view(90,90);
            
            ax = nexttile;
            stem(ax, sF); hold on
            stem(ax, find(mF), sF(mF), 'filled');
            ax.XLim = [0 numel(sF)+1];
            ax.XTick = 1 : numel(sF);
            ax.XTickLabel = sRez.featFam.(sRez.repNames(2));
            ax.YLabel.String = "Similarity (r)";
            MPlot.Axes(ax);
            view(90,90);
            
        end
        
    end
    
end

