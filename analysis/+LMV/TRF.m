classdef TRF < NP.Fit
    
    properties(Constant)
        featSets = ["phone", "strf", "artic3"];
        targets = ["stim", "feedback", "prod"];
    end
    
    methods(Static)
        % Model fitting
        function feats = GetFeatureSet(featSet, target)
            % Get the names of member features in a given feature set
            % 
            %   feats = GetFeatureSet(featSet)
            %   feats = GetFeatureSet(featSet, target)
            % 
            
            switch featSet
                case 'phone'
                    feats = [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.consonants];
                    feats = setdiff(feats, {'AX' 'AXR' 'IX' 'OY' 'SH' 'ZH'});
                case 'strf'
                    if target == "stim"
                        feats = {'speaker1'};
                    else
                        feats = {'mic'};
                    end
                case 'artic3'
                    feats = LMV.TRF.GetFeatureSet('artic');
                    feats = [feats, "rF0", "d_"+feats, "drF0"];
                case 'artic2'
                    feats = [LMV.TRF.GetFeatureSet('artic'), "rF0", "dEnv"];
                case 'artic'
                    feats = ["ja", "la", "pro", "ttcc", "tbcc", "tdcc", "ttcl", "tbcl", "v_x", "v_y"];
                case 'acous'
                    feats = [{'env', 'peakEnv', 'peakRate'}, {'brF0', 'drF0'}];
                case 'brF0'
                    feats = {'brF0'};
                case 'fujisaki'
                    feats = {'voicing', 'phrase', 'accent'};
                    
                case 'pm'
                    feats = {'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', ...
                        'nasals', 'approximants', 'labial', 'velar', 'coronal', 'glottal', 'dental'};
                    
                case 'combo1'
                    feats = [{'sentOn'}, {'peakEnv', 'peakRate'}, {'rF0', 'drF0'}, ...
                        {'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', ...
                        'nasals', 'approximants', 'labial', 'velar', 'coronal', 'glottal', 'dental'}];
                case 'combo2'
                    feats = [{'sentOn'}, {'peakEnv', 'peakRate'}, {'brF0', 'drF0'}, ...
                        {'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', ...
                        'nasals', 'approximants', 'labial', 'velar', 'coronal', 'glottal', 'dental'}];
                case 'combo3pm'
                    feats = [{'peakEnv', 'peakRate'}, {'voicing', 'phrase', 'accent', 'drF0'}, ...
                        {'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', ...
                        'nasals', 'approximants', 'labial', 'velar', 'coronal', 'glottal', 'dental'}, ...
                        {'sentOn', 'wordOn', 'syllOn', 'nSyll', 'blick', 'wordFreq'}];
                case 'combo3akt'
                    feats = [{'peakEnv', 'peakRate'}, {'voicing', 'phrase', 'accent', 'drF0'}, ...
                        {'ja', 'la', 'pro', 'ttcc', 'tbcc', 'tdcc', 'ttcl', 'tbcl', 'v_x', 'v_y'}, ...
                        {'sentOn', 'wordOn', 'syllOn', 'nSyll', 'blick', 'wordFreq'}];
                case 'combo4pm'
                    feats = [{'peakEnv', 'peakRate'}, {'voicing', 'phrase', 'accent', 'drF0'}, ...
                        {'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', ...
                        'nasals', 'approximants', 'labial', 'velar', 'coronal', 'dental'}, ...
                        {'sentOn', 'wordOn', 'syllOn', 'nSyll', 'blick', 'wordFreq'}];
                    
                otherwise
                    error("'%s' is not a defined feature set.", featSet);
            end
        end
        
        function [tShifts, phaseName] = GetTargetParams(targetName, tShort, tLong, stepSize)
            % Get the parameters for time shifts and the task phase of interest for a given target
            % 
            %   [dt, phaseName] = GetTargetParams(targetName, tShort, tLong, stepSize)
            % 
            if nargin < 2
                tShort = 0.1;
            end
            if nargin < 3
                tLong = 0.4;
            end
            if nargin < 4
                stepSize = 0.02;
            end
            switch targetName
                case 'stim'
                    tShifts = -tShort : stepSize : tLong;
                    phaseName = 'stim';
                case 'prod'
                    tShifts = -tLong : stepSize : tShort;
                    phaseName = 'prod';
                case 'feedback'
                    tShifts = -tShort : stepSize : tLong;
                    phaseName = 'prod';
                otherwise
                    error("'%s' is not a supported target name.", targetName);
            end
        end
        
        function clusTb = FitSpeechEncoding(ce, uInd, varargin)
            % Fit speech feature encoding models
            % 
            %   clusTb = FitSpeechEncoding(ce, uInd, 'FeatSet', '', 'Target', '')
            %   clusTb = FitSpeechEncoding(ce, [], 'FeatSet', '', 'Target', '')
            %   clusTb = FitSpeechEncoding(..., 'Lambda', [])
            %   clusTb = FitSpeechEncoding(..., 'NBoot', [])
            % 
            % Inputs
            %   ce          An NP.CodingExplorer object including a clusTb in its metadata.
            %   uInd        The indices of units to fit. Empty value [] includes all units.
            %   'FeatSet'   Names of a set of features.
            %   'Target'    Name of the target time period.
            %   'Lambda'    If empty [], will search the parameter space by fitting models with a series of lambda values.
            %               If provided, this lambda value will be used to fit all models, and models with randomly shifted 
            %               data to bootstrap the null distribution of model weights.
            % Output
            %   clusTb      
            % 
            
            p = inputParser;
            p.addParameter('FeatSet', '', @(x) isstring(x) || ischar(x));
            p.addParameter('Target', '', @(x) isstring(x) || ischar(x));
            p.addParameter('Lambda', [], @isnumeric);
            p.addParameter('NBoot', 100, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Alpha', eps, @isnumeric);
            p.parse(varargin{:});
            featSet = p.Results.FeatSet;
            target = p.Results.Target;
            L = p.Results.Lambda;
            nShifts = p.Results.NBoot;
            alpha = p.Results.Alpha;
            
            % Configure model fitting
            if isempty(uInd)
                uInd = (1 : ce.numResp)';
            end
            feats = LMV.TRF.GetFeatureSet(featSet, target);
            [dt, poi] = LMV.TRF.GetTargetParams(target, 0.1, 0.4, 0.02);
            
            % Get predictors
            [~, X, feats] = ce.GetArray('feat', [], feats, 'TimeShifts', dt);
            X(isnan(X)) = 0;
            
            % Get responses
            [~, Y, resps] = ce.GetArray('resp', [], uInd+1);
            Y(isnan(Y)) = 0;
            
            % Masking
            masks = {'cue1', 'atten', 'stim', 'delay', 'cue3', 'init', 'prod', 'iti'};
            [~, M] = ce.GetArray('feat', [], masks);
            maskTb = array2table(logical(M), 'VariableNames', masks);
            m = any(maskTb{:,poi}, 2);
            X = X(m,:);
            Y = Y(m,:);
            
            % Normalization
            X = MMath.Normalize(X, 'zscore');
            Y = MMath.Normalize(Y, 'zscore');
            X(isnan(X)) = 0;
            Y(isnan(Y)) = 0;
            
            % Print input info
            fprintf("\n%s\n", NP.SE.GetID(ce));
            fprintf("%i predictors, %i responses\n", size(X,2), size(Y,2));
            fprintf("%i raw samples, %i masked, %i valid\n", numel(m), sum(m), sum(~any(isnan([X Y]) | isinf([X Y]), 2)));
            
            % Get clusTb
            NP.Unit.MergeMeta(ce);
            NP.Unit.SetUniqueClusId(ce);
            clusTb = NP.Unit.GetClusTb(ce);
            clusTb = NP.Unit.AddRecMeta(ce, clusTb);
            
            % Model fitting
            mdlName = featSet+"_"+target;
            ss = statset('UseParallel', true);
            if isempty(L)
                % Search lambda
                fprintf("\nSearch lambda for %s encoding models.\n", mdlName);
                lambdaList = logspace(-3, 3, 30);
                mdls = NP.Fit.BatchLasso(X, Y, ...
                    'Alpha', alpha, 'Standardize', false, 'Lambda', lambdaList, 'Options', ss); % {'Alpha', eps} is ridge regressions
                colName = mdlName+"_HS";
            else
                % Fit model with optimal hyperparam
                fprintf("\nFit %s encoding models with optimal hyperparameters.\n", mdlName);
                if size(X,2) > 1000
                    fprintf("Use reduced bootstraping null distributions since the model is too large.\n");
                    nShifts = min(nShifts, 50);
                end
                mdls = LMV.TRF.BatchLasso(X, Y, ...
                    'Alpha', alpha, 'Standardize', false, 'Lambda', L, 'Options', ss, ...
                    'NumShifts', nShifts, 'MinShift', round(0.1*sum(m)));
                colName = mdlName;
            end
            
            % Add additional info
            for i = 1 : numel(mdls)
                mdls{i}.name = mdlName;
                mdls{i}.feats = feats;
                mdls{i}.resps = resps(i);
                mdls{i}.dt = dt;
            end
            
            clusTb.(colName)(uInd) = mdls;
        end
        
        % Utilities
        function [r, p] = GetModelR(mdls)
            % Get r-values from a cell array of model structs
            % 
            %   [r, p] = GetModelR(mdls)
            % 
            r2 = NaN(size(mdls));
            p = NaN(size(mdls));
            for i = 1 : size(mdls,1)
                for j = 1 : size(mdls,2)
                    m = mdls{i,j};
                    if ~isfield(m, 'r2')
                        continue
                    end
                    r2(i,j) = max(m.r2, 0); % avoid negative cv r2
                    if isfield(m, 'null')
                        p(i,j) = m.null.r2Pval;
                    end
                end
            end
            r = sqrt(r2);
        end
        
        function cTbs = LoadModels(mdlNames, varargin)
            % Load model files with tables of the same recording merged into one
            % 
            %   mTbs = LoadModels(mdlNames)
            %   mTbs = LoadModels(mdlNames, recIds)
            %   mTbs = LoadModels(mdlNames, ..., 'LoadTo', []);
            % 
            % Inputs
            %   mdlNames        A list of model names with the format "<feature_set>_<target>".
            %   recIds          Include only recordings based on this list of recording IDs.
            %   'LoadTo'        A cell array of clusTb to load models into.
            % Output
            %   uTbs            A cell array of clusTb with models added to the last columns.
            % 
            
            p = inputParser;
            p.addOptional('recIds', [], @(x) ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('LoadTo', [], @(x) iscell(x) || istable(x));
            p.parse(varargin{:});
            recIds = string(p.Results.recIds);
            cTbs = p.Results.LoadTo;
            
            mdlNames = string(mdlNames);
            
            for i = 1 : numel(mdlNames)
                % Get model directory
                mn = mdlNames(i);
                mnParts = split(mn, "_");
                setName = mnParts{1};
                mdlDir = fullfile(LMV.Data.GetAnalysisDir, "coding", "stRF", setName, "mdls");
                
                % Find recordings
                if isempty(recIds)
                    mSearch = MBrowse.Dir2Table(fullfile(mdlDir, "*_clusTb_" + mn + ".mat"));
                    recIds = unique(cellfun(@(x) string(NP.SE.GetID(x)), mSearch.name));
                end
                recIds = string(recIds);
                
                % Initialize output tables
                if isempty(cTbs)
                    cTbs = repmat({table}, size(recIds));
                end
                if ~iscell(cTbs)
                    cTbs = {cTbs};
                end
                
                % Load model tables
                for j = 1 : numel(recIds)
                    cTb = cTbs{j};
                    load(fullfile(mdlDir, recIds(j)+"_clusTb_"+mn+".mat"), 'clusTb');
                    m = ~ismember(clusTb.Properties.VariableNames, cTb.Properties.VariableNames);
                    if isempty(cTb)
                        cTb = [cTb clusTb(:,m)];
                    else
                        [~, I] = MMath.SortLike(clusTb.clusId, cTb.clusId);
                        cTb = [cTb clusTb(I,m)];
                    end
                    cTbs{j} = cTb;
                end
            end
            
            % Return table if a single recording is requested
            if isscalar(recIds)
                cTbs = cTbs{1};
            end
        end
        
        % Clustering weights
        function s = ClusterModels(mdls)
            % Hierarchical clustering of TRFs
            % 
            %   s = ClusterModels(mdls)
            % 
            % Input
            %   mdls            A cell array of model structs.
            % Outputs
            %   s.mdlNames      Names of the models.
            %   s.feats         Feature names.
            %   s.resps         Unit clusIds as response names.
            %   s.maxBeta       A units-by-feats array of maximal coefficient for each feature and each unit.
            %   s.D             Vector form of pairwise distances.
            %   s.Z             A matrix that encodes a tree containing hierarchical clusters.
            %   s.cc            Cophenetic correlation coefficient.
            %   s.nc            Number of clusters, hardcoded to be 10.
            %   s.hcId          Cluster ID of hierarchical clustering.
            % 
            
            nc = 5; % number of clusters
            
            hasMdl = ~cellfun(@isempty, mdls);
            mdls = mdls(hasMdl);
            
            maxBeta = cell(size(mdls));
            for i = 1 : numel(mdls)
                % Mask insignificant weights
                mdl = mdls{i};
                B = reshape(mdl.Beta, numel(mdl.feats), numel(mdl.dt));
                if isfield(mdl, 'null')
                    B(mdl.null.BetaPval > 0.05) = 0;
                end
                
                % Trim acausal and distant parts of the weight matrix
                if contains(mdl.name, 'prod')
                    mTrim = mdl.dt < 0 & mdl.dt > -0.3;
                else
                    mTrim = mdl.dt > 0 & mdl.dt < 0.3;
                end
                B = B(:, mTrim);
                
                % Take the strongest weight from each feature
                [~, I] = max(abs(B), [], 2);
                ind = sub2ind(size(B), (1:size(B,1))', I);
                maxBeta{i} = B(ind);
            end
            maxBeta = cat(2, maxBeta{:})'; % units-by-feats after transpose
            
            s.hasMdl = hasMdl;
            s.mdlNames = cellfun(@(x) x.name, mdls);
            s.resps = cellfun(@(x) x.resps, mdls);
            s.maxBeta = maxBeta;
            if numel(mdls) < nc
                disp("Insufficient number of models to run the clustering.");
                s.feats = [];
                s.D = [];
                s.Z = [];
                s.cc = [];
                s.nc = [];
                s.hcId = [];
                s.umapCoords = [];
                return
            end
            
            % % Normalize weights
            % maxBeta = maxBeta ./ std(maxBeta, 0, 2);
            % maxBeta = maxBeta ./ max(abs(maxBeta), [], 2);
            % maxBeta(isnan(maxBeta)) = 0;
            
            % Hierarchical clustering
            D = pdist(maxBeta, 'correlation');
            Z = linkage(D, 'average');
            cc = cophenet(Z, D);
            fprintf("Cophenetic correlation coefficient is %.2f\n", cc);
            hcId = cluster(Z, 'maxclust', nc);
            
            % UMAP embedding
            umapCoords = run_umap(maxBeta, 'metric', 'correlation', 'verbose', 'none', 'min_dist', 0.05, 'randomize', false);
            
            % Output
            s.feats = mdl.feats;
            s.D = D;
            s.Z = Z;
            s.cc = cc;
            s.nc = nc;
            s.hcId = hcId;
            s.umapCoords = umapCoords;
        end
        
        function [s, d] = NameClusters(s, verName)
            % Add cluster names to data in the struct
            % 
            %   [s, d] = NameClusters(s, verName)
            % 
            
            switch verName
                case 'combo1-1'
                    % 12/28/2023 10:00pm: xall, nc=15
                    d = struct;
                    d.sent_onset = 5;
                    d.inten = 7;
                    d.voicing = 6;
                    d.pitch = 12;
                    d.fri_den = 2;
                    d.plo_vel = 13;
                    d.nasal = 10;
                    d.apr_lab = 8;
                    d.glottal = 11;
                    d.mid = 14;
                    d.others = [4 1 3 9 15];
                    
                case 'combo3pm-1'
                    % 12/28/2023 10:00pm: xall, nc=15
                    d = struct;
                    d.sent_onset = 12;
                    d.inten = 2;
                    d.voicing = 10;
                    d.pitch = 9;
                    d.accent = 14;
                    d.phrase = 3;
                    d.fri_den = 7;
                    d.blick = 6;
                    d.others = [15 1 8 4 5 11 13];
                    
                otherwise
                    error("%s is not a valid version name", verName);
            end
            s.hcIdDict = d;
            
            catNames = categorical(fieldnames(d), fieldnames(d), 'Ordinal', true);
            catMembers = struct2cell(d);
            iCat = zeros(size(s.hcId));
            for i = 1 : numel(s.hcId)
                iCat(i) = find(cellfun(@(x) any(s.hcId(i)==x), catMembers), 1);
            end
            s.hcName = catNames(iCat);
        end
        
        function avgTb = AverageModels(mdlTb, mdlNames, groupVar, groupList)
            % Compute average model weights across units
            % 
            %   avgTb = AverageModels(mdlTb, groupVar, mdlNames)
            % 
            
            if ~exist('groupList', 'var')
                groupList = unique(mdlTb.(groupVar));
            end
            nG = numel(groupList);
            avgTb = mdlTb(1:nG,:);
            avgTb.clusId = string(groupList);
            
            for i = 1 : nG
                m = mdlTb.(groupVar) == groupList(i);
                if sum(m) == 0
                    for j = 1 : numel(mdlNames)
                        mn = mdlNames{j};
                        avgTb.(mn){i} = [];
                    end
                    continue
                end
                
                % Use info of the first unit in this group to fill the columns
                k = find(m, 1);
                avgTb(i,2:end) = mdlTb(k,2:end);
                
                % Average models
                for j = 1 : numel(mdlNames)
                    mn = mdlNames{j};
                    
                    bb = cellfun(@(x) x.Beta .* (x.null.BetaPval<0.05), mdlTb.(mn)(m), 'Uni', false);
                    bb = cat(2, bb{:});
                    bb = normalize(bb, 1);
                    mb = mean(bb, 2, 'omitnan');
                    
                    r2 = cellfun(@(x) x.r2, mdlTb.(mn)(m));
                    mr2 = mean(r2);
                    
                    avgTb.(mn){i}.Beta = mb;
                    avgTb.(mn){i}.r2 = mr2;
                    avgTb.(mn){i}.resps = groupList(i);
                end
            end
        end
        
        % FPA
        function s = ComputePairwiseCorr(mdlTb)
            % 
            
            % Compute pairwise similarity matrix
            R = corr(mdlTb.maxBeta');
            
            % Compute pairwise distance matrix
            D = abs(mdlTb.depth - mdlTb.depth');
            
            % Take unique triangular portion of the matrices excluding the diagonal
            U = triu(true(size(R)), 1);
            R = R(U);
            D = D(U);
            
            s.region = mdlTb.region(1);
            s.mdlName = mdlTb.mdlName(1);
            s.recId = mdlTb.recId(1);
            s.d = D;
            s.r = R;
        end
        
        function s = ComputeSpatialSimilarity(s, nBoot)
            % Compute TRF similarity as a function of spatial distance
            % 
            %   s = ComputeSpatialSimilarity(mdlTb)
            % 
            % Inputs
            %   mdlTb           A clusTb like table with maxBeta and depth column.
            % 
            
            D = s.d;
            R = s.r;
            
            % Bin the distances
            dEdges = 0 : 100 : 2000;
            dCenters = MMath.BinEdges2Centers(dEdges);
            [~, ~, bins] = histcounts(D, dEdges);
            
            % Exclude out of range samples
            isOut = bins == 0;
            bins(isOut) = [];
            R(isOut) = [];
            D(isOut) = [];
            if isempty(bins)
                s = [];
                return
            end
            
            % Compute mean similarity in the bins
            [G, binInd] = findgroups(bins);
            dCenters = dCenters(binInd);
            sim = splitapply(@MMath.MedianStats, R, G);
            
            s.dEdges = dEdges;
            s.dCenters = dCenters;
            s.r = R;
            s.d = D;
            s.sim = sim;
            
            % Compute null
            if ~nBoot
                return
            end
            rng(61);
            null = zeros(numel(sim), nBoot);
            for i = 1 : nBoot
                null(:,i) = splitapply(@MMath.MeanStats, R, randsample(G, numel(G)));
            end
            nullMed = MMath.MedianStats(null, 2);
            simPval = MMath.EstimatePval(sim, null');
            
            s.null = nullMed;
            s.nullCI = prctile(null, [2.5 97.5], 2);
            s.simPval = simPval;
        end
        
        % Plot inputs
        function PlotInputs(ax, ce, feats, tWin, varargin)
            % 
            
            % p = inputParser;
            % p.addParameter('Normalization', 'none');
            
            feats = string(feats);
            [t, F, featEx] = ce.GetArray('feat', [], feats, varargin{:});
            isWin = t > tWin(1) & t < tWin(2);
            t = t(isWin);
            t = t - t(1);
            F = F(isWin,:);
            y = 1 : size(F,2);
            
            if isscalar(feats) && ismember(feats, ["mic", "speaker1"])
                CData = 1 - F';
                CData(isnan(CData)) = 1;
                
                imagesc(t, y, flip(CData));
                ax.YTick = mean(y);
                ax.YTickLabel = "Mel";
                
            elseif any(ismember(feats, NP.Artic.featNames))
                CData = -F';
                imagesc(t, y, CData);
                ax.CLim = prctile(CData(:), [.5 99.5]);
                ax.YTick = y;
                lb = NP.Artic.GetLabels(featEx);
                lb(2:2:end) = lb(2:2:end) + "        ";
                ax.YTickLabel = lb;
                
            else
                CData = 1 - F';
                CData(isnan(CData)) = 1;
                
                imagesc(t, y, CData);
                ax.YTick = y;
                lb = MLing.ARPA2IPA(featEx);
                lb(2:2:end) = lb(2:2:end) + "        ";
                ax.YTickLabel = lb;
            end
            colormap gray
            ax.XTick = [];
            ax.XTickLabel = [];
            % ax.XLabel.String = "Time (s)";
            MPlot.Axes(ax);
        end
        
        % Plot models
        function PlotWeightsCombo(clusTb, mdlNames, varargin)
            % Plot a grid of model weights
            % 
            %   PlotWeightsCombo(clusTb, mdlNames)
            %   PlotWeightsCombo(..., 'UnitIds', [])
            %   PlotWeightsCombo(..., 'UnitsPerPage', 15)
            %   PlotWeightsCombo(..., 'Page', 0)
            %   PlotWeightsCombo(..., 'Folder', [])
            % 
            % Inputs
            %   clusTb          Cluster table with models.
            %   mdlNames        A m-by-n cell array of model variable names in clusTb.
            % 
            
            p = inputParser;
            p.addParameter('UnitIds', [], @isnumeric);
            p.addParameter('UnitsPerPage', 15, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Page', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('PanelArgs', {}, @iscell);
            p.addParameter('Folder', [], @(x) ischar(x) || isstring(x));
            p.parse(varargin{:});
            uId = p.Results.UnitIds;
            upp = p.Results.UnitsPerPage;
            pg = p.Results.Page;
            panelArgs = p.Results.PanelArgs;
            figDir = p.Results.Folder;
            
            % Determine the units to plot
            if isempty(uId)
                % By page index and # of units per page
                pg = max(pg, 1);
                uInd = (1:upp)+upp*(pg-1);
                uInd(uInd > height(clusTb)) = NaN;
                uId = NP.Unit.Ind2ClusId(uInd, clusTb)';
            else
                % By unit ID
                uInd = NP.Unit.ClusId2Ind(uId, clusTb);
                uInd = sort(uInd);
                upp = numel(uInd);
            end
            uInd = uInd(:)';
            
            % Configure plots
            uIndGrid = NaN(numel(mdlNames), upp);
            uIndGrid(:,1:numel(uInd)) = repmat(uInd, [numel(mdlNames) 1]);
            mdlNameGrid = repmat(mdlNames(:), [1 upp]);
            panelArgsGrid = repmat({panelArgs}, size(mdlNameGrid));
            rowDist = ones(size(mdlNames));
%             for i = 1 : numel(mdlNames)
%                 if mdlNames{i} == "phase"
%                     rowDist(i) = 1;
%                 else
%                     rowDist(i) = 3;
%                 end
%             end
            
            % Plotting
            f = gcf;
            f.WindowState = 'maximized';
            LMV.PlotClus.Grid(clusTb, uIndGrid, mdlNameGrid, rowDist, [], 'PanelArgs', panelArgsGrid);
            
            % Save figure
            if isempty(figDir)
                return
            end
            if ~exist(figDir, 'dir')
                mkdir(figDir);
            end
            uId(isnan(uId)) = [];
            uId = uId(:)';
            figName = strjoin([clusTb.recId(1), "page"+pg, "u"+uId], "_");
            exportgraphics(f, fullfile(figDir, figName+".png"));
        end
        
        function PlotWeightsComboFromCache(clusId, cacheDir, mdlNames)
            % Plot a grid of model weights from unit cache
            % 
            %   PlotWeightsComboFromCache(clusId, cacheDir, mdlNames)
            % 
            
            % Load unit cache
            if nargin < 2 || isempty(cacheDir)
                cacheDir = "m1";
            end
            clusId(isnan(clusId)) = [];
            uCell = NP.Unit.LoadUnitCache(clusId, 'DataSource', cacheDir);
            
            qcTb = struct2table(cellfun(@(x) x.unitInfo, uCell), "AsArray", true);
            
            trfTb = struct2table(cellfun(@(x) x.trf, uCell), "AsArray", true);
            for i = 1 : width(trfTb)
                if ~iscell(trfTb.(i))
                    trfTb.(i) = num2cell(trfTb.(i)); % standardize structs in cell array
                end
            end
            
            uTb = [qcTb trfTb];

            if nargin < 3
                mdlNames = trfTb.Properties.VariableNames(1:3);
            end
            
            % Arrange units by figure pages
            upp = min(numel(clusId), 20);
            nPages = ceil(height(uTb) / upp);
            pInd = repelem(1:nPages, upp);
            uTb.pInd = pInd(1:height(uTb))';
            
            % Plot figures
            for p = 1 : nPages
                MPlot.Figure(1000+p); clf
                LMV.TRF.PlotWeightsCombo(uTb, mdlNames, 'UnitsPerPage', upp, 'Page', p);
            end
        end
        
        function PlotWeights2(mdl, varargin)
            % Plot TRF weights as a heatmap
            % 
            %   PlotWeights2(mdl)
            %   PlotWeights2(mdl, ..., 'WeightAlpha', [])
            %   PlotWeights2(mdl, ..., 'ClusInfo', [])
            %   PlotWeights2(mdl, ..., 'Parent', gca)
            % 
            
            p = inputParser;
            p.addParameter('WeightAlpha', [], @(x) isscalar(x) && isnumeric(x));
            p.addParameter('ClusInfo', [], @(x) istable(x) || isstruct(x));
            p.addParameter('Parent', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
            clusInfo = p.Results.ClusInfo;
            wAlpha = p.Results.WeightAlpha;
            ax = p.Results.Parent;
            
            if isempty(ax)
                ax = gca;
            end
            
            if iscell(mdl)
                mdl = mdl{1};
            end
            if isempty(mdl)
                axis(ax, 'off');
                return
            end
            
            if istable(clusInfo)
                clusInfo = table2struct(clusInfo);
            end
            
            % Make weight matrix
            b = mdl.Beta;
            k = std(b)*5; % for color scaling
            
            % Apply statistics from null
            sigStr = '';
            if isfield(mdl, 'null')
                if ~isempty(wAlpha)
                    b(mdl.null.BetaPval > wAlpha) = 0;
                end
                r2Sig = sum(mdl.null.r2Pval < [0.05 0.01 0.001]);
                sigStr = [' ' repelem('*', r2Sig)];
            end
            
            % mdl.feats(mdl.feats=="time") = []; % posthoc fix
            bm = reshape(b, numel(mdl.feats), []);
            bm = flip(bm, 2);
            if startsWith(mdl.name, "strf")
                bm = flip(bm,1);
            end
            t = -flip(mdl.dt);
            tEdges = MMath.BinCenters2Edges(t);
            
            % Plot weights as a heatmap
            imagesc(ax, t, [], bm); hold(ax, 'on');
            
            if size(bm,1) < 50
                tTxt = repelem(t(1), size(bm,1));
                yTxt = 1 : size(bm,1);
                lb = NP.Artic.GetLabels(mdl.feats);
                text(ax, tTxt, yTxt, lb, 'VerticalAlignment', 'middle', 'FontSize', 8);
            end
            
            if isfield(mdl, 'r2each')
                % % Reselect r2
                % [mdl.r2, mdl.r2t, mdl.r2Idx] = LMV.RF.FindPeakR2(mdl.r2each, mdl.dt);
                
                % Plot r2 time course for time independent model
                h = size(bm, 1);
                r = flip(mdl.r2each) * h; % scale up to the height of the heatmap
                r = h + .5 - r;
                plot(t, r, 'k');
                
                isPk = flip(MMath.Ind2Logical(mdl.r2Idx, numel(r)));
                plot(t(isPk), h+.5 - mdl.r2*h, 'k*');
            end
            
            ax.Box = 'off';
            ax.YTick = [];
            ax.Colormap = flip(brewermap([], 'RdBu'));
            ax.CLim = [-1 1] * max(k, eps);
            ax.XLim = tEdges([1 end]);
            ax.YDir = 'reverse';
            ax.YLim = [0.5 size(bm,1)+0.5];
            ax.YLabel.String = sprintf("%s, %sr2 = %.2f", mdl.name, sigStr, mdl.r2);
            ax.YLabel.Interpreter = 'none';
            ax.Title.String = mdl.resps;
            ax.Title.Interpreter = 'none';
            if sigStr == " "
                ax.Title.FontWeight = 'normal';
            else
                ax.YLabel.FontWeight = 'bold';
            end
        end
        
        function PlotR2Scatter(sets, targets, clusTb, r2Tb, uSig)
            % Plot goodness-of-fit in violin scatter plots
            % 
            %   PlotR2Scatter(sets, targets, clusTb, r2Tb, uSig)
            % 
            
            % Find regions and targets
            mdlNames = sets+"_"+targets;
            regions = unique(clusTb.region);
            
            nRegion = numel(regions);
            nModel = numel(mdlNames);
            
            % Set up plot
            ax = nexttile;
            
            % Preallocate variables to collect the plotted data for labeling and callback
            hh = cell(nModel, nRegion);
            tt = hh;
            xx = hh;
            yy = hh;
            
            for j = 1 : nRegion
                for i = 1 : nModel
                    % Get r2 values
                    mn = mdlNames(i);
                    r2 = r2Tb.(mn);
                    
                    isVal = ~isnan(r2);
                    isRegion = clusTb.region == regions(j);
                    isSig = uSig(:,i) > 0;
                    
                    isUnit = isVal & isRegion & isSig;
                    r = r2(isUnit);
                    
                    x = i + 1/(nRegion+3)*(j-2.5);
                    binEdges = 0:.05:1;
                    hh{i,j} = MPlot.ViolinScatter(x, r, binEdges, 'Span', 0.11, 'Color', NP.Param.GetRegionColors(regions(j)));
                    hold(ax, 'on');
                    
                    % Collect plotting data
                    if isempty(hh{i,j})
                        continue
                    end
                    xx{i,j} = hh{i,j}.XData';
                    yy{i,j} = hh{i,j}.YData';
                    tt{i,j} = clusTb(isUnit,:);
                end
            end
            
            ax.XLim = [.5 nModel+.5];
            ax.YLim = [0 .4];
            ax.XTick = 1 : nModel;
            ax.XTickLabel = targets;
            ax.Title.String = sprintf("Goodness-of-fit of optimal '%s' RF", sets);
            ax.Title.Interpreter = "none";
            ax.YLabel.String = "r-squared";
            MPlot.Axes(ax);
            
            % Set up brushing callback
            X = cat(1, xx{:});
            Y = cat(1, yy{:});
            hCB = plot(X, Y, 'w.');
            hCB.UserData.forBrush = true;
            ax.Children = circshift(ax.Children, -1); % move callback plot to the bottom
            ax.UserData.clusTb = cat(1, tt{:});
            ax.UserData.taskPhase = "full";
            brushObj = brush(gcf);
            brushObj.ActionPostCallback = @LMV.Overview.SessionOnBrush;
            
            % % Formatting
            % ldg = legend([hh{1,:}], regions);
            % ldg.Location = 'eastoutside';
            % ldg.EdgeColor = 'none';
        end
        
        % Plot clustering
        function PlotClusterCombo(ss)
            % Plot clustered max feature weights with dendrogram
            % 
            %   PlotClusterCombo(ss)
            % 
            
            % Create layout
            nMaps = numel(ss);
            rowDist = [12 58 30];
            tl = tiledlayout(sum(rowDist), nMaps);
            tl.Padding = 'compact';
            
            for k = 1 : nMaps
                s = ss{k};
                if numel(s.resps) <= s.nc
                    continue % not enough units
                end
                
                % Get optimal leaf order
                Z = s.Z;
                leafOrder = optimalleaforder(Z, s.D, 'Criteria', 'adjacent');
                cid = s.hcId(leafOrder);
                B = s.maxBeta(leafOrder,:);
                coords = s.umapCoords(leafOrder,:);
                
                % Define cluster colors
                cidList = unique(cid, 'stable');
                nClus = numel(cidList);
                cc = lines(nClus);

                % Plot dendrogram
                ntArgs = MPlot.FindTileInd(rowDist, nMaps, 1, k);
                ax = nexttile(ntArgs{:});
                nClusCut = min(20, size(Z,1)-1);
                cutoff = median([Z(end-nClusCut+1,3) Z(end-nClusCut+2,3)]);
                hh = dendrogram(Z, Inf, 'Reorder', leafOrder);
                for i = 1 : numel(hh)
                    hh(i).Color = 'k';
                end
                hold(ax, 'on');
                wRib = (ax.YLim(2)-cutoff)*.3;
                [cx, cy] = MPlot.GroupRibbon(cid, cutoff-0.1+[-wRib 0], cc, 'Groups', cidList);
                ax.XLim = [0.5 numel(cid)+0.5];
                ax.YLim(1) = cutoff-0.1-wRib;
                ax.XTick = cx;
                ax.XTickLabel = string(cidList);
                ax.XTickLabelRotation = 0;
                ax.YTick = [];
                ax.Box = 'off';
                ax.Title.String = s.regions(1) + ": " + s.mdlNames(1);
                ax.Title.Interpreter = 'none';
                
                % Plot weight heatmap
                ntArgs = MPlot.FindTileInd(rowDist, nMaps, 2, k);
                ax = nexttile(ntArgs{:});
                LMV.TRF.PlotClusterHeatmap(ax, B, s.feats);

                % Plot embedding
                ntArgs = MPlot.FindTileInd(rowDist, nMaps, 3, k);
                ax = nexttile(ntArgs{:});
                LMV.TRF.PlotClusterEmbedding(ax, coords, [], cid, cidList);
            end
        end
        
        function PlotXModelClusterCombo(s)
            % Plot clustered max feature weights
            % 
            %   PlotXModelClusterCombo(s)
            % 
            
            % Get variables with optimal leaf order
            Z = s.Z;
            leafOrder = optimalleaforder(Z, s.D, 'Criteria', 'adjacent');
            cid = s.hcId(leafOrder);
            B = s.maxBeta(leafOrder,:);
            coords = s.umapCoords(leafOrder,:);
            
            % Find target membership
            mInd = find(s.hasMdl);
            mInd = mInd(leafOrder);
            [uInd, tInd] = ind2sub(size(s.hasMdl), mInd);
            
            % Define cluster colors
            cidList = unique(cid, 'stable');
            
            % Create layout
            rowDist = [65 35];
            mdlNames = unique(s.mdlNames, 'stable');
            nMdls = numel(mdlNames);
            tl = tiledlayout(sum(rowDist), nMdls);
            tl.Padding = 'compact';
            
            for k = 1 : nMdls
                % Mask model
                m = tInd == k;
                subB = B(m,:);
                subCID = cid(m);
                subCoords = coords(m,:);
                
                % Plot weight heatmap
                ntArgs = MPlot.FindTileInd(rowDist, nMdls, 1, k);
                ax = nexttile(ntArgs{:});
                LMV.TRF.PlotClusterHeatmap(ax, subB, s.feats, subCID, cidList);
                ax.Title.String = s.regions(1) + ": " + mdlNames(k);
                ax.Title.Interpreter = 'none';
                
                % Plot embedding
                ntArgs = MPlot.FindTileInd(rowDist, nMdls, 2, k);
                ax = nexttile(ntArgs{:});
                LMV.TRF.PlotClusterEmbedding(ax, subCoords, coords, subCID, cidList);
            end
        end
        
        function PlotXRegionClusterCombo(ss, region)
            % Plot clustered max feature weights
            % 
            %   PlotXRegionClusterCombo(ss, region)
            % 
            
            % Create layout
            nMdls = numel(ss);
            rowDist = [65 35];
            tl = tiledlayout(sum(rowDist), nMdls);
            tl.Padding = 'compact';
            
            for k = 1 : nMdls
                % Get variables with optimal leaf order
                s = ss{k};
                Z = s.Z;
                leafOrder = optimalleaforder(Z, s.D, 'Criteria', 'adjacent');
                cid = s.hcId(leafOrder);
                B = s.maxBeta(leafOrder,:);
                coords = s.umapCoords(leafOrder,:);
                regNames = s.regions(leafOrder);
                
                % Define cluster members
                cidList = unique(cid, 'stable');
                
                % Mask region
                m = regNames == region;
                subB = B(m,:);
                subCID = cid(m);
                subCoords = coords(m,:);
                
                % Plot weight heatmap
                ntArgs = MPlot.FindTileInd(rowDist, nMdls, 1, k);
                ax = nexttile(ntArgs{:});
                LMV.TRF.PlotClusterHeatmap(ax, subB, s.feats, subCID, cidList);
                ax.Title.String = region + ": " + s.mdlNames(1);
                ax.Title.Interpreter = 'none';
                
                % Plot embedding
                ntArgs = MPlot.FindTileInd(rowDist, nMdls, 2, k);
                ax = nexttile(ntArgs{:});
                LMV.TRF.PlotClusterEmbedding(ax, subCoords, coords, subCID, cidList);
            end
        end
        
        function PlotXAllClusterCombo(s, region)
            % Plot clustered max feature weights
            % 
            %   PlotXAllClusterCombo(s, region)
            % 
            
            % Get variables with optimal leaf order
            Z = s.Z;
            leafOrder = optimalleaforder(Z, s.D, 'Criteria', 'adjacent');
            cid = s.hcId(leafOrder);
            B = s.maxBeta(leafOrder,:);
            coords = s.umapCoords(leafOrder,:);
            regNames = s.regions(leafOrder);
            
            % Find target membership
            mInd = find(s.hasMdl);
            mInd = mInd(leafOrder);
            [uInd, tInd] = ind2sub(size(s.hasMdl), mInd);
            
            mdlNames = unique(s.mdlNames, 'stable');
            nMdls = numel(mdlNames);
            
            % Define cluster colors
            cidList = unique(cid, 'stable');
            
            % Create layout
            rowDist = [60 40];
            tl = tiledlayout(sum(rowDist), nMdls);
            tl.Padding = 'compact';
            
            for k = 1 : nMdls
                % Mask model
                m = tInd == k & regNames == region;
                subB = B(m,:);
                subCID = cid(m);
                subCoords = coords(m,:);
                
                % Plot weight heatmap
                ntArgs = MPlot.FindTileInd(rowDist, nMdls, 1, k);
                ax = nexttile(ntArgs{:});
                LMV.TRF.PlotClusterHeatmap(ax, subB, s.feats, subCID, cidList);
                ax.Title.String = region + ": " + mdlNames(k);
                ax.Title.Interpreter = 'none';
                
                % Plot embedding
                ntArgs = MPlot.FindTileInd(rowDist, nMdls, 2, k);
                ax = nexttile(ntArgs{:});
                LMV.TRF.PlotClusterEmbedding(ax, subCoords, coords, subCID, cidList);
            end
        end
        
        function PlotClusterHeatmap(ax, B, featNames, cid, cidList)
            % 
            
            % Normalize heatmap
            B = B ./ std(B, 0, 2);
            B(isnan(B)) = 0;
            
            % Plot heatmap
            imagesc(ax, B'); hold(ax, "on");
            
            if nargin > 3
                % Plot cluster indicators
                nClus = numel(cidList);
                cc = lines(nClus);
                [cx, cy] = MPlot.GroupRibbon(cid, [-.5 0], cc, 'Groups', cidList);
                text(ax, cx, cy-.25, string(cidList), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
                ax.YLim(1) = -1.5;
            end
            
            ax.Colormap = flip(brewermap([], 'RdBu'));
            ax.CLim = [-1 1] * prctile(abs(B(:)), 99.8);
            ax.YTick = 1: numel(featNames);
            ax.YTickLabel = featNames;
            ax.TickLabelInterpreter = 'none';
            ax.Box = 'off';
        end
        
        function PlotClusterEmbedding(ax, coords, bkCoords, cid, cidList)
            % 
            % 
            %   PlotClusterEmbedding(ax, coords, bkCoords, cid, cidList)
            % 
            
            % Plot all units as background
            if isempty(bkCoords)
                bkCoords = coords;
            end
            plot(ax, bkCoords(:,1), bkCoords(:,2), '.', 'Color', [0 0 0]+.9); hold(ax, 'on');
            
            % Plot units from each cluster
            nClus = numel(cidList);
            cc = lines(nClus);
            for i = 1 : nClus
                m = cid == cidList(i);
                plot(ax, coords(m,1), coords(m,2), '.', 'Color', cc(i,:));
            end
            
            % Label cluster IDs at cluster centers
            [G, cidText] = findgroups(cid);
            mCoords = splitapply(@(x) median(x,1), coords, G);
            text(ax, mCoords(:,1), mCoords(:,2), string(cidText));
            
            % ax.XLabel.String = 'Units';
            ax.Box = 'off';
            axis(ax, 'off');
        end
        
        % Plot FPA
        function PlotEmbeddingFPA(sHC, clusTb, colorByList, regionList, mdlList)
            % 
            % 
            %   PlotEmbeddingFPA(sHC, clusTb, colorByList, regionList, mdlList)
            % 
            
            uMask = true(size(sHC.hcId));
            if exist('regionList', 'var') && ~isempty(regionList)
                uMask = uMask & ismember(sHC.regions, regionList);
            end
            if exist('mdlList', 'var') && ~isempty(mdlList)
                uMask = uMask & ismember(sHC.mdlNames, mdlList);
            end
            clusTb = clusTb(uMask,:);
            
            % Put clusTb in a ce for input compatibility
            ce = NP.CodingExplorer;
            ce.clusTb = clusTb;
            ce.userData.hcIdList = unique(sHC.hcId, 'stable');
            
            tl = tiledlayout('flow');
            tl.Padding = 'compact';
            
            colorByList = string(colorByList);
            for i = 1 : numel(colorByList)
                ax = nexttile;
                NP.Embed.PlotScatter(ce, colorByList(i));
            end
        end
        
        function PlotSpatialSimilarity(sAll, sCell)
            % 
            % 
            %   PlotSpatialSimilarity(sAll, sCell)
            % 
            
            if ~exist('sCell', 'var')
                sCell = {};
            end
            
            cc = lines(numel(sCell));
%             cc = brewermap(numel(sCell), 'Set1');
            for i = 1 : numel(sCell)
                s = sCell{i};
%                 cc(i,:) = 1 - (1-cc(i,:))*0.8; % make paler
                hh(i) = plot(s.d, s.r, '.', 'Color', cc(i,:)); hold on
            end
            
            s = sAll;
            x = s.dCenters;
%             plot(s.d, s.r, '.', 'Color', [0 0 0]+.8);
            plot(x, s.sim, 'Color', 'k', 'LineWidth', 2);
            plot(x, s.null, 'k');
            MPlot.ErrorShade(x, s.null, s.nullCI(:,1), s.nullCI(:,2), 'IsRelative', false);
            
            ax = gca;
            ax.XLabel.String = "Distance (um)";
            ax.YLabel.String = "Pearson's r";
            ax.Title.String = s.region + ", " + s.mdlName;
            ax.Title.Interpreter = 'none';
            ax.Box = 'off';
            
            if isempty(sCell)
                return
            end
            recIds = cellfun(@(x) x.recId, sCell);
            ldg = legend(hh, recIds);
            ldg.Interpreter = 'none';
            ldg.Location = 'eastoutside';
        end
        
    end
    
end

