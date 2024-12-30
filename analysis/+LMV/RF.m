classdef RF
    
    properties(Constant)
    end
    
    methods(Static)
        % Model fitting
        function clusTb = FitSlidingTimeRF(ce, uInd, varargin)
            % Independtly fit a series of models with different time offsets
            % 
            %   clusTb = FitSlidingTimeRF(ce, uInd, 'FeatSet', '', 'Target', '')
            %   clusTb = FitSlidingTimeRF(ce, [], 'FeatSet', '', 'Target', '')
            %   clusTb = FitSlidingTimeRF(..., 'Lambda', [])
            %   clusTb = FitSlidingTimeRF(..., 'NBoot', [])
            % 
            % Inputs
            %   ce          An NP.CodingExplorer object including a clusTb in its metadata.
            %   uInd        The indices of units to fit. Empty value [] includes all units.
            %   'FeatSet'   Names of a set of features.
            %   'Target'    Name of the target time period.
            %   'Lambda'    L2 regularization strength.
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
            
            ss = statset('UseParallel', true);
            % ss = statset('UseParallel', false);
            
            % Get responses
            [~, Y, resps] = ce.GetArray('resp', [], uInd+1);
            
            % Get sample mask
            [~, M] = ce.GetArray('feat', [], poi, 'TimeShifts', dt);
            isPhase = any(M, 2);
            
            % Normalize responses
            Y(~isPhase,:) = NaN;
            Y = MMath.Normalize(Y, 'zscore');
            Y(isnan(Y)) = 0; % handle 0 devided by 0
            
            % Print info
            mdlName = featSet+"_"+target;
            fprintf("\nFit '%s' models.\n", mdlName);
            fprintf("\n%s\n", NP.SE.GetID(ce));
            fprintf("%i predictors, %i responses\n", numel(feats), numel(uInd));
            fprintf("%i raw samples, %i in selected task phase\n", numel(isPhase), sum(isPhase));
            
            % Fit models with different time shifts in features
            dtMdls = cell(numel(uInd), numel(dt));
            for i = 1 : numel(dt)
                % Get predictors
                [~, X, featEach] = ce.GetArray('feat', [], feats, 'TimeShift', dt(i));
                isVal = any(X, 2);
                
                % Normalize predictors
                X(~isPhase,:) = NaN;
                X = MMath.Normalize(X, 'zscore');
                X(isnan(X)) = 0; % handle 0 devided by 0
                
                % Find samples both in phase and non-zero in any feature
                m = isPhase & isVal;
                
                % Fit models
                if size(X,2) > 30
                    if i == 1
                        fprintf("Use reduced bootstraping null distributions since the model is too large.\n");
                    end
                    nShifts = min(nShifts, 50);
                end
                fprintf("predictor time shift %i/%i (%.2fs)\n", i, numel(dt), dt(i));
                dtMdls(:,i) = LMV.Fit.BatchLasso(X(m,:), Y(m,:), ...
                    'Alpha', alpha, 'Lambda', L, 'Options', ss, 'NumShifts', nShifts, 'MinShift', round(0.1*sum(m)), ...
                    'Verbose', false);
            end
            
            % Combine models
            mdls = cell(numel(uInd), 1);
            for i = 1 : numel(mdls)
                % Find the best model at the largest peak in r2 values
                r2each = cellfun(@(x) x.r2, dtMdls(i,:));
                [~, loc] = findpeaks(r2each, 'SortStr', 'descend', 'NPeaks', 1);
                if isempty(loc)
                    [~, loc] = max(r2each);
                end
                mdl = dtMdls{i,loc};
                
                % Update fields
                Beta = cellfun(@(x) x.Beta, dtMdls(i,:), 'Uni', false);
                mdl.Beta = cat(1, Beta{:});
                mdl.Bias = cellfun(@(x) x.Bias, dtMdls(i,:));
                mdl.r2Idx = loc;
                mdl.r2t = dt(loc);
                mdl.r2each = r2each;
                
                if isfield(mdl, 'null')
                    mdl.null.r2 = [];
                    BetaP = cellfun(@(x) x.null.BetaPval, dtMdls(i,:), 'Uni', false);
                    mdl.null.BetaPval = cat(1, BetaP{:});
                    mdl.null.BiasPval = cellfun(@(x) x.null.BiasPval, dtMdls(i,:));
                    mdl.null.r2PvalEach = cellfun(@(x) x.null.r2Pval, dtMdls(i,:));
                end
                
                % Add additional info
                mdl.name = mdlName;
                mdl.feats = featEach;
                mdl.resps = resps(i);
                mdl.dt = dt;
                mdl.mdls = dtMdls(i,:);
                
                mdls{i} = mdl;
            end
            
            % Get clusTb
            clusTb = NP.Unit.GetClusTb(ce);
            clusTb = NP.Unit.AddRecMeta(ce, clusTb);
            clusTb.(mdlName)(uInd) = mdls;
        end
        
        % Model utilities
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
        
        function mdls = UpdateModelR2(mdls)
            % Find and update peak r2 from for time-independent RF models
            % 
            %   mdls = UpdateModelR2(mdls)
            % 
            if iscell(mdls)
                mdls = cellfun(@LMV.RF.UpdateModelR2, mdls, 'Uni', false);
                return
            end
            if isempty(mdls) || ~all(ismember(["r2each", "dt"], fieldnames(mdls)))
                return
            end
            [mdls.r2, mdls.r2t, mdls.r2Idx] = LMV.RF.FindPeakR2(mdls.r2each, mdls.dt);
        end
        
        function [rPk, tPk, iPk] = FindPeakR2(r, t, tWin)
            % Find peak r2 value with the following criteria
            % 1) Peak time must be within the given window
            % 2) Peak must coincide with the maximum of the timeseries
            % 3) Peak must be greater than 0.01
            % 
            %   [rPk, tPk, iPk] = FindPeakR2(r, dt)
            %   [rPk, tPk, iPk] = FindPeakR2(r, dt, tWin)
            % 
            
            rNull = 0;
            tNull = MMath.Bound(mean(t([1 end])), t);
            
            [rPk, tPk, w, p] = findpeaks(r, t, "SortStr", "descend", "MinPeakProminence", 0.01);
            
            % Check peak times
            if nargin < 3
                [~, I] = max(abs(t));
                tWin = sort([0 sign(t(I))*0.3]);
            end
            isWin = tPk < max(tWin) & tPk > min(tWin);
            rPk = rPk(isWin);
            tPk = tPk(isWin);
            if isempty(rPk)
                rPk = rNull;
                tPk = tNull;
            end
            
            % Check peak value
            isWin = t >= tWin(1) & t <= tWin(2);
            if rPk(1) < 0.01 || rPk(1) ~= max(r(isWin))
                rPk = rNull;
                tPk = tNull;
            end
            
            rPk = rPk(1);
            tPk = tPk(1);
            iPk = find(t==tPk);
        end
        
        % Related to RF models
        function seqTb = MakeSeqTb(M, phSeqCell, se)
            % 
            
            % Get model info
            nameParts = strsplit(char(M.name), '_');
            [set, target] = nameParts{:};
            
            % Find unit index
            uid = str2double(erase(M.resps, "u"));
            uIdx = NP.Unit.ClusId2Ind(uid, se);
            if isnan(uIdx)
                fprintf("%s is not found in se\n", M.resps);
                seqTb = [];
                return
            end
            
            % Define time window
            switch target
                case 'stim'
                    rWin = [-0.1 0.4];
                case 'prod'
                    rWin = [-0.4 0.1];
                case 'feedback'
                    rWin = [-0.1 0.4];
                otherwise
                    error("'%s' is not a supported target name.", target);
            end
            rWin = rWin + [-0.05 0.05]; % add some padding
            fWin = [-0.25 0.25];
            
            % Find the top features
            m = M.mdls{M.r2Idx};
            b = m.Beta;
            b(m.null.BetaPval > 0.05) = 0;
            [~, I] = maxk(abs(b), 3);
            featNames = M.feats(I);
            
            % Extract responses and features
            for k = 1 : numel(phSeqCell)
                seqTb = phSeqCell{k}.seqTb;
                if isempty(seqTb)
                    continue
                end
                seqTb(seqTb.nSample<3,:) = [];
                
                seqTb = LMV.PE.AddResp2SeqTb(seqTb, se, rWin, uIdx, M);
                
                switch set
                    case "phone"
                        [~, ind] = MMath.SortLike(seqTb.seqStr, featNames, false);
                        seqTb = seqTb(ind,:);
                    case "strf"
                        seqTb = LMV.PE.AddFeat2SeqTb(seqTb, se, fWin, "mel", "speaker1");
                        seqTb.binIdx(:) = I(1);
                    case "artic3"
                        for i = 1 : numel(featNames)
                            fn = featNames(i);
                            if ismember(fn, se.GetTable("artic").Properties.VariableNames)
                                tn = "artic";
                            else
                                tn = "pitch";
                            end
                            seqTb = LMV.PE.AddFeat2SeqTb(seqTb, se, fWin, tn, fn, fn);
                        end
                end
                phSeqCell{k}.seqTb = seqTb;
            end
            
            seqTb = cellfun(@(x) x.seqTb, phSeqCell, 'Uni', false);
            seqTb = cat(1, seqTb{:});
        end
        
        % Plotting
        function PlotWeights(mdl, varargin)
            % Plot weight vectors for models fitted by lasso
            % 
            %   PlotWeights(mdl)
            %   PlotWeights(mdl, ..., 'ClusInfo', [])
            %   PlotWeights(mdl, ..., 'Orientation', 'horizontal')
            %   PlotWeights(mdl, ..., 'Parent', gca)
            % 
            
            p = inputParser;
            p.addParameter('ClusInfo', [], @(x) istable(x) || isstruct(x) || isempty(x));
            p.addParameter('Orientation', "horizontal", @(x) any(strcmpi(x, ["horizonal", "vertical"])));
            p.addParameter('SigColor', 'b', @(x) isnumeric(x) || ischar(x));
            p.addParameter('Parent', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
            clusInfo = p.Results.ClusInfo;
            ori = lower(p.Results.Orientation);
            sigColor = p.Results.SigColor;
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
            
            % Order features
            feats = mdl.feats;
            if startsWith(mdl.name, 'phone')
                feats = strrep(feats, 'UX', 'XX'); % such that the unused 'UX' is excluded from plotting
                [~, fInd] = MMath.SortLike(feats, ...
                    [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.labial, NP.Phone.coronal, NP.Phone.palatal, NP.Phone.velar, NP.Phone.glottal], ...
                    false); % by places
            else
                fInd = (1:numel(feats))';
            end
            feats = feats(fInd);
            
            % Make weight vector
            b = mdl.Beta(fInd);
            x = (1:numel(b))';
            
            % Get error
            sig = mdl.null.BetaPval(fInd) < 0.05;
            
            % Plot
            if ori == "horizontal"
                plot(ax, [0.5 numel(b)+0.5], [0 0]', 'Color', [0 0 0 .5]); hold(ax, 'on')
                stem(ax, x(~sig), b(~sig), '.', 'Color', [0 0 0]);
                stem(ax, x(sig), b(sig), '*', 'Color', sigColor);
            else
                plot(ax, [0 0]', [0.5 numel(b)+0.5], 'Color', [0 0 0 .5]); hold(ax, 'on')
                
                % plot(ax, b(~sig), x(~sig), '.', 'Color', [0 0 0]);
                plot(ax, [b(~sig)*0 b(~sig)]', [x(~sig) x(~sig)]', '-', 'Color', [0 0 0]);
                
                plot(ax, b(sig), x(sig), '.', 'Color', sigColor);
                plot(ax, [b(sig)*0 b(sig)]', [x(sig) x(sig)]', '-', 'Color', sigColor);
            end
            
            if ~isempty(clusInfo)
                ax.Title.String = sprintf("u%i, %ium", clusInfo.clusId, clusInfo.depth);
            else
                ax.Title.String = sprintf("%s", mdl.resps);
            end
            
            if ori == "horizontal"
                ax.YLabel.String = "Weights";
                ax.XLim = [0 numel(b)+1];
                ax.XTick = 1 : numel(b);
                ax.XTickLabel = feats;
            else
                ax.XLabel.String = "Weights";
                ax.YLim = [0 numel(b)+1];
                ax.YTick = 1 : numel(b);
                ax.YTickLabel = feats;
                ax.YDir = 'reverse';
            end
            
            ax.Box = 'off';
        end
        
        function uInd = PlotPopWeights(ax, tb, varargin)
            % 
            
            p = inputParser;
            p.addParameter('SortBy', 'none', @(x) isstring(x) || ischar(x));
            p.addParameter('MarkUnits', [], @isnumeric);
            p.parse(varargin{:});
            sortBy = string(p.Results.SortBy);
            uid = p.Results.MarkUnits;
            
            % Get model specs
            mdlName = tb.mdls{1}.name;
            mdlNameParts = strsplit(mdlName, '_');
            setName = mdlNameParts(1);
            targetName = mdlNameParts(2);
            
            % Get weight matrix
            bb = cellfun(@(x) x.mdls{x.r2Idx}.Beta, tb.mdls, 'Uni', false);
            B = cat(2, bb{:});
            B = zscore(B);
            feats = tb.mdls{1}.feats;
            
            % Get feature order
            fInd = (1:numel(feats))';
            switch setName
                case "phone"
                    feats = strrep(feats, 'UX', 'XX'); % such that the unused 'UX' is excluded from plotting
                    [feats, fInd] = MMath.SortLike(feats, ...
                        [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.labial, NP.Phone.coronal, NP.Phone.palatal, NP.Phone.velar, NP.Phone.glottal], ...
                        false); % by places
            end
            B = B(fInd,:);
            
            % Get unit order
            switch sortBy
                case "none"
                    uInd = (1:height(tb))';
                case "hc"
                    % Get unit order by clustering
                    D = pdist(B', "correlation");
                    Z = linkage(D);
                    uInd = optimalleaforder(Z, D, "Criteria", "adjacent");
                otherwise
                    % Get unit order by peak feature
                    [~, uInd] = LMV.NMF.SortPatterns(B, sortBy);
            end
            
            % Sort units
            tb.ind(uInd) = (1:height(tb))';
            tb.region = categorical(tb.region, LMV.Param.regions, 'Ordinal', true);
            [tb, uInd] = sortrows(tb, ["region", "ind"]);
            B = B(:,uInd);
            
            % Flip
            if setName ~= "strf"
                B = flip(B);
                feats = flip(feats);
            end
            
            
            % Plot weights as a heatmap
            imagesc(ax, B); hold(ax, 'on');
            colorbar;
            ax.Colormap = flip(brewermap([], 'RdBu'));
            % ax.CLim = [-1 1] * prctile(abs(B(:)), 99.5);
            ax.CLim = [-1 1] * 4;
            
            % Plot region indicator
            yTop = ax.YLim(2);
            [xc, yc] = MPlot.GroupRibbon(tb.region, yTop*[1.02 1.05], LMV.Param.GetRegionColors(LMV.Param.regions), 'Groups', LMV.Param.regions);
            text(xc, yc+yTop*0.04, LMV.Param.regions, "HorizontalAlignment", "center");
            
            % Highlight units
            cIdx = find(ismember(tb.clusId, uid));
            if ~isempty(cIdx)
                MPlot.Blocks(cIdx+[-.5 .5], [0 numel(feats)+1], 'EdgeColor', 'k', 'LineStyle', "-", 'FaceColor', 'none');
            end
            
            ax.YLim = [-.5 yTop*1.1];
            ax.Title.String = setName + " " + targetName;
            ax.XTick = 1 : size(B,2);
            ax.XTickLabel = cellfun(@(x) x.resps, tb.mdls);
            ax.YTick = 1 : size(B,1);
            ax.YDir = "normal";
            switch setName
                case "phone"
                    ax.YTickLabel = MPlot.StaggerLabels(MLing.ARPA2IPA(feats), [-5 0]);
                case "artic3"
                    ax.YTickLabel = MPlot.StaggerLabels(NP.Artic.GetLabels(feats), [-8 0]);
                case "strf"
                    ax.YTick = [1:20:size(B,1) size(B,1)];
                    ax.YTickLabel = round(feats(ax.YTick)/1e3, 1);
                    ax.YLabel.String = "Frequency (Hz)";
                otherwise
                    ax.YTickLabel = feats;
            end
            ax.TickDir = "out";
            ax.Box = "off";
        end
        
        function PlotR2Timecourse(mdl, varargin)
            % Plot r2 as a function of the time offset between features and responses
            % 
            %   PlotR2Timecourse(mdl, varargin)
            % 
            
            if ~isfield(mdl, 'r2each')
                error("The input is not a sliding time model.");
            end
            
            if contains(mdl.name, "prod")
                m = mdl.dt <= 0 & mdl.dt >= -0.3;
            else
                m = mdl.dt >= 0 & mdl.dt <= 0.3;
            end
            
            dt = mdl.dt;
            r2 = mdl.r2each;
            k = mdl.r2Idx;
            
            null = cellfun(@(x) prctile(x.null.r2, [50 2.5 97.5]), mdl.mdls, 'Uni', false);
            null = cat(1, null{:});
            
            ax = gca;
            plot(dt(m), null(m,1), 'k:'); hold on
            MPlot.ErrorShade(dt(m), null(m,1), null(m,3), null(m,2), 'IsRelative', false);
            plot(dt(m), r2(m), 'k');
            plot(dt(k), r2(k), 'k*');
            dt = dt(m);
            ax.XLim = dt([1 end]);
            ax.YLim(1) = 0;
            ax.XLabel.String = "\Deltat in features (s)";
            ax.YLabel.String = "r2";
            ax.Title.String = mdl.resps+", "+mdl.name;
            ax.Title.Interpreter = "none";
            MPlot.Axes(ax);
        end
        
        function PlotR2Overlay(mdls, varargin)
            % Plot an overlay of r2 time course for sliding time RF models
            % 
            %   PlotR2Overlay(mdls)
            %   PlotR2Overlay(mdls, ..., 'Parent', nexttile)
            % 
            
            p = inputParser;
            p.addParameter('Parent', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
            ax = p.Results.Parent;
            
            % Collect data
            hasMdl = ~cellfun(@isempty, mdls);
            mdls = mdls(hasMdl);
            r2each = cell(size(mdls));
            r2 = NaN(size(mdls));
            r2t = NaN(size(mdls));
            for i = 1 : numel(mdls)
                mdl = mdls{i};
                r2each{i} = mdl.r2each;
                t = mdl.dt;
                r2(i) = mdl.r2;
                r2t(i) = mdl.r2t;
            end
            r2each = cat(1, r2each{:})';
            r2each(r2each<0) = 0;
            
            % Make plot
            if isempty(ax)
                ax = nexttile;
            end
            
            % dy = median(r2);
            dy = 0.05;
            nu = size(r2each,2);
            yPos = cumsum(ones(nu,1)*dy) - dy;
            cc = zeros(nu,3);
            [~, I] = max(r2);
            cc(I,:) = [.8 0 0];
            MPlot.PlotTraceLadder(t, r2each, yPos, 'ColorArray', cc); hold(ax, "on")
            h = plot(r2t, r2+yPos, 'k.');
            h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("clusId", cellfun(@(x) x.resps, mdls));
            
            ax.XLim = [t(1)+.1, t(end)-.1];
            % ax.YLim = [0 max(r2+yPos)]*1.01;
            ax.YLim = [0 dy*50];
            ax.XTick = ax.XLim(1):0.1:ax.XLim(2);
            ax.XTickLabelRotation = 0;
            ax.YTick = [];
            ax.XLabel.String = "Time from feature (s)";
            ax.YLabel.String = "r-squared";
            MPlot.Axes(ax);
        end
        
    end
    
end

