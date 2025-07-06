classdef LM
    methods(Static)
        % Modeling
        function clusTb = FitEncoding(ce, uInd, varargin)
            % Fit linear models that predict prod activity from temporal stim activity for each unit
            % 
            %   clusTb = FitEncoding(ce)
            %   clusTb = FitEncoding(ce, uInd)
            %   clusTb = FitEncoding(..., 'Lambda', [])
            %   clusTb = FitEncoding(..., 'NBoot', [])
            % 
            % Inputs
            %   ce          An NP.CodingExplorer object including a clusTb in its metadata.
            %   uInd        The indices of units to fit. Empty value [] includes all units.
            %   'Lambda'    A scalar for L2 regularization strength.
            % Output
            %   clusTb      
            % 
            
            p = inputParser;
            p.addParameter('Lambda', [], @isnumeric);
            p.addParameter('NBoot', 100, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Alpha', eps, @isnumeric);
            p.parse(varargin{:});
            L = p.Results.Lambda;
            nShifts = p.Results.NBoot;
            alpha = p.Results.Alpha;
            
            % Configure model fitting
            if ~exist('uInd', 'var') || isempty(uInd)
                uInd = (1 : ce.numResp)';
            end
            dt = -0.4 : 0.02 : 0.4;
            
            % ss = statset('UseParallel', true);
            ss = statset('UseParallel', false);
            
            % Get responses
            tv = ce.GetTable("taskValue");
            isProd = tv.phase == "prod";
            [~, Y, resps] = ce.GetArray('resp', isProd, uInd+1);
            
            % Make sample masks
            stim = ce.GetTable("taskTime").stim(~isProd);
            T = ce.GetArray('resp', ~isProd, [], 'DimCat', 0);
            M = cell(size(T));
            for i = 1 : numel(M)
                M{i} = stim(i).MaskTimestamps(T{i});
            end
            M = logical(cat(1, M{:}));
            
            % Normalize responses
            Y(~M,:) = NaN;
            Y = MMath.Normalize(Y, 'zscore');
            Y(isnan(Y)) = 0; % handle 0 devided by 0
            
            % Print info
            mdlName = "linker";
            fprintf("\nFit '%s' models.\n", mdlName);
            fprintf("\n%s\n", NP.SE.GetID(ce));
            fprintf("%i predictors for each of the %i responses\n", numel(dt), numel(uInd));
            fprintf("%i samples\n", sum(M));
            
            % Fit models with different time shifts in features
            mdls = cell(numel(uInd), 1);
            for i = 1 : numel(uInd)
                % Get predictors
                [~, X] = ce.GetArray('resp', ~isProd, uInd(i)+1, 'TimeShift', dt);
                
                % Normalize predictors
                X(~M,:) = NaN;
                [~, c, k] = MMath.Normalize(X(:), 'zscore');
                X = (X-c)./k;
                X(isnan(X)) = 0; % handle 0 devided by 0
                
                % Search for optimal lambda
                if isempty(p.Results.Lambda)
                    s = NP.Fit.BatchSmoothLinear(X(M,:), Y(M,i), ...
                        'Alpha', alpha, 'Lambda', logspace(0, 2, 20), 'Verbose', false);
                    % L = s{1}.mdl.Lambda1SE;
                    L = s{1}.mdl.LambdaMinMSE;
                end
                
                % Fit model
                mdls(i) = NP.Fit.BatchSmoothLinear(X(M,:), Y(M,i), ...
                    'Alpha', alpha, 'Lambda', L, 'NumShifts', nShifts, 'MinShift', round(0.1*sum(M)), 'Verbose', false);
            end
            
            % Add additional info
            for i = 1 : numel(mdls)
                mdls{i}.name = mdlName;
                mdls{i}.feats = string(1:numel(dt));
                mdls{i}.resps = resps(i);
                mdls{i}.dt = dt;
            end
            
            % Get clusTb
            clusTb = NP.Unit.GetClusTb(ce);
            clusTb = NP.Unit.AddRecMeta(ce, clusTb);
            clusTb.(mdlName)(uInd) = mdls;
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
                mdlName = mdlNames(i);
                mdlDir = fullfile(LMV.Data.GetAnalysisDir, "linker", mdlName, "mdls");
                
                % Find recordings
                if isempty(recIds)
                    mSearch = MBrowse.Dir2Table(fullfile(mdlDir, "*_clusTb.mat"));
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
                    load(fullfile(mdlDir, recIds(j)+"_clusTb.mat"), 'clusTb');
                    m = ~ismember(clusTb.Properties.VariableNames, cTb.Properties.VariableNames);
                    cTb = [cTb clusTb(:,m)];
                    cTbs{j} = cTb;
                end
            end
            
            % Return table if a single recording is requested
            if isscalar(recIds)
                cTbs = cTbs{1};
            end
        end
        
        function ComputeLinkingScores(ce)
            % Predict production activity using listening responses
            % 
            %   PredictProdResp(ce)
            % 
            % Output
            %   A timeSeries table named 'linker' will be added to the ce. It contains the following variables.
            %   time        NP.TGEvent objects serving as timestamps.
            %   pred        A t-by-u matrix of predicted production responses.
            %   score       Gaussian weighted moving covariance between prediction and production.
            %   
            
            % Find stim epochs
            [tt, tv] = ce.GetTable('taskTime', 'taskValue');
            isStim = tv.phase == "stim";
            
            % Get PETHs
            [T, Rs] = ce.GetArray('resp', isStim, [], 'DimCat', 0);
            
            % Get kernels
            Beta = cellfun(@(x) x.Beta, ce.clusTb.linker, 'Uni', false);
            Beta = cat(2, Beta{:});
            
            x = ce.clusTb.linker{1}.dt;
            xq = x(1) : ce.userData.rsOps.rsBinSize : x(end);
            Beta = interp1(x', Beta, xq');
            
            Bias = cellfun(@(x) x.Bias, ce.clusTb.linker);
            
            % Predict prod responses from stim responses
            Rh = Rs;
            for j = 1 : numel(Rs)
                r = Rs{j};
                for i = 1 : size(r,2)
                    r(:,i) = conv(r(:,i), Beta(:,i), 'same') + Bias(i);
                end
                Rh{j} = r;
            end
            
            % Resample prod response times to be consistent with stim
            rpTb = table;
            [rpTb.time, rpTb.r] = ce.GetArray('resp', ~isStim, [], 'DimCat', 0);
            tEdges = cellfun(@(x) MMath.BinCenters2Edges(x), T(~isStim), 'Uni', false);
            rpTb = ce.ResampleTimeSeries(rpTb, tEdges);
            Rp = rpTb.r;
            
            % Compute Gaussian weighted moving covariance between prediction and production
            fs = 1/ce.userData.rsOps.rsBinSize;
            S = cellfun(@(x,y) MNeuro.Filter1(x.*y, fs, 'gaussian', 0.1), Rh, Rp, 'Uni', false);
            
            % Make timeseries table
            lkTb = table;
            lkTb.time(isStim) = T;
            lkTb.stim(isStim) = Rs;
            lkTb.prod(isStim) = Rp;
            lkTb.pred(isStim) = Rh;
            lkTb.score(isStim) = S;
            ce.SetTable('linker', lkTb, 'timeSeries');
            
            % Normalize scores
            ce.SetColumn('linker', 'pred', @(x) MMath.Normalize(x, 'minmax'), 'all');
            ce.SetColumn('linker', 'score', @(x) MMath.Normalize(x, 'minmax'), 'all');
            
            % Resample timeseries at phone times
            ph = tt.stimPhn;
            tEdges = cellfun(@(x) [double(x); x(end).T.tmax], ph, 'Uni', false);
            tEdges(~isStim) = cell(sum(~isStim), 1);
            lkTb2 = ce.ResampleTimeSeries("linker", tEdges);
            lkTb2.time(isStim) = ph(isStim);
            ce.SetTable('linkerEvt', lkTb2, 'timeSeries');
        end
        
        function posTb = FindLinkedPositions(ce, uInd)
            % Find peaks in linking scores, and extract peri-peak sequence info for each unit
            % 
            %   posTb = FindLinkedPositions(ce, uInd)
            % 
            % Inputs
            %   ce              NP.CodingExplorer object with sentence PETHs.
            %   uInd            Unit indices.
            % Output
            %   posTb           A table where each row describes a linked position of a unit
            %   posTb.seed      Seed phoneme.
            %   posTb.t0        Seed phoneme time in ce.
            %   posTb.stimId    stimId of the parent sentence.
            %   posTb.stim      NP.TGEvent object of the parent sentence.
            %   posTb.seqStr    A string of text around this position.
            %   posTb.seqTge    A sequence of NP.TGEvent objects time referenced to the seed phoneme.
            %   posTb.score     Linking score.
            %   posTb.clusId    Cluster ID.
            % 
            
            if nargin < 2 || isempty(uInd)
                uInd = (1:ce.numResp)';
            end
            
            [tt, tv, lk, lke] = ce.GetTable('taskTime', 'taskValue', 'linker', 'linkerEvt');
            isStim = tv.phase == "stim";
            tt = tt(isStim,:);
            tv = tv(isStim,:);
            lk = lk(isStim,:);
            lke = lke(isStim,:);
            
            nSen = height(lke);
            nUnit = numel(uInd);
            tbs = cell(nSen, nUnit);
            for i = 1 : nSen
                % Get sentence variables
                phn = lke.time{i};
                tPh = double(phn);
                sen = tt.stim(i);
                wrd = Cut(sen);
                
                % Build a table with all unique sequences
                tb = table;
                tb.seed = erase(phn.GetParentLabel, digitsPattern);
                tb.t0 = tPh;
                tb.stimId(:) = tv.stimId(i);
                tb.stim(:) = tt.stim(i);
                for p = 1 : height(tb)
                    iWrd = find(tPh(p) >= double(wrd), 1, 'last');
                    iWrd = iWrd + (-1:1);
                    iWrd = unique(MMath.Bound(iWrd, [1 numel(wrd)]));
                    tb.seqStr(p) = wrd(iWrd).GetAllParentLabel;
                    tb.seqTge{p} = wrd(iWrd) - tPh(p);
                end
                tb.score = NaN(height(tb), 1);
                tb.isPk = false(height(tb), 1);
                
                for u = 1 : nUnit
                    utb = tb;
                    cid = ce.clusTb.clusId(uInd(u));
                    utb.clusId = repmat(cid, [height(utb) 1]);
                    
                    % Find peaks in raw scores
                    t = lk.time{i};
                    m = t >= tt.stimMatchOn(i) & t <= tt.stimMatchOff(i);
                    [sPk, tPk] = findpeaks(lk.score{i}(m,uInd(u)), t(m));
                    for p = 1 : numel(tPk)
                        [~, I] = min(abs(tPk(p) - tPh));
                        utb.t0(I) = tPk(p);
                        utb.score(I) = sPk(p);
                        utb.isPk(I) = true;
                    end
                    
                    % Remove sequences associated with no peak
                    utb = utb(utb.isPk,:);
                    
                    tbs{i,u} = utb;
                end
            end
            
            posTb = cat(1, tbs{:});
        end
        
        function posTb = IsPosSelective(posTb, clusTb)
            % Determine if the responses at linked positions are sentence-selective
            % 
            %   posTb = IsPosSelective(posTb, clusTb)
            % 
            % Inputs
            %   posTb               Output of FindLinkedPositions, in M11 time.
            %   clusTb              A clusTb with sentence PETH results, in M2 time.
            % Output
            %   posTb.isSelective
            %   posTb.tM2
            % 
            
            posTb.isSelective = false(height(posTb), 1);
            posTb.tM2 = NaN(height(posTb), 2);
            
            for p = 1 : height(posTb)
                % Get sentence response table
                isUnit = clusTb.clusId == posTb.clusId(p);
                if ~any(isUnit)
                    continue
                end
                
                % Find the sentence
                isSen = posTb.stimId(p) == clusTb.stimId{isUnit};
                if ~any(isSen)
                    continue
                end
                t = clusTb.time{isUnit};
                r = clusTb.sigResp{isUnit}(:,isSen);
                
                % Construct interpolants to map time from M11 to M2 stim and prod
                phPk = Cut(Cut(posTb.stim(p)));
                phStim = Cut(Cut(clusTb.stim{isUnit}(isSen)));
                phProd = Cut(Cut(clusTb.prod{isUnit}(isSen)));
                
                [k1, k2] = MLing.FindAlignedTokens(phPk.GetParentLabel, phStim.GetParentLabel);
                gPk2Stim = griddedInterpolant(double(phPk(k1)), double(phStim(k2)));
                
                [k1, k2] = MLing.FindAlignedTokens(phPk.GetParentLabel, phProd.GetParentLabel);
                gPk2Prod = griddedInterpolant(double(phPk(k1)), double(phProd(k2)));
                
                % Check response significance
                tTo = [-0.05 0.05];
                iTo = round(tTo/diff(t(1:2)));
                iTo = iTo(1) : iTo(2);
                
                switch clusTb.hcGroup(isUnit)
                    case "mirror"
                        dtKer = 0.14;
                    case "bridge"
                        dtKer = 0.14;
                    otherwise
                        dtKer = 0;
                end
                
                tStim = gPk2Stim(posTb.t0(p));
                [~, I] = min(abs(t - tStim + dtKer));
                ind = MMath.Bound(I+iTo, [1 numel(r)]);
                isSigStim = any(~isnan(r(ind)));
                
                tProd = gPk2Prod(posTb.t0(p));
                [~, I] = min(abs(t - tProd));
                ind = MMath.Bound(I+iTo, [1 numel(r)]);
                isSigProd = any(~isnan(r(ind)));
                
                posTb.isSelective(p) = isSigStim || isSigProd;
                posTb.tM2(p,:) = [tStim tProd];
            end
        end
        
        % Unit profile
        function PlotUnitProfile(ce, pkTb, triTb, clusId)
            % Plot sentence PETH and score overlay, peak scatter, stRF
            % 
            %   PlotUnitProfile(ce, pkTb, triTb, clusId)
            % 
            
            nRows = 6;
            colDist = [2 2 1 1 1 1];
            tl = tiledlayout(nRows, sum(colDist));
            tl.Padding = "tight";
            
            
            % Rasters
            ntArgs = MPlot.FindTileInd(nRows, colDist, 1:nRows, 1);
            ax = nexttile(ntArgs{:});
            
            sUnit = NP.Unit.LoadUnitCache(clusId, 'DataSource', 'm2');
            se = LMV.SE.UnitCache2SE(sUnit, 'StimIdList', LMV.Param.stimIdList14);
            tt = se.GetTable("taskTime");
            tWin = [tt.stimOn(1) tt.prodMatchOff(1)] + [-.3 .3];
            
            NP.UnitPlot.Heatmap(ax, se, [], tWin, 1);
            NP.TaskBaseClass.PlotBlockBoundary(ax, se, [], tWin, "stimText", 'Label', true, 'Color', [0 0 0 .3]);
            NP.TaskBaseClass.PlotEventBoundary(ax, se, [], tWin, ...
                ["cue1On", "stimOn", "stimOff", "cue3On", "prodMatchOn", "prodMatchOff"], ...
                @(x) LMV.Param.GetTaskPhaseColors(["atten", "stim", "stim", "init", "prod", "prod"]));
            
            
            % Timeseries overaly
            ce = ce.Duplicate;
            tv = ce.GetTable("taskValue");
            isList = ismember(tv.stimId, LMV.Param.stimIdList14);
            ce.RemoveEpochs(~isList);
            
            tv = ce.GetTable("taskValue");
            [~, I] = MMath.SortLike(tv.stimId, flip(LMV.Param.stimIdList14), false);
            I = [I; I+numel(I)];
            ce.SortEpochs(I);
            
            ntArgs = MPlot.FindTileInd(nRows, colDist, 1:nRows, 2);
            ax = nexttile(ntArgs{:});
            LMV.Linker.LM.PlotSentenceOverlay(ce, clusId, 'StimIdList', LMV.Param.stimIdList14, 'TraceList', ["PETH", "pred", "score"], 'Parent', ax);
            
            
            % Score leaderboard
            ntArgs = MPlot.FindTileInd(nRows, colDist, 1:2, 3);
            ax = nexttile(ntArgs{:});
            LMV.Linker.LM.PlotLeaderboard(pkTb, clusId, ax);
            
            % Sequence overlay
            pkTb = pkTb(pkTb.clusId==clusId, :);
            pkTb = sortrows(pkTb, 'score', 'descend');
            
            nOverlay = nRows * 2;
            pkTb = pkTb(1:min(end,nOverlay), :);
            for i = 1 : height(pkTb)
                isSeed = pkTb.seed(i) == triTb.Properties.RowNames;
                triRow = triTb(isSeed,:);
                triRow.Properties.RowNames = {};
                for j = 1 : width(triRow)
                    seqTb = triRow.(j){1}.seqTb;
                    isSeq = seqTb.seqStr == pkTb.seqStr(i);
                    if ~any(isSeq)
                        triRow = [];
                        break
                    end
                    triRow.(j){1}.seqTb = seqTb(isSeq,:);
                end
                pkTb.triRow{i} = triRow;
            end
            triTb = cat(1, pkTb.triRow{:});
            
            axObjs = cell(height(triTb), 1);
            for i = 1 : height(pkTb)
                r = 1 + mod(i-1, nRows);
                c = 3 + ceil(i/nRows);
                ntArgs = MPlot.FindTileInd(nRows, colDist, r, c);
                axObjs{i} = nexttile(ntArgs{:});
            end
            axObjs = cat(1, axObjs{:});
            
            uIdx = find(ce.clusTb.clusId == clusId, 1);
            LMV.Linker.PlotPeriEventSeqArray(triTb, 'UnitIdx', uIdx, 'Parent', axObjs);
            
            
            % Linker kernel
            lm = ce.clusTb.lm{uIdx};
            if isempty(lm)
                return
            end
            ntArgs = MPlot.FindTileInd(nRows, colDist, 3, 3);
            ax = nexttile(ntArgs{:});
            LMV.Linker.LM.PlotWeights(lm, 'Parent', ax);
            
            
            % Sliding time RF models
            stRF = ce.clusTb.phone_stim{uIdx};
            if isempty(stRF)
                return
            end
            RF = stRF.mdls{stRF.r2Idx};
            RF.feats = stRF.feats;
            RF.resps = stRF.resps;
            
            ntArgs = MPlot.FindTileInd(nRows, colDist, 1:3, 6);
            ax = nexttile(ntArgs{:});
            LMV.TRF.PlotWeights2(stRF, 'Parent', ax);
            
            ntArgs = MPlot.FindTileInd(nRows, colDist, 4:6, 6);
            ax = nexttile(ntArgs{:});
            LMV.RF.PlotWeights(RF, 'Orientation', "vertical", 'Parent', ax);
        end
        
        function PlotSentenceOverlay(ce, varargin)
            % Plot overlays of aligned stim and prod responses
            % 
            %   PlotSentenceOverlay(ce)
            %   PlotSentenceOverlay(ce, clusId)
            %   PlotSentenceOverlay(ce, ..., 'StimIdList', LMV.Param.stimIdList4)
            % 
            
            p = inputParser;
            p.addOptional('clusId', ce.clusTb.clusId, @isnumeric);
            p.addParameter('StimIdList', LMV.Param.stimIdList4, @(x) isstring(x) || iscellstr(x) || ischar(x));
            p.addParameter('TraceList', 'peth', @(x) isstring(x));
            p.addParameter('Parent', [], @(x) true);
            p.parse(varargin{:});
            clusId = p.Results.clusId;
            stimIdList = p.Results.StimIdList;
            traceList = lower(string(p.Results.TraceList));
            axObjs = p.Results.Parent;
            
            [~, uInd] = MMath.SortLike(ce.clusTb.clusId, clusId);
            
            [tt, tv, lk] = ce.GetTable('taskTime', 'taskValue', 'linker');
            clusId = ce.clusTb.clusId;
            cc = lines;
            
            for u = uInd(:)'
                if isnan(u)
                    continue
                end
                if isempty(axObjs)
                    ax = nexttile;
                else
                    ax = axObjs;
                end
                
                % Prepare PETHs
                m = tv.phase == "stim" & ismember(tv.stimId, stimIdList);
                [ts, rs] = ce.GetArray('resp', m, u+1, 'DimCat', 0);
                
                m = tv.phase == "prod" & ismember(tv.stimId, stimIdList);
                [tp, rp] = ce.GetArray('resp', m, u+1, 'DimCat', 0);
                
                rr = cell2mat([rs; rp]);
                rrMax = max(rr);
                if rrMax > 0
                    rs = cellfun(@(x) x/rrMax/2, rs, 'Uni', false);
                    rp = cellfun(@(x) x/rrMax/2, rp, 'Uni', false);
                end
                
                m = tv.phase == "stim" & ismember(tv.stimId, stimIdList);
                y = (1 : numel(rs)) - 0.4;
                
                % PETHs
                if any(traceList=="peth")
                    MPlot.PlotTraceLadder(ts, rs, y, 'Color', cc(1,:)); hold on
                    MPlot.PlotTraceLadder(tp, rp, y, 'Color', cc(2,:));
                end
                
                if any(traceList=="pred")
                    t = cellfun(@double, lk.time(m), 'Uni', false);
                    pred = cellfun(@(x) x(:,u)/2, lk.pred(m), 'Uni', false);
                    MPlot.PlotTraceLadder(t, pred, y, 'Color', cc(3,:)); hold on
                end
                
                % Model score timeseries
                if any(traceList=="score")
                    t = cellfun(@double, lk.time(m), 'Uni', false);
                    score = cellfun(@(x) x(:,u), lk.score(m), 'Uni', false);
                    MPlot.PlotTraceLadder(t, score, y, 'Color', cc(4,:)); hold on
                end
                
                % Phonetic markers
                m = tv.phase == "stim" & ismember(tv.stimId, stimIdList);
                m = find(m);
                for i = 1 : numel(m)
                    tge = tt.stim(m(i));
                    tge.PlotTiers(y(i)+0.8, -0.2, 0, 'TierInd', [1 3], 'Parent', ax);
                end
                
                ax.YDir = 'normal';
                ax.Box = 'off';
                ax.YTick = 1:numel(rs);
                ax.XLim = [-.4 max(cellfun(@(x)x(end), ts))];
                ax.YLim = [0.5 numel(rs)+.5];
                ax.Title.String = sprintf("u%i", clusId(u));
                MPlot.Axes(ax);
            end
        end
        
        function PlotLeaderboard(pkTb, clusId, axObjs)
            % 
            %   PlotLeaderboard(pkTb, clusId, seqInd)
            % 
            
            for u = 1 : numel(clusId)
                % Create axes
                if exist('axObjs', 'var')
                    ax = axObjs(u);
                end
                
                % Find the unit index
                m = pkTb.clusId == clusId(u);
                if ~any(m)
                    ax.Visible = 'off';
                    continue
                end
                
                % Get plotting variables
                sc = pkTb.score(m,:);
                x = zeros(size(sc));
                % lb = cellfun(@(x) x(find(double(x)<=0, 1, 'last')).GetParentLabel, pkTb.seqTge(m));
                lb = pkTb.seqStr(m);
                
                % Plot
                text(ax, x, sc, lb, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
                ax.XLim = [-0.1 1.1];
                ax.YLim = [-0.1 1.1];
                ax.YLabel.String = "Score";
                ax.Title.String = "u"+clusId(u);
                MPlot.Axes(ax);
            end
            
            b = brush(gcf);
            b.Enable = 'on';
            b.ActionPostCallback = @LMV.Linker.PlotSeqPethOverlayOnBrush;
        end
        
        function PlotWeights(mdl, varargin)
            % Plot weight vectors for models fitted by lasso
            % 
            %   PlotWeights(mdl)
            %   PlotWeights(mdl, ..., 'ClusInfo', [])
            %   PlotWeights(mdl, ..., 'Orientation', 'horizontal')
            %   PlotWeights(mdl, ..., 'Parent', gca)
            % 
            
            p = inputParser;
            p.addParameter('ClusInfo', [], @(x) istable(x) || isstruct(x));
            p.addParameter('Orientation', "horizontal", @(x) any(strcmpi(x, ["horizonal", "vertical"])));
            p.addParameter('Parent', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
            clusInfo = p.Results.ClusInfo;
            ori = lower(p.Results.Orientation);
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
            
            % Make weight vector
            b = flip(mdl.Beta);
            x = mdl.dt;
            
            % Get error
            sig = mdl.null.BetaPval < 0.05;
            sig = flip(sig);
            
            % Plot
            if ori == "horizontal"
                plot(ax, [0.5 numel(b)+0.5], [0 0]', 'Color', [0 0 0 .5]); hold(ax, 'on')
                plot(ax, x, b, '-', 'Color', [0 0 0]);
                plot(ax, x(sig), b(sig), '.', 'Color', [0 0 1]);
            else
                plot(ax, [0 0]', [0.5 numel(b)+0.5], 'Color', [0 0 0 .5]); hold(ax, 'on')
                
                plot(ax, b(~sig), x(~sig), '.', 'Color', [0 0 0]);
                plot(ax, [b(~sig)*0 b(~sig)]', [x(~sig) x(~sig)]', '-', 'Color', [0 0 0]);
                
                plot(ax, b(sig), x(sig), 'o', 'Color', [0 0 1]);
                plot(ax, [b(sig)*0 b(sig)]', [x(sig) x(sig)]', '-', 'Color', [0 0 1]);
            end
            
            if ~isempty(clusInfo)
                ax.Title.String = sprintf("u%i, %ium", clusInfo.clusId, clusInfo.depth);
            else
                ax.Title.String = sprintf("%s, r2=%.2f (p=%.2f), L=%.2f", mdl.resps, mdl.r2, mdl.null.r2Pval, mdl.mdl.Lambda);
                % ax.Title.String = sprintf("%s, r2=%.2f (p=%.2f)", mdl.resps, mdl.r2, mdl.null.r2Pval);
            end
            
            if ori == "horizontal"
                ax.YLabel.String = "Weight";
                ax.XLim = x([1 end]);
                ax.XLabel.String = "Kernel time (s)";
            else
                ax.XLabel.String = "Weight";
                ax.YLim = [0 numel(b)+1];
                ax.YTick = 1 : numel(b);
                ax.YTickLabel = mdl.feats;
                ax.YDir = 'reverse';
            end
            
            ax.Box = 'off';
        end
        
        % Clustering
        function PlotDendrogram(sHC, cTb, colorBy, varargin)
            % Make an interactive dendrogram of kernels
            % 
            %   PlotDendrogram(sHC, cTb, colorBy, varargin)
            % 
            
            p = inputParser;
            p.addRequired('ColorBy', @(x) isstring(x) || ischar(x));
            p.addParameter('MarkUnit', [], @(x) iscell(x));
            p.addParameter('MarkColor', [], @(x) isnumeric(x));
            p.addParameter('Ribbon', [], @(x) ismember(x, cTb.Properties.VariableNames));
            p.parse(colorBy, varargin{:});
            cidMark = p.Results.MarkUnit;
            ccMark = p.Results.MarkColor;
            ribVar = p.Results.Ribbon;
            
            if ~isempty(cidMark) && isempty(ccMark)
                ccMark = [0 0 0];
            end
            if size(ccMark,1) < numel(cidMark)
                ccMark = repmat(ccMark, numel(cidMark), 1);
            end
            
            % Plot dendrogram
            ax = nexttile;
            hh = dendrogram(sHC.tree, Inf, 'Reorder', sHC.leafOrder); hold on
            for i = 1 : numel(hh)
                hh(i).Color = [0 0 0];
            end
            y0 = ax.YLim(1);
            
            % Plot dots for brushing
            x = 1:height(cTb);
            y = x*0 + y0;
            h = plot(x, y, 'w.');
            h.UserData.forBrush = true;
            
            % Plot leaves by groups
            ccFun = @lines; % default color scheme
            switch colorBy
                case 'waveform'
                    val = cTb.wfId;
                    cc = ccFun(numel(unique(val)));
                case 'isi'
                    val = cTb.isiKmId;
                    colorBy = 'burstiness';
                    cc = ccFun(numel(unique(val)));
                case 'depth'
                    val = MMath.Bound(cTb.depth, 500:1e3:4500);
                    ccFun = @(n) flip(parula(n), 1);
                    cc = ccFun(numel(unique(val)));
                case 'recording'
                    val = cTb.recId;
                    cc = ccFun(numel(unique(val)));
                case 'region'
                    val = cTb.region;
                    val = categorical(val, {'mPrCG', 'vPrCG', 'IFG', 'STG'});
                    cc = LMV.Param.GetRegionColors(unique(val));
                otherwise
                    fprintf("'%s' is not a supported color group option.\n", colorBy);
                    val = cTb.(colorBy);
                    [~, ~, val] = histcounts(val, 5);
                    cc = ccFun(numel(unique(val)));
            end
            
            if ~exist('groups', 'var')
                groups = unique(val);
            end
            for i = numel(groups) : -1 : 1
                m = groups(i) == val;
                if ~any(m)
                    continue
                end
                x = find(m);
                y = x*0 + y0;
                hh(i) = plot(x, y, '.', 'MarkerSize', 16, 'Color', cc(i,:));
            end
            % lgd = legend(hh, string(groups), 'Location', 'eastoutside');
            % lgd.Title.String = colorBy;
            % lgd.Interpreter = 'none';
            
            % Indicate manually selected units
            for i = 1 : numel(cidMark)
                cid = cidMark{i};
                m = ismember(cTb.clusId, cid);
                x = find(m);
                y = x*0 + y0;
                plot(x, y, 'o', 'MarkerSize', 12, 'Color', ccMark(i,:));
            end
            
            % Plot group ribbon
            if ~isempty(ribVar)
                hRib = diff(ax.YLim)*0.1;
                lb = cTb.(ribVar);
                if numel(unique(lb))==4
                    cRib = [0 1 0; 1 0 0; 0 0 1; 1 1 1];
                else
                    cRib = @lines;
                end
                [xc, yc] = MPlot.GroupRibbon(lb, [-2 -1]*hRib, cRib, 'Groups', unique(lb));
                text(xc, yc-hRib, string(unique(lb)));
                ax.YLim(1) = -3*hRib;
            end
            
            % Plot kernel heatmap
            M = cellfun(@(x) flip(x.Beta), cTb.linker, 'Uni', false);
            M = cat(2, M{:});
            M = zscore(M, 1, 1);
            imagesc(1:size(M,2), linspace(-3, -.5, size(M,1)), M, [-2 4]);
            
            ax.YLim(1) = -3;
            % ax.XTick = 1 : height(cTb);
            % ax.XTickLabel = "u"+cTb.clusId;
            ax.XTick = [];
            ax.YTick = [];
            ax.XAxis.Visible = 'off';
            ax.YAxis.Visible = 'off';
            MPlot.Axes(ax);
            
            % Initialize brush callback
            ax.UserData.clusTb = cTb;
            brushObj = brush(gcf);
            brushObj.ActionPostCallback = @LMV.Linker.LM.PlotBrushedKernels;
            brushObj.Enable = 'on';
        end
        
        function FigExampleKernel(mdl, varargin)
            % Plot weight vectors for models fitted by lasso
            % 
            %   FigExampleKernel(mdl)
            %   FigExampleKernel(mdl, ..., 'ClusInfo', [])
            %   FigExampleKernel(mdl, ..., 'Parent', gca)
            % 
            
            p = inputParser;
            p.addParameter('Parent', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
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
            
            % Make weight vector
            b = mdl.Beta;
            % b = flip(b);
            x = mdl.dt;
            
            % Get error
            sig = mdl.null.BetaPval < 0.05;
            sig = flip(sig);
            
            % Plot
            plot(ax, [0.5 numel(b)+0.5], [0 0]', 'Color', [0 0 0 .5]); hold(ax, 'on')
            plot(ax, x, b, '-', 'Color', [0 0 0]);
            plot(ax, x(sig), b(sig), '.', 'Color', [0 0 1]);
            
            ax.Title.String = sprintf("%s", mdl.resps);
            ax.YLabel.String = "Weight";
            ax.XLim = x([1 end]);
            ax.XLabel.String = "Kernel time (s)";
            ax.Box = 'off';
        end
        
        function PlotKernelOverlay(mdls, varargin)
            % Plot weight vectors for models fitted by lasso
            % 
            %   PlotKernelOverlay(mdls)
            %   PlotKernelOverlay(mdls, ..., 'ClusInfo', [])
            %   PlotKernelOverlay(mdls, ..., 'Orientation', 'horizontal')
            %   PlotKernelOverlay(mdls, ..., 'Parent', gca)
            % 
            
            p = inputParser;
            p.addParameter('ClusInfo', [], @(x) istable(x) || isstruct(x));
            p.addParameter('Parent', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
            clusInfo = p.Results.ClusInfo;
            ax = p.Results.Parent;
            
            if isempty(ax)
                ax = gca;
            end
            
            if isempty(mdls)
                axis(ax, 'off');
                return
            end
            if ~iscell(mdls)
                mdls = {mdls};
            end
            
            if istable(clusInfo)
                clusInfo = table2struct(clusInfo);
            end
            
            % Make weight vector
            bb = cell2mat(cellfun(@(x) x.Beta, mdls(:)', 'Uni', false));
            % bb = flip(bb);
            bb = MMath.Normalize(bb, 'zscore');
            x = mdls{1}.dt';
            
            % Plot
            % plot(ax, x([1 end]), [0 0]', 'Color', [0 0 0 .5]); 
            plot(ax, x, bb, '-', 'Color', [0 0 0 .2]); hold(ax, 'on')
            plot(ax, x, mean(bb,2), '-', 'Color', [0 0 0], 'LineWidth', 1.5);
            
            % if ~isempty(clusInfo)
            %     ax.Title.String = sprintf("u%i, %ium", clusInfo.clusId, clusInfo.depth);
            % else
            %     ax.Title.String = sprintf("%s, r2=%.2f (p=%.2f), L=%.2f", mdl.resps, mdl.r2, mdl.null.r2Pval, mdl.mdl.Lambda);
            % end
            
            ax.YLabel.String = "Z-scored weight";
            ax.XLim = x([1 end]);
            ax.XLabel.String = "Kernel time (s)";
            ax.Box = 'off';
        end
        
        function PlotBrushedKernels(fig, axesStruct)
            % 
            %   PlotBrushedKernels(fig, axesStruct)
            % 
            
            % Find brushed units
            clusTb = axesStruct.Axes.UserData.clusTb;
            hh = axesStruct.Axes.Children;
            b = [];
            for i = 1 : numel(hh)
                % Check if handle is the plot for brushing
                s = hh(i).UserData;
                if ~isfield(s, 'forBrush') || ~s.forBrush
                    continue
                end
                
                % Check if brushed any data
                b = logical(hh(i).BrushData);
            end
            if ~any(b)
                return
            end
            
            % Create figure if absent
            if isempty(fig.UserData) || ~isvalid(fig.UserData.cbFig)
                fig.UserData.cbFig = MPlot.Figure( ...
                    'Name', 'Selected Units', ...
                    'NumberTitle', 'off', ...
                    'Menubar', 'none', ...
                    'Toolbar', 'figure', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
            end
            
            mdls = clusTb.linker(b);
            clf(fig.UserData.cbFig);
            tl = tiledlayout("flow", "Parent", fig.UserData.cbFig);
            for i = 1 : numel(mdls)
                ax = nexttile(tl);
                LMV.Linker.LM.PlotWeights(mdls{i}, 'Parent', ax);
            end
            
            f = MPlot.Figure(555); clf(f);
            LMV.Overview.SessionFromCache(clusTb.clusId(b), 'DataSource', 'm2', 'TaskPhase', 'lmv');
        end
        
        
        
        function PlotKernelScatter(cTb, colorBy, varargin)
            % Make an interactive scatter plot of kernels
            % 
            %   PlotKernelScatter(ce, colorBy)
            % 
            
            p = inputParser;
            p.addRequired('ColorBy', @(x) isstring(x) || ischar(x));
            p.addParameter('MarkUnit', [], @(x) iscell(x));
            p.addParameter('MarkColor', [], @(x) isnumeric(x));
            p.parse(colorBy, varargin{:});
            cidMark = p.Results.MarkUnit;
            ccMark = p.Results.MarkColor;
            
            if ~isempty(cidMark) && isempty(ccMark)
                ccMark = [0 0 0];
            end
            if size(ccMark,1) < numel(cidMark)
                ccMark = repmat(ccMark, numel(cidMark), 1);
            end
            
            ax = nexttile;
            
            % Plot dots for brushing
            x = cTb.embedCoords(:,1);
            y = cTb.embedCoords(:,2);
            h = plot(x, y, 'w.'); hold on
            h.UserData.forBrush = true;
            
            % Plot overlays
            ccFun = @lines; % default color scheme
            switch colorBy
                case 'waveform'
                    val = cTb.wfId;
                    cc = ccFun(numel(unique(val)));
                case 'isi'
                    val = cTb.isiKmId;
                    colorBy = 'burstiness';
                    cc = ccFun(numel(unique(val)));
                case 'depth'
                    val = MMath.Bound(cTb.depth, 500:1e3:4500);
                    ccFun = @(n) flip(parula(n), 1);
                    cc = ccFun(numel(unique(val)));
                case 'recording'
                    val = cTb.recId;
                    cc = ccFun(numel(unique(val)));
                case 'region'
                    val = cTb.region;
                    val = categorical(val, {'mPrCG', 'vPrCG', 'IFG', 'STG'});
                    cc = LMV.Param.GetRegionColors(unique(val));
                otherwise
                    fprintf("'%s' is not a supported color group option.\n", colorBy);
                    val = cTb.(colorBy);
                    [~, ~, val] = histcounts(val, 5);
                    cc = ccFun(numel(unique(val)));
            end
            
            % Plot scatter by groups
            if ~exist('groups', 'var')
                groups = unique(val);
            end
            for i = numel(groups) : -1 : 1
                m = groups(i) == val;
                if ~any(m)
                    continue
                end
                x = cTb.embedCoords(m,1);
                y = cTb.embedCoords(m,2);
                hh(i) = plot(x, y, '.', 'MarkerSize', 16, 'Color', cc(i,:));
            end
            
            % Indicate manually selected units
            for i = 1 : numel(cidMark)
                cid = cidMark{i};
                m = ismember(cTb.clusId, cid);
                x = cTb.embedCoords(m,1);
                y = cTb.embedCoords(m,2);
                plot(x, y, 'o', 'MarkerSize', 12, 'Color', ccMark(i,:));
            end
            
            % Format plot
            lgd = legend(hh, string(groups), 'Location', 'eastoutside');
            lgd.Title.String = colorBy;
            lgd.Interpreter = 'none';
            
            x = cTb.embedCoords(:,1);
            y = cTb.embedCoords(:,2);
            ax.XLim = [min(x) max(x)] + 0.05*[-1 1]*(max(x) - min(x));
            ax.YLim = [min(y) max(y)] + 0.05*[-1 1]*(max(y) - min(y));
            ax.XAxis.Visible = 'off';
            ax.YAxis.Visible = 'off';
            MPlot.Axes(ax);
            
            % Initialize brush callback
            ax.UserData.clusTb = cTb;
            brushObj = brush(gcf);
            brushObj.ActionPostCallback = @LMV.Linker.LM.PlotBrushedKernels;
            brushObj.Enable = 'on';
        end
        
        function PlotLinkPositions(ax, posTb, varargin)
            % 
            %   PlotLinkPositions(ax, posTb)
            %   PlotLinkPositions(ax, posTb, ..., "StimIdList", LMV.Param.stimIdList14)
            % 
            
            p = inputParser;
            p.addParameter("StimIdList", LMV.Param.stimIdList14)
            p.parse(varargin{:});
            stimIdList = p.Results.StimIdList;
            
            stimIdList(~ismember(stimIdList, posTb.stimId)) = [];
            senInd = arrayfun(@(x) find(posTb.stimId==x,1), stimIdList);
            sen = posTb.stim(senInd);
            senDur = sen.GetDuration;
            [~, I] = sort(senDur);
            
            cid = unique(posTb.clusId);
            cc = lines(numel(cid));
            mk = repelem(["o", "^", "x", "square", "diamond"], 7);
            senLb = string();
            for i = 1 : numel(stimIdList)
                wrd = Cut(posTb.stim(senInd(i)));
                senLb(i) = wrd(1:2).GetAllParentLabel + " ...";
            end
            
            nSen = numel(stimIdList);
            for i = 1 : nSen
                wrd = Cut(sen(I(i)));
                ph = Cut(wrd);
                
                sen(I(i)).Plot(i+[-.4 .4], FontSize=0);
                % wrd.Plot(i+[-.4 .4], 'Color', [0 0 0]+.5);
                ph.Plot(i, 'Color', [0 0 0]+.9);
                
                isSen = posTb.stimId == stimIdList(I(i));
                for u = 1 : numel(cid)
                    isUnit = posTb.clusId == cid(u);
                    isPk = isSen & isUnit;
                    tPk = posTb.t0(isPk);
                    % plot(tPk, tPk*0+i, 'Color', cc(u,:), LineStyle="none", Marker=mk(u), MarkerSize=8); hold on
                    plot(tPk, tPk*0+i, 'Color', 'k', LineStyle="none", Marker='x', MarkerSize=8); hold on
                end
            end
            ax.XLabel.String = "Aligned prod. time (s)";
            ax.YTick = 1 : nSen;
            ax.YTickLabel = senLb(I);
            ax.YLim = [0 nSen+1];
            ax.YDir = "reverse";
            MPlot.Axes(ax);
        end
        
        function PlotLinkPositionHist(ax, pkTb, varargin)
            % 
            %   PlotLinkPositionHist(ax, pkTb, varargin)
            % 
            
            p = inputParser;
            p.addParameter("StimIdList", LMV.Param.stimIdList14)
            p.parse(varargin{:});
            stimIdList = p.Results.StimIdList;
            
            isStim = ismember(pkTb.stimId, stimIdList);
            tPk = pkTb.t0(isStim) ./ pkTb.stim(isStim).GetDuration;
            tEdges = 0:0.1:1;
            N = histcounts(tPk, tEdges);
            
            b = bar(MMath.BinEdges2Centers(tEdges), N, 'hist');
            b.FaceColor = "none";
            ax.YLabel.String = "# of units";
            ax.XLabel.String = "Relative sentence position (frac.)";
            MPlot.Axes(ax);
        end
        
        function p = TestUniformDist(pkTb, varargin)
            % 
            %   TestUniformDist(pkTb, varargin)
            % 
            
            p = inputParser;
            p.addParameter("StimIdList", LMV.Param.stimIdList14)
            p.parse(varargin{:});
            stimIdList = p.Results.StimIdList;
            
            isStim = ismember(pkTb.stimId, stimIdList);
            tPk = pkTb.t0(isStim) ./ pkTb.stim(isStim).GetDuration;
            tNull = linspace(0, 1, numel(tPk));
            [~, p] = kstest2(tPk, tNull);
        end
        
        function PlotNumLinkHist(ax, pkTb, itemName, varargin)
            % 
            %   PlotNumLinkHist()
            % 
            
            p = inputParser;
            p.addParameter("StimIdList", LMV.Param.stimIdList14)
            p.parse(varargin{:});
            stimIdList = p.Results.StimIdList;
            
            cid = unique(pkTb.clusId);
            sp = NaN(size(cid));
            for i = 1 : numel(cid)
                isUnit = pkTb.clusId==cid(i);
                switch itemName
                    case 'position'
                        sp(i) = sum(isUnit);
                        x = 1:14;
                    case 'sentence'
                        sp(i) = numel(unique(pkTb.stimId(isUnit)));
                        x = 1:14;
                    otherwise
                        error("'%s' is not a valid item name.", itemName);
                end
            end
            N = histcounts(categorical(sp), categorical(x));
            
            b = bar(x, N, 'hist');
            b.FaceColor = "none";
            ax.XTick = 1:2:x(end);
            ax.YLabel.String = "# of units";
            ax.XLabel.String = sprintf("# of %ss", itemName);
            MPlot.Axes(ax);
        end
    end
    
end