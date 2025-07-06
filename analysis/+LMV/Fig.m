classdef Fig
    
    methods(Static)
        function varargout = GetExampleInfo(figName)
            % Return handpicked cluster IDs for a given analysis
            % 
            %   info = GetExampleInfo(figName)
            % 
            % Input
            %   figName     Name of a figure or final panel.
            % Output
            %   uid         Long, globally unique cluster IDs including the recording and probe digits.
            % 
            
            switch lower(figName)
                case 'fig1g'
                    % uid = [ ...
                    %     410100245 450200307 460100248 460100030 460100005 450100286, ...  % phasic units. 
                    %     410100155 410100441 440300575 410100450 450100323 410100463, ...  % non-phasic units. 
                    %     ];
                    
                    uu = cell(0);
                    
                    % Stim
                    uu{end+1} = struct('clusId', 410100249, 'startText', "it was nobody");
                    
                    % Delay
                    % uu{end+1} = struct('clusId', 560110633, 'startText', "the girl");
                    % uu{end+1} = struct('clusId', 410100138, 'startText', "were you in love");
                    % uu{end+1} = struct('clusId', 410100408, 'startText', "will you tell me");
                    
                    % Prod
                    uu{end+1} = struct('clusId', 410100245, 'startText', "they even");
                    
                    % Init
                    uu{end+1} = struct('clusId', 450200287, 'startText', "we've got");
                    
                    % % Stim-delay
                    % uu{end+1} = struct('clusId', 410100454, 'startText', "will you tell me");
                    
                    % % Init-prod
                    % uu{end+1} = struct('clusId', 410100115, 'startText', "you took me");
                    
                    % Stim-prod
                    uu{end+1} = struct('clusId', 460100003, 'startText', "junior");
                    
                    % Stim-delay-init-prod
                    % uu{end+1} = struct('clusId', 410100454, 'startText', "something pulled");
                    % uu{end+1} = struct('clusId', 450100286, 'startText', "you took me");
                    uu{end+1} = struct('clusId', 410100155, 'startText', "have you got");
                    
                    uu = struct2table([uu{:}]);
                    varargout{1} = uu;
                    
                case 'selectivity'
%                     uid = [ ...
%                         460100240, ... % 0-1mm
%                         460100098 440200462 450100286 460100013 460100005, ... % 1-2mm
%                         410100249 410100260 440200479 460100003 410100245 460100248, ... % 2-3mm
%                         440300575 410100211 440300585 440300575 410100367 450200307, ... % 3-4mm
%                         410100463 410100441 410100170 450100323 440300587, ... % 4-5mm
%                         410100441 410100454 440300601 410100434 470300139, ... % 5-6mm
%                         ];
                    uid = [ ...
                        410100245 450200307 460100248 460100030 460100005 450100286, ...  % phasic units. 
                        410100155 410100441 440300575 410100450 450100323 410100463, ...  % non-phasic units. 
                        ];
                    
                otherwise
                    varargout{1} = [];
            end
        end
        
        % Figure 1
        function UnitRaster(ax, unitId, varargin)
            % Plot spike raster for a single unit
            % 
            %   UnitRaster(unitId)
            %   UnitRaster(unitId, ...)
            % 
            % See also LMV.Overview.Sentences
            
            % Process inputs
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('StimIdList', LMV.Param.stimIdList4, @isstring);
            p.parse(varargin{:});
            stimIdList = p.Results.StimIdList;
            
            % Load cache
            sUnit = NP.Unit.LoadUnitCache(unitId, varargin{:});
            
            % Put data into se
            se = LMV.SE.UnitCache2SE(sUnit, varargin{:});
            
            % Add unit attributes to clusTb
            clusTb = NP.Unit.GetClusTb(se);
            se.userData.ksMeta.clusTb = clusTb;
            
            % Make sentence table
            senTb = se.SplitConditions('stimId', 'taskValue');
            [~, ind] = MMath.SortLike(senTb.stimId, stimIdList);
            senTb = senTb(ind,:);
            
            nSen = height(senTb);
            for i = 1 : nSen
                se = senTb.se(i);
                [tt, tv] = se.GetTable('taskTime', 'taskValue');
                tWin = [tt.cue1On(1) tt.matchOff(1)] + [-1 1]*.2;
                
                % Raster
                NP.UnitPlot.RasterStack(ax, se, [], tWin, 1, 'PETH', false, 'PeakSpikeRates', [], 'LineWidth', 1.5);
                
                % 
                tge = [tt.stim(1), tt.prod{1}];
                % tge = Cut(tge);
                y0 = 0;
                dy = .3;
                tge.PlotTiers(y0, dy, 0, 'TierInd', 2, 'FontSize', 10, 'Parent', ax);
                
                ax.XGrid = 'off';
                ax.XLabel.String = [];
                ax.XTick = [];
                ax.YLim = [y0 1.5];
                % ax.YLabel.String = "Repeats";
            end
        end
        
        function SessionRaster(ax, uid, varargin)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('PlotFun', 'raster', @(x) ischar(x) || isa(x, "function_handle"));
            p.parse(varargin{:});
            fPlot = p.Results.PlotFun;
            
            % Load cache
            sUnit = NP.Unit.LoadUnitCache(uid, varargin{:});
            
            % Put data into se
            se = LMV.SE.UnitCache2SE(sUnit, varargin{:});
            
            % 
            tWin = [-0.3 6];
            
            NP.UnitPlot.Heatmap(ax, se, [], tWin, 1);
            
            NP.TaskBaseClass.PlotEventBoundary(ax, se, [], tWin, ["cue1On", "stimOn", "stimOff", "cue3On", "prodMatchOn", "prodMatchOff"], ...
                @(x) LMV.Param.GetTaskPhaseColors(["none", "stim", "stim", "none", "prod", "prod"]));
            
            NP.TaskBaseClass.PlotBlockBoundary(ax, se, [], tWin, "stimText", 'Label', false, 'Color', [0 0 0 .3]);
            
            grid(ax, 'off');
            ax.TickLength = [.01 .01];
            ax.YTick = [];
            ax.XLabel.String = [];
            ax.XTickLabel = [];
            MPlot.Axes(ax);
        end
        
        function SessionFromCache(uid, varargin)
            % Plot session heatmap rasters
            % 
            %   SessionFromCache(uid)
            % 
            % Inputs
            %   unitId          A vector of cluster IDs. IDs should be consistent with the ID format (i.e. origianl 
            %                   KS short ID or unique long ID) in se.userData.ksMeta.clusTb.clusId. NaN IDs are 
            %                   legal and will be ignored.
            %   
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('NRows', 1, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('PlotFun', 'raster', @(x) ischar(x) || isa(x, "function_handle"));
            p.parse(varargin{:});
            nr = p.Results.NRows;
            fPlot = p.Results.PlotFun;
            
            uid = unique(uid, 'stable');
            uInd = 1:numel(uid);
            
            % Load cache
            sUnit = NP.Unit.LoadUnitCache(uid, varargin{:});
            
            % Put data into se
            se = LMV.SE.UnitCache2SE(sUnit, varargin{:});
            
            % Configure plots
            upp = numel(uid);
            ppr = ceil(upp/nr) + 1; % plot per row = #unit / 3 rows + one label plot
            panels1 = repmat({fPlot}, [nr ppr]);
            panels1args = repmat({{}}, [nr ppr]);
            
            ebArgs = { ...
                ["cue1On", "stimOn", "stimOff", "cue3On", "prodMatchOn", "prodMatchOff"], ...
                @(x) LMV.Param.GetTaskPhaseColors(["none", "stim", "stim", "none", "prod", "prod"]) ...
                };
            panels2 = repmat({'event_boundary'}, size(panels1));
            panels2args = repmat({ebArgs}, size(panels1));
            
            bbArgs = {"stimText", 'Label', false, 'Color', [0 0 0 .3]};
            panels3 = repmat({'block_boundary'}, size(panels1));
            panels3args = repmat({bbArgs}, size(panels1));
            
            k = 1;
            for i = 1 : size(panels1,1)
                % first column for labels
                panels1{i,1} = ''; % not plotting spikes in the first column
                panels2{i,1} = ''; % not plotting event_boundary in the first column
                panels3args{i,1}{3} = true; % show labels in the first column
                
                % following columns for rasters
                for j = 2 : size(panels1,2)
                    if k > numel(uInd) || isnan(uInd(k))
                        panels1{i,j} = '';
                        panels2{i,j} = '';
                        panels3{i,j} = '';
                        continue
                    end
                    panels1args{i,j} = {uInd(k)};
                    k = k + 1;
                end
            end
            
            tWin = [-0.3 6];
            
            % Plotting
            f = gcf; clf
            NP.TaskBaseClass.PlotSE(se, panels1, 'PanelArgs', panels1args, 'TimeWindow', tWin);
            NP.TaskBaseClass.PlotSE(se, panels2, 'PanelArgs', panels2args, 'TimeWindow', tWin);
            NP.TaskBaseClass.PlotSE(se, panels3, 'PanelArgs', panels3args, 'TimeWindow', tWin);
        end
        
        % Figure 2, 3
        function LinkerRaster(ax, unitId, varargin)
            % Plot spike raster for a single unit
            % 
            %   LinkerRaster(unitId)
            %   LinkerRaster(unitId, ...)
            % 
            % See also LMV.Overview.Sentences
            
            % Process inputs
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('StimIdList', LMV.Param.stimIdList4, @isstring);
            p.parse(varargin{:});
            stimIdList = p.Results.StimIdList;
            
            % Load cache
            sUnit = NP.Unit.LoadUnitCache(unitId, varargin{:});
            
            % Put data into se
            se = LMV.SE.UnitCache2SE(sUnit, varargin{:});
            
            % Add unit attributes to clusTb
            clusTb = NP.Unit.GetClusTb(se);
            se.userData.ksMeta.clusTb = clusTb;
            
            % Make sentence table
            senTb = se.SplitConditions('stimId', 'taskValue');
            [~, ind] = MMath.SortLike(senTb.stimId, stimIdList);
            senTb = senTb(ind,:);
            
            nSen = height(senTb);
            for i = 1 : nSen
                se = senTb.se(i);
                [tt, tv] = se.GetTable('taskTime', 'taskValue');
                tWin = [tt.stimOn(1) tt.matchOff(1)] + [-.3 .3];
                
                % Raster
                NP.UnitPlot.RasterStack(ax, se, [], tWin, 1, 'HeightScale', 1.5, 'PETH', false, 'PeakSpikeRates', [], 'LineWidth', 1);
                
                % Labels
                tgePhase = {tt.stim(1), tt.prod{1}};
                y0 = 0;
                dy = -.4;
                for p = 1 : numel(tgePhase)
                    tge = Cut(tgePhase{p});
                    tge.Plot(y0+[-.2 0], 'FontSize', 0, 'Parent', ax, 'Style', 'patch');
                    for r = 1 : 3
                        tge(r:3:end).Plot(y0-.3+(4-r)*dy, 'FontSize', 12, 'Parent', ax);
                    end
                end
                
                ax.XGrid = 'off';
                ax.XLabel.String = [];
                ax.XTick = [];
                ax.YLim = [-1.8 2];
                % ax.YLabel.String = "Repeats";
            end
        end
        
        function SeqRespOverlay(ax, C, varargin)
            % Plot responses separated by different speech sequences
            % 
            %   SeqRespOverlay(ax, C)
            % 
            
            p = inputParser;
            p.addParameter('SeqStr', [], @(x) isstring(x) || iscellstr(x) || ischar(x) || isempty(x));
            p.addParameter('UnitIdx', 1, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('MaxSpikeRate', [], @(x) isscalar(x) && isnumeric(x));
            p.parse(varargin{:});
            seqStr = string(p.Results.SeqStr);
            uIdx = p.Results.UnitIdx;
            maxR = p.Results.MaxSpikeRate;
            
            cc = LMV.Param.GetTaskPhaseColors(["stim", "prod"]);
            ls = {':', '-'};
            
            for i = 1 : numel(C)
                % Organize table
                s = C{i};
                tb = s.seqTb;
                if ~isempty(seqStr)
                    m = ismember(tb.seqStr, seqStr);
                    tb = tb(m,:);
                end
                
                % PETH
                if isempty(tb)
                    continue
                end
                t = tb.tResp{1};
                hh = cat(3, tb.hh{:});
                hh = squeeze(hh(:,uIdx,:));
                ee = cat(3, tb.ee{:});
                ee = squeeze(ee(:,uIdx,:));
                if isempty(p.Results.MaxSpikeRate)
                    [~, c, k] = MMath.Normalize(hh, 'minmax');
                    hh = (hh-c) ./ k * 0.75;
                    ee = ee ./ k * 0.75;
                else
                    hh = hh ./ maxR;
                    ee = ee ./ maxR;
                end
                
                MPlot.PlotHistStack(t, hh, ee, 'Scaling', .8, 'Color', cc(i,:), 'LineStyle', ls{i}, 'LineWidth', 1.5);
            end
            
            % Labels
            tWin = [-.4 .4];
            for i = 1 : height(tb)
                tge = tb.seqTge{i};
                isIn = tge.GetTfield('tmax') > tWin(1) & tge < tWin(2);
                tge = tge(isIn);
                tge.PlotTiers(i-1, 0.25, 0, 'TierInd', [1 3], 'Color', {'k', 'k'}, 'FontSize', 12);
                
                hold(ax, 'on');
            end
            
            ax.XLim = tWin;
            ax.YLim = [0 1+0.4];
            ax.XTickLabel = [];
            ax.YTick = [];
            
            % ax.XLabel.String = sprintf("Time from /%s/ (s)", MLing.ARPA2IPA(s.seed));
            MPlot.Axes(ax);
        end
        
        function SentenceRespOverlay(ax, sPETH, posTb, varargin)
            % Plot overlay of mean spike rates in different sentences
            % 
            %   SentenceRespOverlay(ax, sPETH, posTb)
            %   SentenceRespOverlay(ax, sPETH, posTb, ..., "MinScore", 0)
            %   SentenceRespOverlay(ax, sPETH, posTb, ..., "SentenceMask", true(size(sPETH.stimId)))
            %   SentenceRespOverlay(ax, sPETH, posTb, ..., "ShowNonLinked", true)
            %   SentenceRespOverlay(ax, sPETH, posTb, ..., "ShowText", true)
            %   SentenceRespOverlay(ax, sPETH, posTb, ..., "ShowScore", true)
            % 
            
            p = inputParser;
            p.addParameter("MinScore", LMV.Linker.scoreTh, @isnumeric);
            p.addParameter("SentenceMask", true(size(sPETH.stimId)), @islogical);
            p.addParameter("ShowNonLinked", true, @islogical);
            p.addParameter("ShowText", true, @islogical);
            p.addParameter("ShowScore", true, @islogical);
            p.parse(varargin{:});
            minScore = p.Results.MinScore;
            isSen = p.Results.SentenceMask(:)';
            isNL = p.Results.ShowNonLinked;
            isText = p.Results.ShowText;
            isScore = p.Results.ShowScore;
            
            if istable(sPETH)
                sPETH = table2struct(sPETH);
            end
            
            % Find the number of links in each sentence
            if exist("posTb", "var") && ~isempty(posTb)
                nLinks = arrayfun(@(x) sum(posTb.stimId==x & posTb.score>minScore), sPETH.stimId);
            else
                nLinks = zeros(size(sPETH.stimId));
            end
            nLinks = nLinks(:)';
            nLinks(~isSen) = 0;
            
            % Get sentence colors
            cc = LMV.Param.GetSentenceColors(sPETH.stimId) * 0.9; % make them slightly darker
            
            % Plot sentence responses
            t = sPETH.time;
            r0 = sPETH.resp;
            rSig = sPETH.sigResp;
            if isNL
                plot(t, r0, 'Color', [0 0 0 .15]); hold on
                plot(t, rSig, 'Color', [0 0 0 .15], 'LineWidth', 1.5);
            end
            
            rLk = r0;
            rLk(:,~nLinks) = NaN;
            hh = plot(t, rLk, 'LineWidth', 1.5); hold on
            for i = 1 : numel(hh)
                hh(i).Color = [cc(i,:) 0.3];
            end
            
            rSig(:,~nLinks) = NaN;
            hh = plot(t, rSig, 'LineWidth', 1.5);
            for i = 1 : numel(hh)
                hh(i).Color = cc(i,:);
            end
            
            yRange = [0 max(r0(:))*1.5];
            ax.YLim = yRange;
            
            % Plot phase shading
            se = LMV.SE.UnitCache2SE({sPETH.uCache});
            colNames = {'stim', 'prod'};
            phaseRange = yRange(2)*[0.95 1];
            NP.TaskBaseClass.PlotEventWindows(ax, se, 1, t([1 end])', colNames, ...
                'YRange', phaseRange, 'Colors', LMV.Param.GetTaskPhaseColors(colNames), 'Alpha', 0.1, 'Text', false);
            
            % Label linked sentences
            senInd = find(nLinks);
            yText = linspace(yRange(2)*0.95, yRange(2)*0.7, sum(nLinks)+2);
            N = 2;
            for i = senInd(:)'
                % Label the linked sentence
                if isText
                    stimText = sPETH.stim(i).GetParentLabel;
                    text(-0.1, yText(N), stimText, 'Color', cc(i,:), 'FontSize', 6, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
                end
                
                % Mark links
                lkInd = find(posTb.stimId==sPETH.stimId(i) & posTb.score > minScore);
                for j = 1 : numel(lkInd)
                    k = lkInd(j);
                    plot(posTb.tM2(k,:), yText([N N]), '*-', 'Color', cc(i,:));
                    if isScore
                        text(posTb.tM2(k,2), yText(N), sprintf("%.2g", posTb.score(k)), 'FontSize', 6);
                    end
                    N = N + 1;
                end
            end
            
            ax.XLim = t([1 end]);
            ax.YLabel.String = "Spk/s";
            ax.Title.String = "u"+sPETH.clusId;
            MPlot.Axes(ax);
        end
        
        % Figure 6
        function PlotPeriEventSeqArray(C, varargin)
            % Plot an array of panels each contains rasters and PETHs separated by different speech sequences
            % 
            %   PlotPeriEventSeqArray(C)
            % 
            
            p = inputParser;
            p.addParameter('UnitIdx', 1, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('MaxSpikeRate', [], @(x) isscalar(x) && isnumeric(x));
            p.addParameter('FeatureName', [], @(x) ischar(x) || isstring(x) || iscellstr(x) || isempty(x));
            p.parse(varargin{:});
            uIdx = p.Results.UnitIdx;
            maxR = p.Results.MaxSpikeRate;
            feat = p.Results.FeatureName;
            
            if isempty(maxR) && isempty(feat)
                % Find the maximum PETH peak spike rate
                maxR = zeros(size(C));
                for i = 1 : numel(C)
                    if isempty(C{i}.seqTb)
                        continue
                    end
                    maxR(i) = max(cellfun(@(x) x.pkVal(uIdx), C{i}.seqTb.stats));
                end
                maxR = max(maxR(:));
            end
            
            % Plot
            [nRow, nCol] = size(C);
            tl = tiledlayout(nRow, nCol);
            tl.Padding = "compact";
            for i = 1 : nRow
                for j = 1 : nCol
                    ax = nexttile;
                    LMV.Fig.PlotSeqResp(ax, C{i,j}, 'MaxSpikeRate', maxR, 'UnitIdx', uIdx);
                end
            end
        end
        
        function PlotSeqResp(ax, s, varargin)
            % Plot rasters and PETHs separated by different speech sequences
            % 
            %   PlotSeqResp(ax, s)
            % 
            
            p = inputParser;
            p.addParameter('UnitIdx', 1, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('MaxSpikeRate', [], @(x) isscalar(x) && isnumeric(x));
            p.addParameter('SortByRate', true, @islogical);
            p.parse(varargin{:});
            uIdx = p.Results.UnitIdx;
            maxR = p.Results.MaxSpikeRate;
            isSort = p.Results.SortByRate;
            
            % Organize table
            tb = s.seqTb;
            if isempty(tb)
                return
            end
            tb(tb.nSample<3,:) = []; % remove seqs that have too few samples
            if isempty(tb)
                return
            end
            r = round(cellfun(@(x) x.optiResp(uIdx), tb.stats));
            if isSort
                [r, I] = sort(r, 'descend');
                tb = tb(I,:);
            end
            
            % PETH
            hh = cat(3, tb.hh{:});
            hh = squeeze(hh(:,uIdx,:));
            ee = cat(3, tb.ee{:});
            ee = squeeze(ee(:,uIdx,:));
            if isempty(maxR)
                maxR = max(hh(:));
            end
            hh = hh ./ maxR;
            ee = ee ./ maxR;
            MPlot.PlotHistStack(tb.tResp{1}, hh, ee, 'Scaling', 0.8);
            
            % Raster
            st = cellfun(@(x) x(:,uIdx), tb.spikeTime, 'Uni', false);
            MPlot.PlotRasterStack(st, 1:height(tb), 'HeightScale', 0.8, 'Parent', ax);
            
            % Labels
            [~, I] = sort(s.positions);
            seqLb = cellfun(@(x) join(MLing.ARPA2IPA(x(I).GetParentLabel), " "), tb.seqTge);
            labelArray = [seqLb'; r'+"Hz"];
            labelArray(ismissing(labelArray)) = "";
            tickLabels = sprintf('%s\\newline%s\n', labelArray(:));
            ax.YTickLabel = strtrim(tickLabels);
            
            if s.nElem == 1
                ax.Title.String = sprintf("/%s/ ", MLing.ARPA2IPA(seqLb));
            else
                ax.Title.String = sprintf("/%s/ %+i %s", MLing.ARPA2IPA(s.seed), s.nElem-1, s.level);
            end
            MPlot.Axes(ax);
        end
        
        function PlotArticResp(mdl, seqTb, featName, itvl, axArray)
            % 
            % 
            %   PlotArticResp(mdl, seqTb, featName, itvl, axArray)
            % 
            
            A = cat(1, seqTb.(featName){:});
            A = cat(2, A{:})';
            
            tFeat = seqTb.tFeat{1}{1};
            [~, I] = min(abs(tFeat - 0));
            [a, I] = sort(A(:,I), 1, "ascend");
            A = A(I,:);
            
            tResp = seqTb.tResp{1};
            st = cat(1, seqTb.spikeTime{:});
            st = cellfun(@(x) x - mdl.r2t, st, 'Uni', false);
            st = st(I);
            
            m = ~isoutlier(a, "median", "ThresholdFactor", 2);
            a = a(m);
            A = A(m,:);
            st = st(m);
            
            if ~exist('itvl', 'var')
                itvl = 1;
            end
            iSub = 1:itvl:numel(a);
            N = numel(iSub);
            
            % Plotting
            if exist('axArray', 'var')
                ax = axArray(1);
            else
                ax = nexttile;
            end
            MPlot.PlotRaster(st(iSub), [], N/100, 'ColorArray', [0 0 0], 'Parent', ax);
            % ax.XLim = tResp([1 end]);
            ax.XLim = [-1 1]*0.1; % + mdl.r2t;
            ax.YLim = [0 N+1];
            ax.YTick = [];
            ax.XLabel.String = "Time + \Deltat (s)";
            ax.YLabel.String = sprintf("Samples (n = %i)", N);
            ax.Title.String = mdl.resps;
            MPlot.Axes(ax);
            
            if exist('axArray', 'var')
                ax = axArray(2);
            else
                ax = nexttile;
            end
            plot(ax, a(iSub), 1:N);
            ax.YLim = [0 N+1];
            ax.YTick = [];
            ax.XLabel.String = "Feat. value (AU)";
            ax.Title.String = featName;
            ax.Title.Interpreter = "none";
            MPlot.Axes(ax);
            
            if exist('axArray', 'var')
                ax = axArray(3);
            else
                ax = nexttile;
            end
            imagesc(ax, tFeat, [], flip(A(iSub,:),1));
            ax.XLim = [-1 1]*0.1;
            ax.YTick = [];
            ax.XLabel.String = "Time (s)";
            ax.Title.String = featName;
            ax.Title.Interpreter = "none";
            MPlot.Axes(ax);
        end
        
        function PlotSpectralResp(mdl, seqTb, binIdx, binFreq, itvl, axArray)
            % 
            % 
            %   PlotSpectralResp(mdl, seqTb, binIdx, binFreq, itvl, axArray)
            % 
            
            A = cat(1, seqTb.mel{:});
            A = cat(3, A{:});
            A = squeeze(A(:,binIdx,:))';
            
            tFeat = seqTb.tFeat{1}{1};
            [~, I] = min(abs(tFeat - 0));
            [a, I] = sort(A(:,I), 1, "ascend");
            A = A(I,:);
            
            tResp = seqTb.tResp{1};
            st = cat(1, seqTb.spikeTime{:});
            st = cellfun(@(x) x - mdl.r2t, st, 'Uni', false);
            st = st(I);
            
            m = ~isoutlier(a, "median", "ThresholdFactor", 2);
            a = a(m);
            A = A(m,:);
            st = st(m);
            
            if ~exist('itvl', 'var')
                itvl = 1;
            end
            iSub = 1:itvl:numel(a);
            N = numel(iSub);
            
            % Plotting
            if exist('axArray', 'var')
                ax = axArray(1);
            else
                ax = nexttile;
            end
            MPlot.PlotRaster(st(iSub), [], N/100, 'ColorArray', [0 0 0], 'Parent', ax);
            % ax.XLim = tResp([1 end]);
            ax.XLim = [-1 1]*0.1; % + mdl.r2t;
            ax.YLim = [0 N+1];
            ax.YTick = [];
            ax.XLabel.String = "Time + \Deltat (s)";
            ax.YLabel.String = sprintf("Samples (n = %i)", N);
            ax.Title.String = mdl.resps;
            MPlot.Axes(ax);
            
            if exist('axArray', 'var')
                ax = axArray(2);
            else
                ax = nexttile;
            end
            plot(ax, a(iSub), 1:N);
            ax.YLim = [0 N+1];
            ax.YTick = [];
            ax.XLabel.String = "Feat. value (AU)";
            ax.Title.String = round(binFreq)+" Hz";
            ax.Title.Interpreter = "none";
            MPlot.Axes(ax);
            
            if exist('axArray', 'var')
                ax = axArray(3);
            else
                ax = nexttile;
            end
            imagesc(ax, tFeat, [], flip(A(iSub,:),1));
            ax.XLim = [-1 1]*0.1;
            ax.YTick = [];
            ax.XLabel.String = "Time (s)";
            ax.Title.String = round(binFreq)+" Hz";
            ax.Title.Interpreter = "none";
            MPlot.Axes(ax);
        end
        
    end
    
end