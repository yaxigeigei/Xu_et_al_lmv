classdef TS
    % Analysis of response time shifts
    
    methods(Static)
        function clusTb = ComputeRespTimeShift(clusTb, posTb, seArray, targetPhase)
            % Compute time shifts of response timing between early and late repeats
            % 
            %   clusTb = ComputeRespTimeShift(clusTb, posTb, seArray, targetPhase)
            %
            % Inputs
            %   clusTb          A table of unit information
            %   posTb           A table of linked positions
            %   seArray         An array of SE objects
            %   targetPhase     A string of target phase name ("stim" or "prod")
            %
            % Outputs
            %   clusTb          A table of unit information with time shift information
            % 

            recIds = arrayfun(@(x) string(NP.SE.GetID(x)), seArray);
            
            tWin = LMV.Linker.respTimeWin;
            tRs = (tWin(1) : 0.01 : tWin(2))';
                
            % Loop through units
            for u = 1 : height(clusTb)
                % Find linked positions of this unit
                cid = clusTb.clusId(u);
                pTb = posTb(posTb.clusId==cid, :);
                
                % Get recording data
                isRec = recIds == clusTb.recId(u);
                se = seArray(isRec);
                [tt, tv, sr] = se.GetTable("taskTime", "taskValue", "spikeRate");
                
                % Loop through linked positions of this unit
                for p = 1 : height(pTb)
                    % Find target trials
                    stimId = pTb.stimId(p);
                    epInd = find(tv.stimId==stimId & tv.phase==targetPhase);
                    
                    % Find time windows
                    if targetPhase == "stim"
                        tCenter = pTb.t0(p) + double(tt.stim(epInd)); % re-reference to cue onset
                    else
                        tCenter = pTb.t0(p) + double(tt.stim(epInd-1)); % re-reference to cue onset
                    end
                    sliceWins = NaN(height(tv), 2);
                    sliceWins(epInd,:) = tCenter + tWin;
                    
                    % Slice spike rates
                    srSlices = se.SliceTimeSeries(sr(:, ["time", "u"+cid]), sliceWins);
                    srSlices = srSlices(epInd, :);
                    srSlices.time = cellfun(@(t,t0) t-t0, srSlices.time, num2cell(tCenter), 'Uni', false); % re-reference to linking position
                    
                    % Get responses
                    rTb = table;
                    rTb.repIdx = (1:numel(epInd))';
                    rTb.trialNum = tv.trialNum(epInd);
                    rTb.rawResp = NaN(height(rTb), numel(tRs));
                    rTb.smoothResp = NaN(height(rTb), numel(tRs));
                    rTb.pkTime(:) = NaN;
                    rTb.pkRate(:) = 0;
                    for i = 1 : height(srSlices)
                        t = srSlices.time{i};
                        r = srSlices.(2){i};
                        if numel(t) < 40 % ignore trials with less than 0.4s of data
                            continue
                        end
                        r = interp1(t, r, tRs, "linear", 0);
                        rTb.rawResp(i,:) = r;
                        r = MNeuro.Filter1(r, 100, 'gaussian', 0.05);
                        rTb.smoothResp(i,:) = r;
                        [rPk, iPk] = max(r);
                        rTb.pkTime(i) = tRs(iPk);
                        rTb.pkRate(i) = rPk;
                    end
                    isSilent = rTb.pkRate == 0;  % trials with no response
                    isEdge = rTb.pkTime == tRs(1) | rTb.pkTime == tRs(end); % max value found at the edges
                    rTb(isSilent | isEdge, :) = [];
                    nResp = height(rTb);
                    
                    pTb.rTime{p} = tRs;
                    pTb.rTb{p} = rTb;
                    pTb.numRepeat(p) = numel(epInd);
                    pTb.numResp(p) = nResp;
                    if nResp < 4
                        continue
                    end
                    
                    % Test difference in responses between early and late repeats
                    iHalf = ceil(nResp/2);
                    x1 = rTb.pkTime(1:iHalf);
                    x2 = rTb.pkTime(iHalf+1:end);
                    [pv, ~, stats] = ranksum(x1, x2);
                    
                    nx1 = numel(x1);
                    nx2 = numel(x2);
                    W  = stats.ranksum;         % Wilcoxon W from ranksum
                    U  = W - nx1*(nx1+1)/2;     % Mann-Whitney U
                    cliff = 1 - 2*U/(nx1*nx2);  % rank-biserial correlation (== Cliff's δ)
                    
                    pTb.diff(p) = mean(x2) - mean(x1);
                    pTb.diffPval(p) = pv;
                    pTb.diffSize(p) = cliff;
                    
                    % Test correlation between interval length and change in response time
                    dx = zscore(diff(rTb.trialNum));
                    dy = zscore(diff(rTb.pkTime));
                    [pTb.rho(p), pTb.rhoPval(p)] = corr(dx, dy);
                end
                
                clusTb.posTb{u} = pTb;
            end
        end
        
        function axs = PlotTimeShiftCombo(posTb, varargin)
            % Plot scatter plot of time shifts vs effect sizes
            % 
            %   PlotTimeShiftCombo(posTb)
            %   PlotTimeShiftCombo(posTb, ..., "Axes", [])
            % 

            p = inputParser;
            p.addParameter("XLim", [-100 100], @isnumeric);
            p.parse(varargin{:});
            xlims = p.Results.XLim;

            % Set up tiles
            rowDist = [1 3];
            colDist = [3 1];
            tl = tiledlayout(sum(rowDist), sum(colDist));
            tl.Padding = "compact";

            % Plot time shift histogram on top
            ntArgs = MPlot.FindTileInd(rowDist, colDist, 1, 1);
            ax1 = nexttile(ntArgs{:});
            LMV.Linker.TS.PlotTimeShiftHistogram(posTb, "Axes", ax1);
            ax1.YLim = [0 .4];
            ax1.XTickLabel = [];
            ax1.XLabel.String = [];
            ax1.Title.String = [];

            % Plot time shift vs effect size
            ntArgs = MPlot.FindTileInd(rowDist, colDist, 2, 1);
            ax2 = nexttile(ntArgs{:});
            LMV.Linker.TS.PlotTimeShiftEffectScatter(posTb, "Axes", ax2, "XLim", xlims);
            ax2.Title.String = [];
            
            % Plot time shift effect histogram on right
            ntArgs = MPlot.FindTileInd(rowDist, colDist, 2, 2);
            ax3 = nexttile(ntArgs{:});
            LMV.Linker.TS.PlotTimeShiftEffectHistogram(posTb, "Axes", ax3, "Orientation", "horizontal");
            ax3.XLim = [0 .4];
            ax3.YTickLabel = [];
            ax3.YLabel.String = [];
            ax3.Title.String = [];

            axs = [ax1, ax2, ax3];
        end

        function PlotTimeShiftEffectScatter(posTb, varargin)
            % Plot scatter plot of time shifts vs effect sizes
            % 
            %   PlotTimeShiftEffectScatter(posTb)
            %   PlotTimeShiftEffectScatter(posTb, ..., "Axes", [])
            %   PlotTimeShiftEffectScatter(posTb, ..., "XLim", [-100 100])
            % 
            
            p = inputParser;
            p.addParameter("Axes", [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.addParameter("XLim", [-100 100], @isnumeric);
            p.addParameter("LineArgs", {}, @iscell);
            p.parse(varargin{:});
            ax = p.Results.Axes;
            xlims = p.Results.XLim;
            lineArgs = p.Results.LineArgs;

            x = posTb.diff * 1e3;
            y = posTb.diffSize;
            
            if isempty(ax)
                ax = nexttile;
            end
            
            h = plot(ax, x, y, 'o', lineArgs{:}); hold on
            plot(ax, xlims, [0 0], '-', 'Color', [0 0 0 .2]);
            ax.XLim = xlims;
            ax.YLim = [-1 1];
            ax.XLabel.String = "Time shift (ms)";
            ax.YLabel.String = "Effect size (r)";
            ax.Title.String = sprintf("n = %d", numel(x));
            MPlot.Axes(ax);
            
            % Set up brushing callback
            h.UserData.forBrush = true;
            ax.UserData.posTb = posTb;
            brushObj = brush(gcf);
            brushObj.ActionPostCallback = @LMV.Linker.TS.LadderOnBrush;
        end
        
        function PlotTimeShiftHistogram(posTb, varargin)
            % Plot histogram of time shifts
            % 
            %   PlotTimeShiftHistogram(posTb)
            %   PlotTimeShiftHistogram(posTb, ..., "Axes", [])
            %   PlotTimeShiftHistogram(posTb, ..., "Edges", -100:20:100)
            % 
            
            p = inputParser;
            p.addParameter("Axes", [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.addParameter("Edges", -100:20:100, @isnumeric);
            p.addParameter("HistArgs", {}, @iscell);
            p.parse(varargin{:});
            ax = p.Results.Axes;
            tEdges = p.Results.Edges;
            xlims = tEdges([1 end]);
            histArgs = p.Results.HistArgs;
            x = posTb.diff * 1e3;
            m = x > xlims(1) & x < xlims(2);
            x = x(m);
            [p, ~, stats] = signrank(x);
            
            if isempty(ax)
                ax = nexttile;
            end
            histogram(ax, x, tEdges, "Normalization", "probability", "EdgeColor", "none", histArgs{:});
            ax.Title.String = sprintf("p = %.2g", p);
            ax.XLim = xlims;
            ax.XLabel.String = "Time shift (ms)";
            ax.YLabel.String = "Frac.";
            MPlot.Axes(ax);
        end
        
        function PlotTimeShiftEffectHistogram(posTb, varargin)
            % Plot histogram of time shifts
            % 
            %   PlotTimeShiftEffectHistogram(posTb)
            %   PlotTimeShiftEffectHistogram(posTb, ..., "Axes", [])
            %   PlotTimeShiftEffectHistogram(posTb, ..., "Edges", -1:0.2:1)
            %   PlotTimeShiftEffectHistogram(posTb, ..., "Orientation", "vertical")
            %   PlotTimeShiftEffectHistogram(posTb, ..., "HistArgs", {})
            
            p = inputParser;
            p.addParameter("Axes", [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.addParameter("Edges", -1.1:0.2:1.1, @isnumeric);
            p.addParameter("HistArgs", {}, @iscell);
            p.addParameter("Orientation", "vertical", @(x) any(strcmp(x, ["vertical", "horizontal"])));
            p.parse(varargin{:});
            ax = p.Results.Axes;
            edges = p.Results.Edges;
            orient = p.Results.Orientation;
            histArgs = p.Results.HistArgs;

            y = posTb.diffSize;
            [p, ~, stats] = signrank(y);
            
            if isempty(ax)
                ax = nexttile;
            end
            histogram(ax, y, edges, "Normalization", "probability", "Orientation", orient, 'EdgeColor', 'none', histArgs{:});
            centers = edges(1:end-1) + diff(edges)/2;
            xlims = centers([1 end]);
            if orient == "vertical"
                ax.XLim = xlims;
                ax.XLabel.String = "Effect size (r)";
                ax.YLabel.String = "Frac.";
            else
                ax.YLim = xlims;
                ax.YLabel.String = "Effect size (r)";
                ax.XLabel.String = "Frac.";
            end
            ax.Title.String = sprintf("p = %.2g", p);
            MPlot.Axes(ax);
        end
        
        function LadderOnBrush(fig, axesStruct)
            % 
            %   LadderOnBrush(fig, axesStruct)
            % 
            
            % Iterate through plotted handles
            hh = axesStruct.Axes.Children;
            b = [];
            for i = 1 : numel(hh)
                % Check if handle is labeled for brushing
                s = hh(i).UserData;
                if ~isfield(s, 'forBrush') || ~s.forBrush
                    continue
                end
                
                % Get the brushing mask
                b = logical(hh(i).BrushData);
            end
            
            % Check if brushed any data
            if ~any(b)
                return
            end
            
            % Sample a subset if too many
            b = find(b);
            if numel(b) > 8
                b = randsample(b, 8);
            end
            
            % Create figure if absent
            if isempty(fig.UserData) || ~isvalid(fig.UserData.cbFig)
                fig.UserData.cbFig = MPlot.Figure( ...
                    'Name', 'Selected positions', ...
                    'NumberTitle', 'off', ...
                    'Menubar', 'none', ...
                    'Toolbar', 'figure', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
            end
            
            % Highlight and clear figure
            figure(fig.UserData.cbFig);
            clf(fig.UserData.cbFig);
            tl = tiledlayout(fig.UserData.cbFig, 2, 4);
            tl.Padding = "compact";
            
            % Make plots
            pTb = axesStruct.Axes.UserData.posTb(b,:);
            LMV.Linker.TS.PlotRespLadder(pTb);
        end
        
        function PlotRespLadder(posTb)
            % Plot response ladder for selected positions
            % 
            %   PlotRespLadder(posTb)
            % 
            
            for i = 1 : height(posTb)
                cid = posTb.clusId(i);
                tge = posTb.seqTge{i};
                t = posTb.rTime{i};
                rTb = posTb.rTb{i};
                rRaw = rTb.rawResp';
                rSmooth = rTb.smoothResp';
                yBase = 1 : size(rRaw, 2);
                tPk = rTb.pkTime;

                % Normalize response amplitudes
                k = 1./max(rSmooth(:)) * 0.8;
                rRaw = rRaw * k;
                rSmooth = rSmooth * k;
                
                % Assign colors for the first and second half of the repeats
                nRep = numel(tPk);
                nHalf = ceil(nRep/2);
                cMap = lines(2);
                ccArray = repmat(cMap(1,:), nRep, 1);
                ccArray(nHalf+1:end,:) = repmat(cMap(2,:), nRep-nHalf, 1);
                
                ax = nexttile;
                % MPlot.PlotTraceLadder(t, -rRaw, yBase, 'Color', [0 0 0 .15], 'Parent', ax); hold(ax, "on")
                MPlot.PlotTraceLadder(t, -rSmooth, yBase, 'ColorArray', ccArray, 'Parent', ax); hold(ax, "on")
                MPlot.PlotPointAsLine(tPk, yBase-0.5, 1, 'Color', [0 0 0 .3], 'Parent', ax);
                tge.PlotTiers(-1.25, 0.5);
                text(ax, 0, -1.5, sprintf("r = %.2f", posTb.diffSize(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
                ax.XLim = t([1 end]);
                ax.XTick = t(1) : 0.2 : t(end);
                ax.XTickLabel([2 4]) = {'', ''};
                ax.XTickLabelRotation = 0;
                ax.YLim = [-2 8];
                ax.YTick = 1:8;
                ax.XLabel.String = "Time from link (s)";
                ax.YLabel.String = "Repeats";
                ax.Title.String = sprintf("u%d", cid);
                MPlot.Axes(ax);
            end
        end
        
        function PlotIntervalRelations(posTb, varargin)
            % Plot scatter plot of interval vs response time changes
            % 
            %   PlotIntervalRelations(posTb)
            %   PlotIntervalRelations(posTb, ..., "Axes", [])
            % 
            
            p = inputParser;
            p.addParameter("Axes", [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
            ax = p.Results.Axes;
            
            if isempty(ax)
                ax = nexttile;
            end
            
            for i = 1 : height(posTb)
                x = diff(posTb.rTb{i}.trialNum);
                y = diff(posTb.rTb{i}.pkTime) * 1e3;
                [x, I] = sort(x);
                y = y(I);
                plot(ax, x, y, '.-', 'Color', [0 0 0 .15]); hold on
            end
            ax.XLim(1) = 0;
            ax.XLabel.String = "Interval (# of trials)";
            ax.YLabel.String = "\DeltaResponse time (ms)";
            MPlot.Axes(ax);
        end
        
        function PlotIntervalCorrHistogram(posTb, varargin)
            % Plot histogram of interval correlation coefficients
            % 
            %   PlotIntervalCorrHistogram(posTb)
            %   PlotIntervalCorrHistogram(posTb, ..., "Axes", [])
            % 
            
            p = inputParser;
            p.addParameter("Axes", [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
            ax = p.Results.Axes;
            
            if isempty(ax)
                ax = nexttile;
            end
            
            r = posTb.rho;
            [p, ~, stats] = signrank(r);
            histogram(ax, r, 15);
            ax.XLim = [-1 1];
            text(ax, ax.XLim(1)*0.9, ax.YLim(2), sprintf("p = %.2g", p), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
            ax.XLabel.String = "Pearson's r";
            ax.YLabel.String = "# of positions";
            MPlot.Axes(ax);
        end

    end
end