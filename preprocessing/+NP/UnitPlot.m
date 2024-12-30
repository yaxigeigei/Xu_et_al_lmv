classdef UnitPlot < NP.PlotBaseClass
    
    properties(Constant)
        suColor = [.2 .7 .2];
        suAlpha = 0.5;
    end
    
    methods(Static)
        % Quality control
        function WaveformArray(clusTb, figPath)
            % Plot an array of mean spike waveform
            
            nRows = 9;
            nCols = 15;
            npp = nRows * nCols; % number of plots per page
            p = 0; % number of page plotted
            
            while p*npp < height(clusTb)
                a = p*npp+1;
                b = min(p*npp+npp, height(clusTb));
                subTb = clusTb(a:b,:);
                
                f = MPlot.Figure(214040); clf
                f.WindowState = 'maximized';
                pause(0.1);
                
                panels = cell(nCols, nRows);
                panelArgs = panels;
                for i = 1 : height(subTb)
                    panels{i} = 'waveform';
                    panelArgs{i} = {i};
                end
                NP.UnitPlot.Array(subTb, panels', 'PanelArgs', panelArgs');
                
                if nargin > 1
                    exportgraphics(f, figPath + "_page" + (p+1) + ".png");
                end
                
                p = p + 1;
            end
        end
        
        function Array(clusTb, panels, varargin)
            % Plot various speech features, unit raster, and PETHs aligned in time
            % 
            %   QC(clusTb, panels)
            %   QC(clusTb, panels, rowDist)
            %   QC(clusTb, panels, rowDist, colDist)
            %   QC(..., 'PanelArgs', {{}})
            % 
            % Inputs
            %   clusTb          The common data to plot from in each panel.
            %   panels          A m-by-n cell array. The shape specifies that the figure has m rows and n 
            %                   columns of plots (or Axes). The variable in each cell specifies what to 
            %                   plot, and can be one of the following:
            %                   1) A char string of plot name. 'cla' clears the Axes.
            %                   2) A positional index of units in the spikeTime table.
            %                   3) Empty, then nothing will be done to the Axes.
            %   rowDist         A m-element vector specifying integer height (i.e. row) ratios of the panels.
            %   colDist         A n-element vector specifying integer width (i.e. column) ratios of the panels.
            %   'PanelArgs'     A cell array of cell arrays. The size of the outer array should match that of 
            %                   panels. Each inside array contains the additional arguments for plotting.
            % 
            
            % Parse user inputs
            p = inputParser();
            p.addOptional('rowDist', [], @isnumeric);
            p.addOptional('colDist', [], @isnumeric);
            p.addParameter('PanelArgs', {{}}, @iscell);
            p.parse(varargin{:});
            rowDist = p.Results.rowDist;
            colDist = p.Results.colDist;
            panelArgs = p.Results.PanelArgs;
            
            % Resolve the figure layout
            if isempty(rowDist)
                rowDist = ones(size(panels,1), 1);
            end
            if isempty(colDist)
                colDist = ones(size(panels,2), 1);
            end
            nRow = sum(rowDist);
            nCol = sum(colDist);
            tl = tiledlayout(nRow, nCol);
            tl.Padding = 'tight';
            
            % Match argument size to panels
            if isscalar(panelArgs)
                panelArgs = repmat(panelArgs, size(panels));
            end
            
            % Plot pannels
            for i = 1 : size(panels, 1)
                for j = 1 : size(panels, 2)
                    % Create Axes
                    ntArgs = MPlot.FindTileInd(rowDist, colDist, i, j);
                    ax = nexttile(ntArgs{:});
                    hold(ax, 'on');
                    
                    % Plot a panel
                    p = panels{i,j};
                    args = panelArgs{i,j};
                    if isempty(p)
                        axis(ax, 'off'); % hide axes
                    elseif ischar(p)
                        switch p
                            case 'waveform'
                                NP.UnitPlot.Waveform(ax, clusTb, args{:});
                            case 'isi'
                                NP.UnitPlot.ISI(ax, clusTb, args{:});
                            case 'logISI'
                                NP.UnitPlot.LogISI(ax, clusTb, args{:});
                            case 'cla'
                                cla(ax);
                        end
                    elseif isnumeric(p)
                        % plot something
                    end
                    
                    % Omit x-label if it's not the last row
                    if i ~= size(panels, 1)
                        ax.XLabel.String = [];
                    end
                end
            end
        end
        
        function Waveform(ax, clusTb, uIdx)
            % Plot mean waveform as trace and heatmap overlay of a given unit
            % 
            %   Waveform(ax, clusTb, uIdx)
            % 
            
            if nargin < 3
                uIdx = 1;
            end
            if istable(clusTb)
                s = table2struct(clusTb(uIdx,:));
            else
                s = clusTb(uIdx);
            end
            
            % Scale waveform such that the mean SD equals 0.5
            m = s.waveformMed;
            e = s.waveformSD;
            u = mean(e(:)) / 0.5;
            m = m ./ u;
            e = e ./ u;
            mMax = max(abs(m(:)));
            
            % Make coordinates
            [nChan, nTime] = size(m);
            c0 = floor(nChan/2);
            c = (1 : nChan) - c0;
            t = 1 : nTime;
            
            % Scale small waveform up
            k = size(m,1)/mMax/4;
            k = 1;
            m0 = -m(c0,:)*k;
            e0 = e(c0,:)*k;
            m0Max = max(abs(m0));
            k = min(3, 1/m0Max*3);
            m0 = m0 * k;
            e0 = e0 * k;
            
            % Plot heatmap
            imagesc(ax, t, c, m); hold on
            
            % Plot trace
            MPlot.ErrorShade(t, m0, e0, 'Parent', ax);
            plot(ax, t, m0, 'Color', [0 0 0]);
            
            axis(ax, 'tight');
            axis(ax, 'off');
            ax.YLim = c([1 end]);
            colormap(ax, MPlot.PolarMap());
            ax.CLim = [-1 1] * mMax;
            ax.Title.String = sprintf('u%i@%ium', s.clusId, s.depth);
        end
        
        function WaveformAtDepth(ax, cTb, x, varargin)
            % Plot waveforms at unit depth
            % 
            %   WaveformAtDepth(ax, clusTb, x)
            % 
            
            p = inputParser;
            p.addParameter('TimeScaling', 0.05, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('WaveformScaling', 0.2, @(x) isnumeric(x) && isscalar(x));
            p.parse(varargin{:});
            tScale = p.Results.TimeScaling;
            wScale = p.Results.WaveformScaling;
            
            % Get waveforms at middle channel
            W = cat(3, cTb.waveformMed{:});
            [nCh, nTm, nW] = size(W);
            W = squeeze(W(round(nCh/2),:,:));
            W = W./max(abs(W),[],1);
            W = -W*wScale;
            t = ((0:nTm-1)'-nTm/2) / 30; % convert to millisecond
            t = t * tScale;
            
            % Get depths
            d = cTb.depth / 1e3;
            
            % Find center coordinates
            if isscalar(x)
                x = MPlot.ViolinScatter(x, d, 0:.5:8, 'IsPlot', false);
            end
            
            % Get colors
            % cc = lines(max(height(cTb), 4));
            cc = brewermap(max(height(cTb), 4), 'Dark2');
            cc = cc([3 1 2 4],:);
            if ismember("wfId", cTb.Properties.VariableNames)
                cc = cc(double(cTb.wfId),:);
            end
            
            % Plot waveforms
            for i = 1 : nW
                plot(t+x(i), W(:,i)+d(i), 'Color', cc(i,:), 'Parent', ax);
                hold(ax, 'on');
            end
            
            ax.YDir = "reverse";
            MPlot.Axes(gca);
        end
        
        function ISIArray(clusTb, figPath)
            % Plot an array of mean spike waveform
            
            nRows = 9;
            nCols = 15;
            npp = nRows * nCols; % number of plots per page
            p = 0; % number of page plotted
            
            while p*npp < height(clusTb)
                a = p*npp+1;
                b = min(p*npp+npp, height(clusTb));
                subTb = clusTb(a:b,:);
                
                f = MPlot.Figure(214041); clf
                f.WindowState = 'maximized';
                pause(0.1);
                tl = tiledlayout('flow');
                tl.Padding = 'compact';
                
                panels = cell(nCols, nRows);
                panelArgs = panels;
                for i = 1 : height(subTb)
                    panels{i} = 'isi';
                    panelArgs{i} = {i};
                end
                NP.UnitPlot.Array(subTb, panels', 'PanelArgs', panelArgs');
                
                if nargin > 1
                    exportgraphics(f, figPath + "_page" + (p+1) + ".png");
                end
                
                p = p + 1;
            end
        end
        
        function ISI(ax, clusTb, uIdx)
            % Plot ISI histogram of a given unit (using no-duplicate stats)
            
            if nargin < 3
                uIdx = 1;
            end
            if istable(clusTb)
                s = table2struct(clusTb(uIdx,:));
            else
                s = clusTb(uIdx);
            end
            
            x = s.isiEdges;
            x = x(1:end-1) + diff(x)/2;
            y = s.isiCountND;
            
            h = bar(ax, x, y, 'histc');
            h.FaceColor = [0 .7 0];
            h.EdgeColor = 'none';
            
            MPlot.Blocks([0 NP.Param.RP], [0 max(y)], [1 0 0], 'FaceAlpha', .3, 'Parent', ax);
            
            axis(ax, 'tight');
            axis(ax, 'off');
            ax.XLim = [0 0.02];
            ax.Title.String = sprintf('u%i: RPV %.1f, CR %.1f', s.clusId, s.RPV_ND, s.contamND);
        end
        
        function LogISI(ax, clusTb, uIdx)
            % Plot ISI histogram of a given unit (using no-duplicate stats)
            
            if nargin < 3
                uIdx = 1;
            end
            if istable(clusTb)
                s = table2struct(clusTb(uIdx,:));
            else
                s = clusTb(uIdx);
            end
            
            x = s.isiEdges;
            x = x(1:end-1) + diff(x)/2;
            y = s.isiCount;
            
            % h = bar(ax, x, y, 'histc');
            % h.FaceColor = [0 .7 0];
            % h.EdgeColor = 'none';
            area(ax, x, y);
            axis(ax, 'tight');
            axis(ax, 'off');
            ax.XLim = [0 0.02];
            ax.Title.String = sprintf('u%i', s.clusId);
        end
        
        function Span(se, varargin)
            % Plot the period when the units are present
            % 
            %   NP.UnitPlot.Span(se)
            %   NP.UnitPlot.Span(se, uInd)
            %   NP.UnitPlot.Span(se, uInd, figBasePath)
            % 
            
            p = inputParser();
            p.addOptional('uInd', [], @isnumeric);
            p.addOptional('figBasePath', [], @(x) ischar(x) || isstring(x));
            p.parse(varargin{:});
            uInd = p.Results.uInd;
            figBasePath = p.Results.figBasePath;
            
            % Get data
            clusTb = NP.Unit.GetClusTb(se);
            if isempty(uInd)
                uInd = (1:height(clusTb))';
            end
            clusTb = clusTb(uInd,:);
            
            % Vectorize se
            seLite = se.Duplicate({'spikeSpan', 'spikeTime'}, false);
            rt = se.GetReferenceTime;
            seLite.SliceSession([rt(1)-10 rt(end)+10], 'absolute');
            seLite.RemoveEpochs(2);
            
            % Get unit span
            sTb = seLite.GetTable('spikeSpan');
            S = cell2mat(sTb{:,uInd+1});
            t = cell2mat(sTb.time) + seLite.GetReferenceTime;
            y = (1:numel(uInd))';
            
            % Get spike rate
            rTb = seLite.ResampleEventTimes('spikeTime', MMath.BinCenters2Edges(sTb.time{1}));
            R = cell2mat(rTb{:,uInd+1});
            R = R ./ max(R,[],1);
            
            % Combine span and spike rate
            R(~S) = -1;
            
            % Initialize figure
            f = MPlot.Figure(5321); clf
            f.WindowState = 'maximized';
            pause(0.1);
            tl = tiledlayout('flow');
            tl.Padding = 'compact';
            ax = nexttile;
            
            % Plot heatmap
            imagesc(ax, t, y, R');
            colormap(ax, flip(brewermap([], 'RdGy')));
            hold(ax, 'on');
            
            % Add cluster info
            tLabel = repmat(t(end), size(y)) + 1;
            labels = "u"+clusTb.clusId + ", " + clusTb.depth+"um";
            text(ax, tLabel, y, labels, 'VerticalAlignment', 'middle');
            
            % Plot trial marks
            plot([rt rt]', repmat([0 numel(uInd)+1]', [1 numel(rt)]), '-', 'Color', [0 0 0]+.2); hold on
            yLabel = zeros(size(rt));
            text(ax, rt, yLabel, string(1:numel(rt))', 'VerticalAlignment', 'bottom');
            
            % Label recording
            text(ax, tLabel(1), 0, NP.SE.GetID(se), 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'Interpreter', 'none');
            
            MPlot.Axes(ax);
            ax.XLim = [t(1) t(end)+0.05*diff(t([1 end]))];
            ax.YLim = [0 y(end)+1];
            ax.YDir = 'reverse';
            ax.XLabel.String = 'Time (s)';
            ax.YLabel.String = 'Clusters';
            ax.Title.String = " ";
            
            if ~isempty(figBasePath)
                exportgraphics(f, figBasePath + "_cluster_span.png");
            end
        end
        
        function RPV_CDF(clusTb, isSU)
            % 
            
            r = clusTb.RPV;
            rLim = 3;
            r = min(r, rLim); % this will put all units beyond in the last bin
            edges = 0 : .02 : rLim;
            if ~exist('isSU', 'var')
                isSU = false(size(r));
            end
            NP.UnitPlot.MetricCDF(r, edges, isSU);
            ax = MPlot.Axes(gca);
            ax.YLim = [0 numel(r)];
            ax.XTick = 0 : .5 : rLim;
            ax.YTick = linspace(0, numel(r), numel(0:.2:1));
            ax.YTickLabel = 0:.2:1;
            xlabel('Refractory period violation rate (%)');
            ylabel('Fraction of units');
            
            th = NP.Param.maxRPV;
            plot([th th]', ax.YLim', 'Color', 'k');
        end
        
        function ContamCDF(clusTb, isSU)
            % 
            
            r = clusTb.contam;
            rLim = 50;
            edges = 0 : 1 : rLim;
            if ~exist('isSU', 'var')
                isSU = false(size(r));
            end
            NP.UnitPlot.MetricCDF(r, edges, isSU);
            ax = MPlot.Axes(gca);
            ax.YLim = [0 numel(r)];
            ax.XTick = 0 : 5 : rLim;
            ax.YTick = linspace(0, numel(r), numel(0:.2:1));
            ax.YTickLabel = 0:.2:1;
            xlabel('Contamination rate (%)');
            ylabel('Fraction of units');
            
            th = NP.Param.maxContam;
            plot([th th]', ax.YLim', 'Color', 'k');
        end
        
        function SpikeRateCDF(clusTb, isSU)
            % 
            
            r = clusTb.meanActiveRate;
            r = MMath.Bound(r, [0.1 100]); % this will put all units beyond in the last bin
            edges = 10.^(-1 : .02 : 2);
            if ~exist('isSU', 'var')
                isSU = false(size(r));
            end
            NP.UnitPlot.MetricCDF(r, edges, isSU);
            ax = MPlot.Axes(gca);
            ax.XScale = 'log';
            ax.XTick = 10.^(-1:2);
            ax.YLim = [0 numel(r)];
            ax.YTick = linspace(0, numel(r), numel(0:.2:1));
            ax.YTickLabel = 0:.2:1;
            xlabel('Mean spike rate (spk/s)');
            ylabel('Fraction of units');
        end
        
        function NumSpikesCDF(clusTb, isSU)
            % 
            
            N = clusTb.numSpikes;
            N = MMath.Bound(N, [1e2 1e5]); % this will put all units beyond in the last bin
            edges = 10.^(2 : .05 : 5);
            if ~exist('isSU', 'var')
                isSU = false(size(N));
            end
            NP.UnitPlot.MetricCDF(N, edges, isSU);
            ax = MPlot.Axes(gca);
            ax.XScale = 'log';
            ax.XTick = 10.^(2:5);
            ax.YLim = [0 numel(N)];
            ax.YTick = linspace(0, numel(N), numel(0:.2:1));
            ax.YTickLabel = 0:.2:1;
            xlabel('Number of spikes');
            ylabel('Fraction of units');
        end
        
        function SNR_CDF(clusTb, isSU)
            % 
            
            r = clusTb.SNR;
            rLim = 8;
            r = MMath.Bound(r, [0 rLim]); % this will put all units beyond in the last bin
            edges = 0 : .1 : rLim;
            if ~exist('isSU', 'var')
                isSU = false(size(r));
            end
            NP.UnitPlot.MetricCDF(r, edges, isSU);
            ax = MPlot.Axes(gca);
            ax.XTick = 0:1:rLim;
            ax.YLim = [0 numel(r)];
            ax.YTick = linspace(0, numel(r), numel(0:.2:1));
            ax.YTickLabel = 0:.2:1;
            xlabel('Waveform SNR');
            ylabel('Fraction of units');
            
            th = NP.Param.minSNR;
            plot([th th]', ax.YLim', 'Color', 'k');
        end
        
        function MetricCDF(val, edges, isSU)
            % 
            histogram(val, edges, 'Normalization', 'cumcount', ...
                'EdgeColor', 'none', 'FaceColor', [0 0 0], 'FaceAlpha', 0.2); hold on
            histogram(val(isSU), edges, 'Normalization', 'cumcount', ...
                'EdgeColor', 'none', 'FaceColor', NP.UnitPlot.suColor, 'FaceAlpha', NP.UnitPlot.suAlpha); hold on
        end
        
        function ContamVsRPV(clusTb, isSU)
            % 
            
            vLim = 3;
            cLim = 50;
            v = min(clusTb.RPV, vLim);
            c = clusTb.contam;
            
            MPlot.Blocks([0, vLim*1.01], [NP.Param.maxContam, cLim*1.01], [0 0 0], 'FaceAlpha', .1); hold on
            MPlot.Blocks([NP.Param.maxRPV, vLim*1.01], [0, cLim*1.01], [0 0 0], 'FaceAlpha', .1);
            
            if exist('isSU', 'var')
                plot(v(~isSU), c(~isSU), 'k.', 'MarkerSize', 3);
                plot(v(isSU), c(isSU), '.', 'Color', [0 .7 0], 'MarkerSize', 5);
            else
                plot(v, c, 'k.', 'MarkerSize', 3);
            end
            
            ax = MPlot.Axes(gca);
            ax.XTick = 0:.5:vLim;
            ax.YTick = 0:10:cLim;
            axis tight
            xlabel('Refractory period violation rate (%)');
            ylabel('Contamination rate (%)');
        end
        
        function ContamVsSNR(clusTb, isSU)
            % 
            
            rLim = 20;
            cLim = 50;
            r = min(clusTb.SNR, rLim);
            c = clusTb.contam;
            
            MPlot.Blocks([0, rLim*1.01], [NP.Param.maxContam, cLim*1.01], [0 0 0], 'FaceAlpha', .1); hold on
            MPlot.Blocks([0, NP.Param.minSNR*1.01], [-2, cLim*1.01], [0 0 0], 'FaceAlpha', .1);
            
            if exist('isSU', 'var')
                plot(r(~isSU), c(~isSU), 'k.', 'MarkerSize', 3);
                plot(r(isSU), c(isSU), '.', 'Color', [0 .7 0], 'MarkerSize', 5);
            else
                plot(r, c, 'k.', 'MarkerSize', 3);
            end
            
            ax = MPlot.Axes(gca);
            ax.XTick = 0:4:rLim;
            ax.YTick = 0:10:cLim;
            axis tight
            xlabel('SNR');
            ylabel('Contamination rate (%)');
        end
        
        % Physiology
        function WaveformEmbedding(uTb, colorBy, cbName)
            % Plot an interactive embedding map
            % 
            %   WaveformEmbedding(uTb, colorBy)
            % 
            
            if nargin < 2
                colorBy = '';
            end
            
            % Plot dots for brushing
            x = uTb.embedCoords(:,1);
            y = uTb.embedCoords(:,2);
            h = plot(x, y, 'w.'); hold on
            h.UserData.forBrush = true;
            
            % Plot overlays
            ccFun = @lines; % default color scheme
            switch colorBy
                case {'nlrId', 'nmfcId'}
                    val = uTb.(colorBy);
                    ccFun = @LMV.Param.GetTaskPhaseColors;
                    colorBy = 'cluster';
                case 'depth'
                    val = MMath.Bound(uTb.depth, 500:1e3:4500);
                    ccFun = @(n) flip(parula(n), 1);
                case 'recording'
                    val = uTb.recId;
                case 'region'
                    val = uTb.region;
                    val = categorical(val, {'mPrCG', 'vPrCG', 'IFG', 'STG'});
                case {'wfId', 'isiKmId'}
                    val = uTb.(colorBy);
                    ccFun = @(n) brewermap(n, 'Set1');
                    colorBy = 'cluster';
                case ''
                    val = ones(height(uTb), 1);
                otherwise
                    val = uTb.(colorBy);
            end
            groups = unique(val);
            cc = ccFun(numel(groups));
            for i = numel(groups) : -1 : 1
                m = groups(i) == val;
                x = uTb.embedCoords(m,1);
                y = uTb.embedCoords(m,2);
                hh(i) = plot(x, y, '.', 'MarkerSize', 6, 'Color', cc(i,:));
            end
            if numel(groups) > 1
                lgd = legend(hh, string(groups), 'Location', 'eastoutside');
                lgd.Title.String = colorBy;
                lgd.Interpreter = 'none';
            end
            
            x = uTb.embedCoords(:,1);
            y = uTb.embedCoords(:,2);
            ax = gca;
            ax.XLim = [min(x) max(x)] + 0.05*[-1 1]*(max(x) - min(x));
            ax.YLim = [min(y) max(y)] + 0.05*[-1 1]*(max(y) - min(y));
            ax.XAxis.Visible = 'off';
            ax.YAxis.Visible = 'off';
            ax.Title.String = 'Embedding of spatiotemporal spike waveform';
            MPlot.Axes(ax);
            
            % Initialize brush callback
            brushObj = brush(gcf);
            if nargin < 3
                cbName = 'waveform';
            end
            brushObj.ActionPostCallback = @(x,y) NP.UnitPlot.ArrayOnBrush(x, y, uTb, cbName);
        end
        
        function ArrayOnBrush(fig, axesStruct, uTb, cbName)
            % 
            %   ArrayOnBrush(fig, axesStruct, uTb, cbName)
            % 
            
            % Find brushed units
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
            
            b = find(b);
            if numel(b) > 50
                b = randsample(b, 50);
            end
            uTbSub = uTb(b,:);
            
            figure(fig.UserData.cbFig);
            panels = cell(10, 5);
            panelArgs = panels;
            for i = 1 : height(uTbSub)
                panels{i} = cbName;
                panelArgs{i} = {i};
            end
            NP.UnitPlot.Array(uTbSub, panels', 'PanelArgs', panelArgs');
        end
        
        % Response
        function Raster(ax, se, trialInd, tWin, uIdx)
            % Plot raster for a given unit and trials
            % 
            %   Raster(ax, se, trialInd, tWin, uIdx)
            % 
            % Inputs
            %   uIdx            The positional index of the unit in spikeTime table.
            % 
            
            % Slice out spike times data
            if isempty(se.GetReferenceTime('spikeTime'))
                st = se.SliceEventTimes('spikeTime', tWin, trialInd, uIdx, 'Fill', 'none');
            else
                st = se.SliceEventTimes('spikeTime', tWin, trialInd, uIdx, 'Fill', 'bleed');
            end
            
            % Get unit info
            uID = NP.Unit.Ind2ClusId(uIdx, se);
            d2tip = se.userData(1).ksMeta.clusTb.depth(uIdx);
            
            % Plot raster
            hold(ax, 'on');
            [t, y] = NP.UnitPlot.ConvertEventTimesForRasters(st.(1));
            MPlot.PlotPointAsLine(t, y, .6, 'Color', [0 0 0], 'LineWidth', 0.5, 'Parent', ax);
            
            ax.XLim = tWin;
            ax.XLabel.String = 'Time (s)';
            ax.YDir = 'reverse';
            ax.YTick = [];
            ax.YLabel.String = 'Trials';
            title(ax, sprintf('u%i: %ium', uID, d2tip));
        end
        
        function Heatmap(ax, se, trialInd, tWin, uIdx)
            % Plot spike rate for a given unit and trials as a heatmap
            % 
            %   Heatmap(ax, se, trialInd, tWin, uIdx)
            % 
            % Inputs
            %   uIdx            The positional index of the unit in spikeRate table.
            % 
            
            % Get unit info
            uid = NP.Unit.Ind2ClusId(uIdx, se);
            d = se.userData(1).ksMeta.clusTb.depth(uIdx);
            
            if ~ismember('spikeRate', se.tableNames)
                % Compute spikeRates if not present
                t = tWin(1) : 2.5e-3 : tWin(2);
                tEdges = MMath.BinCenters2Edges(t);
                sr = se.ResampleEventTimes('spikeTime', tEdges, trialInd, uIdx);
                sr.(2) = cellfun(@(x) MNeuro.Filter1(x, 1/0.0025, 'gaussian', 0.015), sr.(2), 'Uni', false);
                
            else
                % Slice out spike times data
                sr = se.SliceTimeSeries('spikeRate', tWin, trialInd, uIdx+1, 'Fill', 'none');
                
                % Resample timeseries
                dt = diff(sr.time{1}(1:2));
                t = tWin(1) : dt : tWin(2);
                tEdges = MMath.BinCenters2Edges(t);
                sr = se.ResampleTimeSeries(sr, tEdges);
            end
            R = cat(2, sr.(2){:})';
            
            
            % Mask dropout periods
            if ismember('spikeSpan', se.tableNames)
                ss = se.SliceTimeSeries('spikeSpan', tWin+[-.1 .1], trialInd, uIdx+1, 'Fill', 'none'); % pad a bin of 0.1s
                ss = se.ResampleTimeSeries(ss, tEdges);
                M = cat(2, ss.(2){:})' ~= 1;
                R(M) = NaN;
            end
            
            % Plot raster
            hold(ax, 'on');
            y = 1 : height(sr);
            imagesc(t, y, -R, 'Parent', ax);
            colormap(ax, 'gray');
            
            ax.XLim = tWin;
            ax.XLabel.String = 'Time (s)';
            ax.YDir = 'reverse';
            ax.YLim = [0 height(sr)+1];
            ax.YTick = y(1:2:end);
            ax.YLabel.String = 'Trial';
            ax.TickLength = [0 0];
            grid(ax, 'on');
            title(ax, sprintf('u%i: %ium', uid, d));
        end
        
        function RasterStack(ax, se, trialInd, tWin, unitInd, varargin)
            % Plot a stack of raster for the given units and trials
            % 
            %   NP.UnitPlot.RasterStack(ax, se, trialInd, tWin, unitInd)
            %   NP.UnitPlot.RasterStack(..., 'PETH', false)
            %   NP.UnitPlot.RasterStack(..., 'PETHArgs', {})
            %   NP.UnitPlot.RasterStack(..., 'PeakSpikeRates', [])
            % 
            % Inputs
            %   unitInd         The positional indices of units in the spikeTime table.
            % 
            
            unitInd = unitInd(:);
            
            % Handle optional inputs
            p = inputParser();
            p.KeepUnmatched = true;
            p.addParameter('PETH', false, @islogical);
            p.addParameter('PETHArgs', {}, @iscell);
            p.addParameter('PeakSpikeRates', [], @isnumeric);
            p.parse(varargin{:});
            isPETH = p.Results.PETH;
            pethArgs = p.Results.PETHArgs;
            pkSpkRate = p.Results.PeakSpikeRates(:)'; % make it a row vector
            
            % Remove NaN units
            nu = numel(unitInd);
            isEmp = isnan(unitInd);
            unitInd(isEmp) = [];
            if numel(pkSpkRate) == numel(isEmp)
                pkSpkRate(isEmp) = [];
            end
            if isempty(unitInd)
                disp("No valid unit is specified");
                return
            end
            
            % Slice out spike times data
            st = se.GetTable('spikeTime');
            st = se.SliceEventTimes(st, tWin, trialInd, unitInd, 'Fill', 'none');
            M = ~cellfun(@isempty, st{:,:});
            for i = width(st) : -1 : 1
                stNes{i} = st.(i)(M(:,i));
            end
            
            % Plot PETH stack as background
            if isPETH
                % Get spike rates
                fs = 100;
                tEdges = tWin(1) : 1/fs : tWin(2);
                if ismember('spikeRate', se.tableNames)
                    sr = se.SliceTimeSeries('spikeRate', tWin, trialInd, unitInd+1, 'Fill', 'none');
                    if ~isscalar(unique(cellfun(@numel, sr.time)))
                        sr = se.ResampleTimeSeries(sr, tEdges);
                    end
                else
                    sr = se.ResampleEventTimes(st, tEdges);
                    sr{:,2:end} = cellfun(@(x) MNeuro.Filter1(x, fs, 'gaussian', 0.025), sr{:,2:end}, 'Uni', false);
                end
                t = sr.time{1};
                sr = sr{:,2:end};
                
                % Compute PETH
                [hh, ee] = MNeuro.MeanTimeSeries(sr);
                
                % Normalize
                if isempty(pkSpkRate)
                    [~, ~, pkSpkRate] = MMath.Normalize(hh, 'max');
                end
                hh = hh ./ pkSpkRate;
                ee = ee ./ pkSpkRate;
                ee(:) = NaN;
                
                % Plot PETH stack
                MPlot.PlotHistStack(t, hh, ee, 'Parent', ax, pethArgs{:});
            end
            
            % Plot raster stack
            MPlot.PlotRasterStack(stNes, 'Parent', ax, p.Unmatched);
            
            % Construct unit labels
            clusTb = NP.Unit.GetClusTb(se);
            uDepth = clusTb.depth(unitInd);
            labelArray = ["u"+NP.Unit.Ind2ClusId(unitInd, se)'; uDepth'+"um"];
            tickLabels = sprintf('%s\\newline%s\n', labelArray(:));
            ax.YTick = 1 : sum(~isEmp);
            ax.YTickLabel = strtrim(tickLabels);
            ax.XLim = [min(tWin(:,1)) max(tWin(:,2))];
            ax.YLim = [.5, nu+.5];
            ax.XLabel.String = 'Time (s)';
            % ax.FontSize = 8;
            MPlot.Axes(ax);
        end
        
        function HeatmapStack(ax, se, trialInd, tWin, unitInd, varargin)
            % Plot a stack of raster for the given units and trials
            % 
            %   NP.UnitPlot.HeatmapStack(ax, se, trialInd, tWin, unitInd)
            %   NP.UnitPlot.HeatmapStack(..., 'PeakSpikeRates', [])
            % 
            % Inputs
            %   unitInd         The positional indices of units in the spikeTime table.
            % 
            
            % Remove NaN units
            plotHeight = numel(unitInd);
            unitInd(isnan(unitInd)) = [];
            
            % Slice out an resample spike rate data
            sr = se.SliceTimeSeries('spikeRate', tWin, trialInd, unitInd+1, 'Fill', 'bleed');
            tEdges = MMath.BinCenters2Edges(sr.time{1});
            sr = se.ResampleTimeSeries(sr, tEdges);
            
            n_units = width(sr) - 1;
            n_trials = height(sr);
            
            % Handle custom inputs
            p = inputParser();
            p.KeepUnmatched = true;
            p.addParameter('PeakSpikeRates', [], @(x) numel(x)==n_units);
            p.parse(varargin{:});
            pkSpkRate = p.Results.PeakSpikeRates(:)'; % make it a row vector
            
            % Plot raster stack
            hold(ax, 'on');
            dy_trial = 1 / (n_trials+2);
            t = sr.time{1};
            for i = 1 : n_units
                y_trials = cumsum(repelem(dy_trial, n_trials));
                y_trials = y_trials - mean(y_trials) + i;
                
                R = cat(2, sr.(i+1){:})';
                if isempty(pkSpkRate)
                    pkSpkRate(i) = max(R(:));
                end
                R = R / pkSpkRate(i);
                
                imagesc(ax, t, y_trials, R);
            end
            ax.CLim = [0 2];
            ax.Colormap = brewermap([], 'Greys');
            ax.YLim = [.5, n_units+.5];
            ax.YTick = 1 : n_units;
            ax.YDir = 'reverse';
            ax.XGrid = 'on';
            
            % Construct unit labels
            clusTb = NP.Unit.GetClusTb(se);
            uDepth = clusTb.depth(unitInd);
            labelArray = ["u"+NP.Unit.Ind2ClusId(unitInd, se)'; round(pkSpkRate)+"Hz"; uDepth'+"um"];
            tickLabels = sprintf('%s\\newline%s\\newline%s\n', labelArray(:));
%             labelArray = ["u"+NP.Unit.Ind2ClusId(unitInd, se)';uDepth'+"um"];
%             tickLabels = sprintf('%s\\newline%s\n', labelArray(:));
            ax.YTickLabel = strtrim(tickLabels);
            ax.XLim = [min(tWin(:,1)) max(tWin(:,2))];
            ax.YLim = [.5, plotHeight+.5];
            ax.XLabel.String = 'Time (s)';
            ax.FontSize = 8;
            MPlot.Axes(ax);
        end
        
        function PethRespTests(ce, varargin)
            % Plot PETH and optionally signrank, t-test, Zeta, NLR p-value heatmaps where row are ordered by NMF clusters
            % 
            %   PethRespTests(ce)
            %   PethRespTests(ce, ..., 'signrank', [])
            %   PethRespTests(ce, ..., 'ttest', [])
            %   PethRespTests(ce, ..., 'zeta', [])
            %   PethRespTests(ce, ..., 'nlr', [])
            %   PethRespTests(ce, ..., 'Regions', 'all')
            % 
            
            p = inputParser;
            p.addParameter('signrank', [], @istable);
            p.addParameter('ttest', [], @istable);
            p.addParameter('zeta', [], @istable);
            p.addParameter('nlr', [], @istable);
            p.addParameter('Regions', 'all', @(x) isstring(x) || iscellstr(x) || ischar(x));
            p.parse(varargin{:});
            rTb = p.Results.nlr;
            tTb = p.Results.ttest;
            sTb = p.Results.signrank;
            zTb = p.Results.zeta;
            regions = string(p.Results.Regions);
            
            nPanels = 1;
            nPanels = nPanels + ~isempty(tTb);
            nPanels = nPanels + ~isempty(zTb);
            nPanels = nPanels + ~isempty(rTb);
            nPanels = nPanels + ~isempty(sTb);
            
            tl = tiledlayout(numel(regions), nPanels);
            tl.Padding = "compact";
            
            for i = 1 : numel(regions)
                % Find units in the region
                isRegion = strcmp(regions(i), ce.clusTb.region);
                ceSub = ce.Duplicate;
                ceSub.RemoveUnits(~isRegion);
                
                % PETHs sorted by NMF clusters
                ax = nexttile;
                NP.UnitPlot.PethHeatmap(ceSub, 'SortVars', {'nmfcId', 'nmfcW'}, 'SortOrder', {'ascend', 'descend'});
                titleStr = sprintf("%s (n = %i)", regions(i), ceSub.numResp);
                title(ax, titleStr, 'Interpreter', 'none');
                MPlot.Axes(ax);
                
                % Sort units based on NMF clusters and weighting
                [~, I] = sortrows(ceSub.clusTb, {'nmfcId', 'nmfcW'}, {'ascend', 'descend'});
                
                % t-test p-values
                if ~isempty(sTb)
                    sTbSub = sTb(isRegion,:);
                    sTbSub = sTbSub(I,:);
                    phaseNames = ["atten", "stim", "delay", "init", "prod"];
                    nPhase = numel(phaseNames);
                    S = sign(sTbSub{:,phaseNames+"Resp"});
                    P = -log10(sTbSub{:,phaseNames}) .* S;
                    nStim = size(sTbSub.(phaseNames(1)), 2);
                    
                    ax = plotPvalMap(P);
                    ax.XTick = 1+(0:nPhase-1)*nStim;
                    ax.XTickLabel = phaseNames;
                    ax.Title.String = "signrank -log10(pval)";
                    MPlot.Axes(ax);
                end
                
                % t-test p-values
                if ~isempty(tTb)
                    tTbSub = tTb(isRegion,:);
                    tTbSub = tTbSub(I,:);
                    phaseNames = ["atten", "stim", "delay", "init", "prod"];
                    nPhase = numel(phaseNames);
                    S = sign(tTbSub{:,phaseNames+"Resp"});
                    P = -log10(tTbSub{:,phaseNames}) .* S;
                    nStim = size(tTbSub.(phaseNames(1)), 2);
                    
                    ax = plotPvalMap(P);
                    ax.XTick = 1+(0:nPhase-1)*nStim;
                    ax.XTickLabel = phaseNames;
                    ax.Title.String = "t-test -log10(pval)";
                    MPlot.Axes(ax);
                end
                
                % Zeta p-values
                if ~isempty(zTb)
                    zTbSub = zTb(isRegion,:);
                    zTbSub = zTbSub(I,:);
                    phaseNames = ["atten", "stim", "delay", "init", "prod", "iti"];
                    nPhase = numel(phaseNames);
                    P = -log10([zTbSub{:,phaseNames}]);
                    nStim = size(zTbSub.(phaseNames(1)), 2);
                    
                    ax = plotPvalMap(P);
                    ax.XTick = 1+(0:nPhase-1)*nStim;
                    ax.XTickLabel = phaseNames;
                    ax.Title.String = "Zeta -log10(pval)";
                    MPlot.Axes(ax);
                end
                
                % NLR p-values
                if ~isempty(rTb)
                    rTbSub = rTb(isRegion,:);
                    rTbSub = rTbSub(I,:);
                    phaseNames = ["atten", "stim", "delay", "init", "prod", "iti"];
                    nPhase = numel(phaseNames);
                    P = -log10([rTbSub{:,phaseNames}]);
                    nStim = size(rTbSub.(phaseNames(1)), 2);
                    
                    ax = plotPvalMap(P);
                    ax.XTick = 1+(0:nPhase-1)*nStim;
                    ax.XTickLabel = phaseNames;
                    ax.Title.String = "NLR -log10(pval)";
                    MPlot.Axes(ax);
                end
            end
            
            function ax = plotPvalMap(P)
                ax = nexttile;
                P(isnan(P)) = 0;
                P(abs(P) < -log10(0.05)) = 0;
                imagesc(ax, P);
                colormap(ax, flip(brewermap(255, 'RdBu')));
                colorbar(ax);
                ax.CLim = [log10(0.001) -log10(0.001)];
                ax.YTickLabel = [];
                ax.YTick = [];
            end
        end
        
        function PethHeatmap(ce, varargin)
            % Plot PETHs of as a heatmap, optionally sorting units by given variables
            % 
            %   PethHeatmap(ce)
            %   PethHeatmap(ce, compInd)
            %   PethHeatmap(ce, ..., 'SortVars', [])
            %   PethHeatmap(ce, ..., 'SortOrder', [])
            %   PethHeatmap(ce, ..., 'Ribbon', true)
            %   PethHeatmap(ce, ..., 'Events', true)
            %   PethHeatmap(ce, ..., 'GroupColor', {})
            % 
            
            p = inputParser;
            p.addOptional('compInd', [], @(x) isnumeric(x));
            p.addParameter('SortVars', [], @(x) isstring(x) || iscellstr(x));
            p.addParameter('SortOrder', [], @(x) isstring(x) || iscellstr(x));
            p.addParameter('Ribbon', true, @(x) isscalar(x) && islogical(x));
            p.addParameter('Events', true, @(x) isscalar(x) && islogical(x));
            p.addParameter('GroupColor', {}, @(x) isnumeric(x) || iscell(x));
            p.parse(varargin{:});
            compInd = p.Results.compInd;
            sortVars = cellstr(p.Results.SortVars);
            sortOrder = cellstr(p.Results.SortOrder);
            isRibbon = p.Results.Ribbon;
            isEvents = p.Results.Events;
            gc = p.Results.GroupColor;
            if isnumeric(gc)
                gc = {gc};
            end
            if isscalar(gc) && numel(sortVars)>1
                gc = repmat(gc, size(sortVars));
            end
            
            % Get PETHs with optional sorting
            uTb = ce.clusTb;
            if isempty(sortVars)
                [t, mm] = ce.GetArray('resp');
            else
                if isempty(sortOrder)
                    sortOrder = repmat({'ascend'}, size(sortVars));
                end
                [uTb, I] = sortrows(uTb, sortVars, sortOrder);
                [t, mm] = ce.GetArray('resp', [], I+1);
            end
            
            % Baseline subtraction and normalization
            isBaseline = t < 0;
            mm = mm - mean(mm(isBaseline,:), 1);
            mm = mm ./ max(abs(mm),[],1);
            
            % Heatmap
            ax = gca;
            y = 1 : size(mm,2);
            h = imagesc(ax, t, y, mm'); hold on
            
            % Task event markers
            if isEvents
                colNames = {'cue1', 'stim', 'cue3On', 'prod'};
                NP.TaskBaseClass.PlotEventWindows(ax, ce, 1, [], colNames, 'YRange', y([1 end])+[-1 1]*0.5, ...
                    'Colors', LMV.Param.GetTaskPhaseColors(["none", "stim", "none", "prod"]), 'Style', 'line', 'Alpha', 1);
            end
            
            % Unit group labels
            nSorts = numel(sortVars);
            dt = t(2) - t(1);
            xLims = [t(1)-dt/2, t(end)+dt/2];
            
            for i = nSorts : -1 : 1
                if ~isRibbon
                    continue
                end
                if i > 1 && ~any(sortVars{i} == ["nlrId2", "tId2", "tId3"])
                    continue
                end
                if isempty(compInd)
                    compInd = unique(uTb.(sortVars{i}));
                end
                
                if isempty(gc)
                    cc = lines(numel(compInd));
                else
                    cc = gc{i};
                end
                
                wRib = xLims(1) - .05 - [.2 0];
                xLims(1) = wRib(1);
                [xSeg, ySeg] = MPlot.GroupRibbon(wRib, uTb.(sortVars{i}), cc, ...
                    'Groups', compInd, 'Style', 'patch', 'PlotArgs', {'Parent', ax});
                if i == 1
                    text(ax, xSeg-.1, ySeg, string(compInd), 'HorizontalAlignment', 'right');
                end
            end
            
            if exist('wRib', 'var')
                ax.XLim = xLims;
            end
            
            colormap(ax, flip(brewermap([],'RdBu')));
            ax.CLim = [-1 1];
            ax.XLabel.String = "Aligned time (s)";
            ax.YTick = [];
            ax.YTickLabel = [];
            MPlot.Axes(ax);
            
            h.ButtonDownFcn = @(x,y) NP.UnitPlot.ClickHeatmap(x, y, uTb);
            h.BusyAction = 'cancel';
        end
        
        function ClickHeatmap(hImg, eventData, clusTb)
            % 
            
            % Clears the selection box
%             disp(eventData);
            ax = hImg.Parent;
            if isfield(ax.UserData, 'hRect')
                delete(ax.UserData.hRect);
            end
            if eventData.Button ~= 1
                return
            end
            
            % Find the range of units
            k = round(eventData.IntersectionPoint(2));
            ind = (-4 : 5) + k;
            % ind = (-2 : 2) + k;
            ind = MMath.Bound(ind, [1 height(clusTb)]);
            ind = unique(ind);
            
            % Plot a selection box on the heatmap
            ax.UserData.hRect = rectangle(ax, 'Position', [ax.XLim(1) ind(1)-0.5 diff(ax.XLim) ind(end)-ind(1)+1], 'EdgeColor', 'r');
            
            % Plot raster
            f = MPlot.Figure(123); clf
            % LMV.Overview.SessionFromCache(clusTb.clusId(ind), 'DataSource', 'm2');
            LMV.Overview.SessionFromCache(clusTb.clusId(ind), 'DataSource', 'm2', 'StimIdList', LMV.Param.stimIdList4, 'PlotFun', 'raster');
        end
        
        function RasterOnBrush(fig, axesStruct)
            % 
            %   RasterOnBrush(fig, axesStruct)
            % 
            
            % Find brushed units
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
            
            b = find(b);
            if numel(b) > 50
                b = randsample(b, 50);
            end
            uTb = axesStruct.Axes.UserData.clusTb;
            phaseName = axesStruct.Axes.UserData.taskPhase;
            clusId = uTb.clusId(b);
            
            figure(fig.UserData.cbFig);
            LMV.Overview.SessionFromCache(clusId, 'DataSource', "m2");
            % LMV.Overview.SentencesFromCache(clusId, 'TaskPhase', phaseName);
        end
        
        % Utilities
        function arrId = MakeArrayID(vecId, nRows)
            % Arrange a vector of cluster IDs into an array
            % 
            %   arrId = MakeArrayID(vecId, nRows)
            % 
            nIds = numel(vecId);
            nPages = ceil(nIds / nRows);
            arrId = NaN(nRows, nPages);
            arrId(1:nIds) = vecId;
        end
        
        % Figure specific
        function ProbePETH(ceArray, egStimId)
            % 
            
            tl = tiledlayout(100,1);
            tl.Padding = 'compact';
            
            nRec = numel(ceArray);
            nUnits = arrayfun(@(x) x.numResp, ceArray);
            panelRatio = [8 10];
            panelRatio(3:2+nRec) = floor((100-sum(panelRatio)) / sum(nUnits) * nUnits);
            
            % Get parameters
            ce = ceArray(1); % use example trial from the first recording
            [tt, tv] = ce.GetTable('taskTime', 'taskValue');
            k = find(tv.stimId == egStimId);
            tWin = [tt.cue1On(k) tt.prodMatchOff(k)] + [-1 1]*.2;
            
            % Plot task events
            ax = nexttile([panelRatio(1),1]);
            yRange = [.5 3.5];
            NP.TaskBaseClass.PlotEventWindows(ax, ce, k, tWin, {'cue1', 'cue3'}, 'YRange', yRange);
            NP.TaskBaseClass.PlotTGEHier(ax, ce, k, tWin, {'stim', 'prod'});
            ax.XLabel.String = [];
            ax.XTickLabel = [];
            ax.YLim = yRange;
            
            % Plot spectrogram
            ax = nexttile([panelRatio(2),1]);
            NP.TaskBaseClass.PlotMelSpectrogram(ax, ce, k, [], 'Waveform', 0.3);
            ax.XLim = tWin;
            ax.XLabel.String = [];
            ax.XTickLabel = [];
            
            % Plot neural
            for i = 1 : nRec
                ce = ceArray(i);
                [tt, tv] = ce.GetTable('taskTime', 'taskValue');
                k = find(tv.stimId == egStimId);
                
                ax = nexttile([panelRatio(2+i),1]);
                NP.UnitPlot.SpikeRateHeatmap(ax, ce, k, [], tWin);
                ax.XLabel.String = [];
                ax.XTickLabel = [];
            end
        end
        
        function SpikeRateHeatmap(ax, ce, trialIdx, varargin)
            % Plot spike rate heatmap of a given epoch from the 'resp' table
            % 
            %   PlotSpikeRateHeatmap(ax, se, trialIdx)
            %   PlotSpikeRateHeatmap(ax, se, trialIdx, tWin)
            %   PlotSpikeRateHeatmap(ax, se, trialIdx, unitInd, tWin)
            % 
            % Inputs
            %   tWin                1-by-2 numeric array for the time window to plot.
            % 
            
            p = inputParser;
            p.addOptional('unitInd', [], @isnumeric);
            p.addOptional('tWin', [], @isnumeric);
            p.parse(varargin{:});
            tWin = p.Results.tWin;
            uInd = p.Results.unitInd;
            
            % Get data tables
            [t, R] = ce.GetArray('resp', trialIdx, uInd+1, 'Normalization', 'maxsoft');
            
            % Plot heatmap
            imagesc(ax, t, [], 1-R');
            colormap(ax, 'gray');
            
            % Label depth
            N = ce.numResp;
            if isempty(uInd)
                uInd = (1:N)';
            end
            d = ce.clusTb.depth(uInd) / 1e3;
            dTicks = round(linspace(1, N, 8));
            dLabels = round(d(dTicks), 1);
            
            if ~isempty(tWin)
                ax.XLim = tWin;
            end
            ax.XLabel.String = 'Time (s)';
            ax.YTick = dTicks;
            ax.YTickLabel = dLabels;
            ax.YLabel.String = 'Cortical depth (mm)';
            ax.Title.String = sprintf("%s, n = %i units", NP.SE.GetRegion(ce), N);
            MPlot.Axes(ax);
        end
        
        function MeanWaveformArray(uTb, groupVar)
            % 
            
            % Compute mean waveforms
            groups = unique(uTb.(groupVar));
            k = numel(groups);
            mTb = uTb(1:k,:);
%             mTb.(groupVar) = groups;
            
            for i = 1 : k
                % Collect waveforms of the cluster
                m = uTb.(groupVar) == groups(i);
                W = cat(3, uTb.waveformMed{m});
                
                % Normalize to positive or negative peak
                Wmax = max(abs(W), [], 1);
                Wmax = max(Wmax, [], 2);
                W = W ./ Wmax;
                
                % Compute mean stats
                [mTb.waveformMed{i}, mTb.waveformSD{i}] = MMath.MeanStats(W, 3);
            end
            
            % Plot
            if groupVar == "depthType"
%                 cc = flip(parula(k*3));
%                 cc = cc(2:2:end,:);
                cc = zeros(k,3);
            else
                cc = brewermap(k, 'Set1');
            end
            ccStr = MPlot.Color2Str(round(cc,2));
            
            tl = tiledlayout('flow');
            tl.Padding = 'compact';
            
            for i = 1 : k
                % Scale waveform such that the mean SD equals 0.3
                m = mTb.waveformMed{i};
                e = mTb.waveformSD{i};
                u = mean(e(:)) / 0.3;
                m = m ./ u;
                e = e ./ u;
                mMax = max(abs(m(:)));
                
                % Make coordinates
                [nChan, nTime] = size(m);
                c0 = round(nChan/2);
                c = (1 : nChan) - c0;
                t = 1 : nTime;
                
                % Plot heatmap
                ax = nexttile;
                imagesc(ax, t, c, m); hold on
                
                % Plot trace
                MPlot.ErrorShade(t, m(c0,:), e(c0,:), 'Parent', ax);
                plot(ax, t, m(c0,:), 'Color', [0 0 0]);
                
                axis(ax, 'tight');
                axis(ax, 'off');
                ax.YDir = 'normal';
                ax.YLim = c([1 end]);
                colormap(ax, MPlot.PolarMap());
                ax.CLim = [-1 1] * mMax;
                ax.Title.String = sprintf('\\color[rgb]{%s}%s', ccStr{i}, string(groups(i)));
            end
        end
        
        function MeanIsiArray(uTb, groupVar)
            % 
            
            % Compute mean waveforms
            groups = unique(uTb.(groupVar));
            k = numel(groups);
            mTb = uTb(1:k,:);
%             mTb.(groupVar) = groups;
            
            for i = 1 : k
                % Collect waveforms of the cluster
                m = uTb.(groupVar) == groups(i);
                H = cat(1, uTb.isiCount{m});
                
                % Normalize to peak
                H = H ./ max(H, [], 2);
                
                % Compute mean stats
                [mTb.isiMed{i}, mTb.isiSD{i}, mTb.isiSE{i}] = MMath.MeanStats(H, 1);
            end
            
            % Plot
            if groupVar == "depthType"
%                 cc = flip(parula(k*3));
%                 cc = cc(2:2:end,:);
                cc = zeros(k,3);
            else
                cc = brewermap(k, 'Set1');
            end
            ccStr = MPlot.Color2Str(round(cc,2));
            
            tl = tiledlayout('flow');
            tl.Padding = 'compact';
            
            for i = 1 : k
                % Scale waveform such that the mean SD equals 0.3
                t = MMath.BinEdges2Centers(mTb.isiEdges{i})*1e3;
                m = mTb.isiMed{i};
                e = mTb.isiSD{i};
                
                % Plot
                ax = nexttile;
                MPlot.ErrorShade(t, m, e, 'IsRelative', true, 'Parent', ax);
                hold(ax, 'on');
                plot(ax, t, m, 'Color', cc(i,:));
                ax.XScale = 'log';
                ax.XTick = [2 5 10 20 50 100];
                ax.XLabel.String = "ISI (ms)";
                ax.YTick = [];
                ax.YLabel.String = "Norm. # of ISIs";
                ax.Title.String = sprintf('\\color[rgb]{%s}%s', ccStr{i}, string(groups(i)));
                axis(ax, 'tight');
                MPlot.Axes(ax);
            end
        end
        
    end
    
end
