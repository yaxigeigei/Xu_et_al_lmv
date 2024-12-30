classdef Embed
    
    methods(Static)
        function BatchUMAP(ce)
            % Compute UMAP embedding from PETHs
            % 
            %   BatchUMAP(ce)
            % 
            % Input
            %   ce              One or a vector of NP.CodingExplorer object(s).
            % 
            
            for i = 1 : numel(ce)
                % Print progress
                fprintf("\nCompute UMAP embedding for %s\n", NP.SE.GetID(ce(i)));
                
                % Get PETHs
                [~, X] = ce(i).GetArray('resp', 'Normalization', 'minmaxsoft');
                
                % Embedding
                rng(61);
                ce(i).clusTb.embedCoords = run_umap(X', 'metric', 'correlation', 'verbose', 'text', ...
                    'min_dist', 0.1, 'n_neighbors', 15, 'randomize', false);
            end
        end
        
        function PlotScatter(ce, colorBy)
            % Plot an interactive embedding map
            % 
            %   PlotScatter(ce, colorBy)
            % 
            
            if nargin < 2
                colorBy = 'cluster';
            end
            
            ax = nexttile;
            
            % Plot dots for brushing
            uTb = ce.clusTb;
            x = uTb.embedCoords(:,1);
            y = uTb.embedCoords(:,2);
            h = plot(x, y, 'w.'); hold on
            h.UserData.forBrush = true;
            
            % Plot overlays
            ccFun = @lines; % default color scheme
            switch colorBy
                case 'waveform'
                    val = uTb.wfId;
                    colorBy = 'waveform';
                case 'isi'
                    val = uTb.isiKmId;
                    colorBy = 'burstiness';
                case 'depth'
                    val = MMath.Bound(uTb.depth, 1500:1e3:4500);
                case 'recording'
                    val = uTb.recId;
                case 'region'
                    val = uTb.region;
                    val = categorical(val, {'mPrCG', 'vPrCG', 'IFG', 'STG'});
                case {'nlrId', 'nmfcId', 'zetaId'}
                    val = uTb.(colorBy);
                    ccFun = @LMV.Param.GetTaskPhaseColors;
                    colorBy = 'responsiveness';
                case {'tId1'}
                    val = uTb.(colorBy);
                    groups = unique(val);
                    ccFun = @(x) LMV.Param.GetTaskPhaseColors(string(groups));
                    colorBy = 'responsiveness';
                case 'hcId'
                    val = uTb.(colorBy);
                    groups = ce.userData.hcIdList;
                otherwise
                    warning("'%s' is not a supported color group option.", colorBy);
                    val = uTb.(colorBy);
            end
            
            if ~exist('groups', 'var')
                groups = unique(val);
            end
            
            cc = ccFun(numel(groups));
            for i = numel(groups) : -1 : 1
                m = groups(i) == val;
                if ~any(m)
                    continue
                end
                x = uTb.embedCoords(m,1);
                y = uTb.embedCoords(m,2);
                hh(i) = plot(x, y, '.', 'MarkerSize', 16, 'Color', cc(i,:));
            end
            lgd = legend(hh, string(groups), 'Location', 'eastoutside');
            lgd.Title.String = colorBy;
            lgd.Interpreter = 'none';
            
            x = uTb.embedCoords(:,1);
            y = uTb.embedCoords(:,2);
            ax.XLim = [min(x) max(x)] + 0.05*[-1 1]*(max(x) - min(x));
            ax.YLim = [min(y) max(y)] + 0.05*[-1 1]*(max(y) - min(y));
            ax.XAxis.Visible = 'off';
            ax.YAxis.Visible = 'off';
%             ax.Title.String = 'Embedding of unit responses';
            MPlot.Axes(ax);
            
            % Initialize brush callback
            ax.UserData.ce = ce;
            brushObj = brush(gcf);
            brushObj.ActionPostCallback = @NP.Embed.PlotBrushedUnits;
        end
        
        function PlotBrushedUnits(fig, axesStruct)
            % 
            %   PlotBrushedUnits(fig, axesStruct)
            % 
            
            % Find brushed units
            ce = axesStruct.Axes.UserData.ce;
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
            
            ax = fig.UserData.cbFig.Children;
            if isempty(ax) || ~isvalid(ax)
                ax = MPlot.Axes(fig.UserData.cbFig);
            end
            hold(ax, 'off');
            NP.Embed.PlotPETH(ax, ce, b);
        end
        
        function PlotPETH(ax, ce, unitInd)
            % Plot PETH overlay, cluster means, and histograms of unit properties
            % 
            %   PlotPETH(ax, ce, unitInd)
            % 
            
            uTb = ce.clusTb;
            cidList = unique(uTb.nmfcId);
%             cc = lines(numel(cidList));
            cc = LMV.Param.GetTaskPhaseColors(numel(cidList));
            
            if ~islogical(unitInd)
                unitInd = MMath.Ind2Logical(unitInd, ce.numResp);
            end
            unitInd = unitInd(:);
            
            cla(ax);
            hold(ax, 'on');
            
            colNames = {'cue1', 'stim', 'cue3', 'prod'};
            yRange = [0.95 1];
            NP.TaskBaseClass.PlotEventWindows(ax, ce, 1, [], colNames, 'YRange', yRange);
            
            for i = 1 : numel(cidList)
                isClus = uTb.nmfcId == cidList(i);
                mask = isClus & unitInd;
                if ~any(mask)
                    continue
                end
                [t, mm] = ce.GetArray('resp', [], find(mask)+1, 'Normalization', 'minmaxsoft');
                m = MMath.MeanStats(mm, 2);
%                 [~, m] = ce.ComputeMean(isClus, 'Normalization', 'minmaxsoft');
                plot(ax, t, m, 'Color', cc(i,:), 'LineWidth', 1.5);
                plot(ax, t, mm, 'Color', [cc(i,:) MPlot.AlphaForOverlap(sum(unitInd))]);
            end
            
            ax.XLim = t([1 end]);
            ax.XLabel.String = "Time (s)";
        end
        
    end
    
end

