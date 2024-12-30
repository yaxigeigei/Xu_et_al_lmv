classdef ST
    
    methods(Static)
        
        function [stF, N] = TriggerFeatures(R, F, varargin)
            % Get spike triggered features
            % 
            %   [stF, N] = TriggerFeatures(R, F, varargin)
            % 
            
            % Get feature vectors associated with spikes
            isSpk = any(R, 2);
            stF = F(isSpk,:);
            N = R(isSpk,:);
            
            % Replicate feature vectors by spike counts
            stF = repelem(stF, N, 1);
        end
        
        function PlotFeatureUMAP(Z, stF, featNames, t)
            % 
            
            x = Z(:,1);
            y = Z(:,2);
            
            hCB = plot(x, y, 'k.');
            hCB.UserData.forBrush = true;
            
            ax = gca;
            ax.XLim = [min(x) max(x)] + 0.05*[-1 1]*(max(x) - min(x));
            ax.YLim = [min(y) max(y)] + 0.05*[-1 1]*(max(y) - min(y));
            ax.XAxis.Visible = 'off';
            ax.YAxis.Visible = 'off';
%             ax.Title.String = 'Embedding of unit responses';
            MPlot.Axes(ax);
            
            % Initialize brush callback
            ax.UserData.stF = stF;
            ax.UserData.featNames = featNames;
            ax.UserData.t = t;
            brushObj = brush(gcf);
            brushObj.ActionPostCallback = @LMV.ST.PlotBrushedFeatures;
        end
        
        function PlotBrushedFeatures(fig, axesStruct)
            % 
            %   PlotBrushedFeatures(fig, axesStruct)
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
                    'Name', 'Selected Features', ...
                    'NumberTitle', 'off', ...
                    'Menubar', 'none', ...
                    'Toolbar', 'figure', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
            end
            
            ax = axesStruct.Axes;
            stF = ax.UserData.stF;
            featNames = ax.UserData.featNames;
            t = ax.UserData.t;
            
            stF = mean(stF(b,:), 1);
            
            figure(fig.UserData.cbFig);
            LMV.ST.PlotFeatureHeatmap(stF, featNames, t);
            ax = gca;
            ax.Title.String = sprintf("n = %i/%i", sum(b), numel(b));
        end
        
        function PlotFeatureHeatmap(F, featNames, t)
            % 
            
            if isvector(F)
                F = reshape(F, numel(featNames), []);
            end
            
            if nargin < 3
                t = [];
            end
            
            imagesc(t, [], F);
            ax = gca;
            if startsWith(featNames{1}, ["mic", "speaker"])
                ax.YTick = [];
                ax.YDir = "normal";
            else
                ax.YTick = 1 : numel(featNames);
                ax.YTickLabel = featNames;
            end
            ax.XLabel.String = "Time from spike (s)";
            MPlot.Axes(ax);
        end
        
    end
    
end