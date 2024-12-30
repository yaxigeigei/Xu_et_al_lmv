classdef PlotBaseClass
    
    methods(Static)
        function [evtTimes, evtY, yTickVal] = ConvertEventTimesForRasters(evtTimes)
            
            if isnumeric(evtTimes)
                evtTimes = num2cell(evtTimes, 2);
            end
            evtTimes = evtTimes(:);
            
            evtNumber = cellfun(@length, evtTimes);
            yTickVal = (1:length(evtNumber))';
            evtY = arrayfun(@(n,i) ones(n,1)*i, evtNumber, yTickVal, 'Uni', false);
            
            evtTimes = cell2mat(evtTimes);
            evtY = cell2mat(evtY);
        end
        
        function MovableType(data, panels, varargin)
            % Plot various speech features, unit raster, and PETHs aligned in time
            % 
            %   MovableType(data, panels)
            %   MovableType(data, panels, rowDist)
            %   MovableType(data, panels, rowDist, colDist)
            %   MovableType(..., 'PanelArgs', {{}})
            % 
            % Inputs
            %   data            The common data to plot from in each panel.
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
            
            % Match argument size to panels
            if numel(panelArgs) == 1
                panelArgs = repmat(panelArgs, size(panels));
            end
            
            % Plot pannels
            for i = 1 : size(panels, 1)
                for j = 1 : size(panels, 2)
                    % Create Axes
                    g = mat2cell(zeros(nRow, nCol), rowDist, colDist);
                    g{i,j}(:) = 1;
                    g = cell2mat(g);
                    gInd = find(g');
                    ax = MPlot.Axes(nRow, nCol, gInd);
                    hold(ax, 'on');
                    
                    % Plot a panel
                    p = panels{i,j};
                    args = panelArgs{i,j};
                    if isempty(p)
                        % plot nothing
                    elseif ischar(p)
                        switch p
                            case 'plot_type_1'
                                % plot something
                            case 'plot_type_2'
                                % plot something
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
        
    end
end

