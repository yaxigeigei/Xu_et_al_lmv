classdef PlotClus
    
    methods(Static)
        
        function Grid(clusTb, uInd, colNames, varargin)
            % Plot various speech features, unit raster, and PETHs aligned in time
            % 
            %   Grid(clusTb, uInd, colNames)
            %   Grid(clusTb, uInd, colNames, rowDist)
            %   Grid(clusTb, uInd, colNames, rowDist, colDist)
            %   Grid(..., 'PanelArgs', {{}})
            % 
            % Inputs
            %   clusTb          Cluster table with models.
            %   uInd            A m-by-n numeric array of clusTb row indices. The shape of this array tells 
            %                   the figure to plot models in m rows and n columns of subplots.
            %   mdlNames        A m-by-n cell array of model variable names in clusTb.
            %   rowDist         A m-element vector specifying integer height (i.e. row) ratios of the panels.
            %   colDist         A n-element vector specifying integer width (i.e. column) ratios of the panels.
            %   'PanelArgs'     A cell array of structs. The size of this array should match that of panels. 
            %                   Each struct contains the additional arguments for plotting. Therefore, the 
            %                   recepient function must support struct expansion of Parameter-Value pairs.
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
                rowDist = ones(size(uInd,1), 1);
            end
            if isempty(colDist)
                colDist = ones(size(uInd,2), 1);
            end
            tl = tiledlayout(sum(rowDist), sum(colDist));
            tl.Padding = 'tight';
            
            % Match arguments size to panels
            if isscalar(panelArgs)
                panelArgs = repmat({{}}, size(uInd));
            end
            
            % Plot models
            for i = 1 : size(uInd, 1)
                for j = 1 : size(uInd, 2)
                    % Create Axes
                    ntArgs = MPlot.FindTileInd(rowDist, colDist, i, j);
                    ax = nexttile(ntArgs{:});
                    hold(ax, 'on');
                    
                    % Plot a panel
                    ui = uInd(i,j);
                    vn = colNames{i,j};
                    args = panelArgs{i,j};
                    if isnan(ui) || isempty(vn)
                        cla(ax);
                        axis(ax, 'off');
                        continue
                    end
                    switch vn
                        case 'phase'
                            mdl = clusTb.(vn){ui};
                            LMV.RF.PlotWeights(mdl, 'ClusInfo', clusTb(ui,:), 'Parent', ax);
                        otherwise
                            % fprintf("Row %i column %i: '%s' is not a valid model name.\n", i, j, mn);
                            mdl = clusTb.(vn){ui};
                            LMV.TRF.PlotWeights2(mdl, 'ClusInfo', clusTb(ui,:), 'Parent', ax, args{:});
                    end
                end
            end
        end
        
        
        
        
    end
    
end