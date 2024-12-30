classdef Recon
    %RECON Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        iElvisPath = fullfile(getenv('NP_ROOT'), 'code', 'third_party_matlab', 'iElvis');
    end
    
    methods(Static)
        function cfgOut = PlotSitesOnAvgBrain(src, varargin)
            % Plot recording sites on average brain
            % 
            %   cfgOut = PlotSitesOnAvgBrain(src)
            %   cfgOut = PlotSitesOnAvgBrain(seArray)
            %   cfgOut = PlotSitesOnAvgBrain(recMeta)
            %   cfgOut = PlotSitesOnAvgBrain(..., userCfg)
            %   cfgOut = PlotSitesOnAvgBrain(..., 'ShowNames', false)
            %   cfgOut = PlotSitesOnAvgBrain(..., 'Color', [])
            % 
            
            p = inputParser;
            p.addOptional('UserCfg', [], @isstruct);
            p.addParameter('ShowNames', false, @(x) isscalar(x) && islogical(x));
            p.addParameter('Color', [], @(x) isempty(x) || x=='k');
            p.parse(varargin{:});
            userCfg = p.Results.UserCfg;
            isShowNames = p.Results.ShowNames;
            colorOpt = p.Results.Color;
            
            % Standardize metadata to table
            if isa(src, 'MSessionExplorer')
                src = arrayfun(@(x) x.userData.recMeta, src);
            end
            if isstruct(src)
                src = struct2table(src, 'AsArray', true);
                src.recId = string(src.Subject) + "_" + string(src.Block);
            end
            if ~iscell(src.Coords)
                src.Coords = mat2cell(src.Coords, ones(height(src),1));
            end
            
            % Set missing coordinates to empty
            for i = 1 : height(src)
                if isempty(src.Coords{i}) || ischar(src.Coords{i})
                    src.Coords{i} = NaN(0,3);
                end
            end
            
            % Rename STG subregions as just STG
            m = contains(src.Region, 'STG');
            src.Region(m) = {'STG'};
            
            % Make a matrix of site coordinates
            m = ~cellfun(@isempty, src.Coords);
            if all(~m)
                error("No coordinates are available to plot.");
            end
            coords = cat(1, src.Coords{m});
            coords(:,4) = 1;
            coords(:,1) = -abs(coords(:,1));
            
            % Get site names
            siteNames = src.recId(m);
            siteNames = replace(siteNames, '_', '-');
            
            % Get site colors
            if isempty(colorOpt)
                regionColors = NP.Param.GetRegionColors(src.Region(m));
            elseif colorOpt == 'k'
                regionColors = zeros(sum(m), 3);
            end
            
            % Configure plotting
            cfg = [];
            cfg.fsurfSubDir = fullfile(NP.Recon.iElvisPath, 'freesurfer_brains');
            
            cfg.figId = gcf;
            cfg.view = 'l';
            % cfg.elecCoord = 'n';
            cfg.elecCoord = coords;
            cfg.title = "";
            
            cfg.elecShape = 'sphere';
            cfg.elecColors = regionColors;
            cfg.elecColorScale = [0 1 1];
            cfg.showLabels = 'n';
            cfg.elecUnits = [];
            cfg.elecCbar = 'n';
            cfg.elecNames = siteNames;
            cfg.elecSize = 2;
            
            % Overwrite default cfg by user values
            if exist('userCfg', 'var') && isstruct(userCfg)
                fn = fieldnames(userCfg);
                for i = 1 : numel(fn)
                    cfg.(fn{i}) = userCfg.(fn{i});
                end
            end
            
            % Plot
            cfgOut = plotPialSurf("cvs_avg35_inMNI152", cfg);
            
            if isShowNames
                coords = cfgOut.electrodeCoords;
                text(coords(:,1)-20, coords(:,2)-3, coords(:,3), src.recId(m), 'FontSize', 6, 'FontWeight', 'bold', 'Interpreter', 'none');
            end
            
            % Make the brain brighter than default
            lightLevel = 1; % originally 1
            hl = light('Color',[1 1 1]*0.15*lightLevel);
            lightangle(hl, -90, -20);
            hl = light('Color',[1 1 1]*0.2*lightLevel);
            lightangle(hl, -30, 70);
            hl = light('Color',[1 1 1]*0.2*lightLevel);
            lightangle(hl, -150, 70);
        end
        
        
    end
end

