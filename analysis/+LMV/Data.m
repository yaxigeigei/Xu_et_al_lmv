classdef Data
    
    methods(Static)
        function srcTb = FindSource(kw, varargin)
            % Find source data for the specified analysis
            % 
            %   [fileTb, metaTb] = FindSource(preset)
            %   [fileTb, metaTb] = FindSource(recId)
            %   [fileTb, metaTb] = FindSource(..., 'PathPatterns', fullfile(NP.Data.GetProjectRoot, "se", "*", "NP*_B*_se.mat"))
            %   [fileTb, metaTb] = FindSource(..., 'ExcludePatterns', [])
            %   [fileTb, metaTb] = FindSource(..., 'Tasks', [])
            %   [fileTb, metaTb] = FindSource(..., 'Regions', [])
            % 
            % Inputs
            %   preset              1) The name of a preset, such as one for a specific analysis.
            %                       2) Recording ID(s).
            %                       3) Include all files if pass [].
            %                       Note that preset can overwrite 'PathPatterns' argument.
            %   'PathPatterns'      One or more patterns of file path to search with.
            %   'ExcludePatterns'   One or more patterns of file path to exclude files from the match using 'PathPatterns'.
            %   'Tasks'             Narrow down to specific task(s). This can overwrite preset tasks.
            %   'Regions'           Narrow down to specific regions(s). This can overwrite preset regions.
            % 
            
            r = fullfile(NP.Data.GetProjectRoot, "se");
            
            exPatterns = fullfile(r, "*", [ ...
                "NP80_B1_se.mat"; ... % IFG pars triangularis
                "NP106_B1_se.mat"; ... % later recordings
                "NP113_B1_se.mat"; ... % later recordings
                "NP119_B1_se.mat"; ... % later recordings
                ]);
            
            % Handle inputs
            p = inputParser;
            p.addParameter('Tasks', "lmv", @(x) ischar(x) || isstring(x) || iscellstr(x));
            p.addParameter('Regions', ["mPrCG", "vPrCG", "IFG", "STG"], @(x) ischar(x) || isstring(x) || iscellstr(x));
            p.addParameter('PathPatterns', fullfile(r, "*", "NP*_B*_se.mat"), @(x) ischar(x) || isstring(x));
            p.addParameter('ExcludePatterns', exPatterns, @(x) ischar(x) || isstring(x));
            p.parse(varargin{:});
            tasks = string(p.Results.Tasks);
            regions = string(p.Results.Regions);
            pathPatterns = string(p.Results.PathPatterns);
            exPatterns = string(p.Results.ExcludePatterns);
            
            % Get preset filters
            kw = lower(string(kw));
            if isempty(kw)
                % Not using preset
                
            elseif ~isempty(regexpi(kw(1), '^NP[0-9]{2}_B[0-9]{1,2}$', 'once'))
                % Use recording ID
                pathPatterns = fullfile(r, "*", upper(kw)+"_se.mat");
                
            elseif kw == "probe_peth"
                % Example recordings used to plot probe PETHs
                pathPatterns = fullfile(r, "*", [ ...
                    "NP41_B1_se.mat"; ... % mPrCG
                    % "NP52_B1_se.mat"; ... % vPrCG
                    "NP69_B1_se.mat"; ... % IFG
                    "NP43_B1_se.mat"; ... % STG
                    ]);
                
            elseif kw == "phone_stats"
                % Example recording used in 
                pathPatterns = [ ...
                    fullfile(r, "mPrCG", "NP41_B1_se.mat"); ...
%                     fullfile(r, "mPrCG", "NP44_B2_se.mat"); ...
                    ];
            end
            
            % Find se files that match the pattern
            srcTb = NP.Data.FindSource([], 'Tasks', tasks, 'Regions', regions, 'PathPatterns', pathPatterns, 'ExcludePatterns', exPatterns);
        end
        
        function anaDir = GetAnalysisDir(varargin)
            % Return the root directory of LMV analysis
            % 
            %   anaRoot = GetAnalysisDir()
            %   anaDir = GetAnalysisDir(subDir, subsubDir, ...)
            % 
            anaDir = fullfile(NP.Data.GetProjectRoot, "analysis_lmv", varargin{:});
            if ~exist(anaDir, 'dir')
                mkdir(anaDir);
            end
        end
        
    end
end
