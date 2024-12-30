classdef Data
    
    properties(Constant)
        metaSheetName = "np_meta_sheet.xlsx";
    end
    
    methods(Static)
        function srcTb = FindSource(kw, varargin)
            % Find source data for the specified analysis
            % 
            %   [fileTb, metaTb] = FindSource(preset)
            %   [fileTb, metaTb] = FindSource(recId)
            %   [fileTb, metaTb] = FindSource(..., 'PathPatterns', fullfile(NP.Data.GetProjectRoot, "se", "*", "NP*_B*_se.mat"))
            %   [fileTb, metaTb] = FindSource(..., 'ExcludePatterns', fullfile(NP.Data.GetProjectRoot, "se", "auto", "NP*_B*_se.mat"))
            %   [fileTb, metaTb] = FindSource(..., 'Tasks', [])
            %   [fileTb, metaTb] = FindSource(..., 'Regions', [])
            % 
            % Inputs
            %   preset, recId       1) The name of a preset, such as one for a specific analysis.
            %                       2) Recording ID(s).
            %                       3) Include all files if pass [].
            %                       Note that preset can overwrite 'PathPatterns' argument.
            %   'PathPatterns'      One or more patterns of file path to search with.
            %   'ExcludePatterns'   One or more patterns of file path to exclude files from the match using 'PathPatterns'.
            %   'Tasks'             Narrow down to specific task(s). This can overwrite preset tasks.
            %   'Regions'           Narrow down to specific regions(s). This can overwrite preset regions.
            % 
            
            r = fullfile(NP.Data.GetProjectRoot, "se");
            
            % Handle inputs
            p = inputParser;
            p.addParameter('Tasks', [], @(x) ischar(x) || isstring(x) || iscellstr(x));
            p.addParameter('Regions', [], @(x) ischar(x) || isstring(x) || iscellstr(x));
            p.addParameter('PathPatterns', fullfile(r, "*", "NP*_B*_se.mat"), @(x) ischar(x) || isstring(x));
            p.addParameter('ExcludePatterns', fullfile(r, "auto", "NP*_B*_se.mat"), @(x) ischar(x) || isstring(x));
            p.parse(varargin{:});
            tasks = string(p.Results.Tasks);
            regions = string(p.Results.Regions);
            pathPatterns = string(p.Results.PathPatterns);
            exPatterns = string(p.Results.ExcludePatterns);
            
            % Check for recording ID
            kw = lower(string(kw));
            if isempty(kw)
                % No selection
            elseif ~isempty(regexpi(kw(1), '^NP[0-9]{2}_B[0-9]{1,2}$', 'once'))
                % Recognize recording ID
                pathPatterns = fullfile(r, "*", upper(kw)+"_se.mat");
            end
            
            % Find se files that match the pattern
            fileTbs = arrayfun(@(x) NP.Data.Dir2Table(x), pathPatterns, 'Uni', false);
            fileTb = cat(1, fileTbs{:});
            
            % File se files that match the exclusion patterns
            exTbs = arrayfun(@(x) NP.Data.Dir2Table(x), exPatterns, 'Uni', false);
            exTb = cat(1, exTbs{:});
            if ~isempty(exTb)
                isEx = ismember(fullfile(fileTb.folder, fileTb.name), fullfile(exTb.folder, exTb.name));
                fileTb(isEx,:) = [];
            end
            
            % Load metadata from spreadsheet
            metaTb = NP.Preproc.ImportMetaSheet(fileTb.recId);
            
            % Filter recordings
            if isempty(tasks)
                hasTask = true(size(metaTb.Task));
            else
                taskList = arrayfun(@(x) strsplit(string(x), {' ', ','}), string(metaTb.Task), 'Uni', false);
                hasTask = cellfun(@(x) any(lower(tasks(:))==lower(x), 'all'), taskList);
            end
            if isempty(regions) || any(strcmpi(regions, "all"))
                isRegion = true(size(metaTb.Region));
            else
                isRegion = contains(metaTb.Region, regions);
            end
            srcTb = [fileTb metaTb];
            srcTb = srcTb(hasTask & isRegion, :);
        end
        
        function tb = Dir2Table(searchPattern)
            % A wrapper of MBrowse.Dir2Table that also adds recording info
            
            tb = MBrowse.Dir2Table(searchPattern);
            
            varNames = {'recId', 'folder', 'name', 'path'};
            if isempty(tb)
                tb = array2table(zeros(0, numel(varNames)), 'VariableNames', varNames);
            else
                tb.path = fullfile(tb.folder, tb.name);
                tb.recId = cellfun(@NP.SE.GetID, tb.name, 'Uni', false);
                tb = tb(:,varNames);
            end
        end
        
        function p = GetDatastoreRoot()
            % Return the root directory of Neuropixels data on server in a machine dependent way
            
            serverDir = '/data_store2/neuropixels';
            mountDir = fullfile('S:', serverDir);
            
            if exist(serverDir, 'dir')
                p = serverDir;
            elseif exist(mountDir, 'dir')
                p = mountDir;
            else
                error('Cannot find server data root');
            end
        end
        
        function p = GetRawRoot()
            % Return the root directory of raw recordings in a machine dependent way
            p = fullfile(NP.Data.GetDatastoreRoot, 'raw');
            if ~exist(p, 'dir')
                error('Cannot find raw data root');
            end
        end
        
        function p = GetPreprocRoot()
            % Return the root directory of preprocessing in a machine dependent way
            
            serverPath = fullfile(NP.Data.GetDatastoreRoot, 'preproc');
            npRoot = getenv('NP_ROOT');
            localPath = fullfile(npRoot, 'preproc');
            
            if exist(serverPath, 'dir')
                p = serverPath;
            elseif ~isempty(npRoot) && exist(localPath, 'dir')
                p = localPath;
            else
                error('Cannot find preprocessing root');
            end
        end
        
        function p = GetProjectRoot()
            % Return the environment variable for the NP project root directory
            p = getenv('NP_ROOT');
        end
        
        function p = GetAnalysisRoot()
            % Return the root directory of analysis in a machine dependent way
            
            npRoot = getenv('NP_ROOT');
            localPath = fullfile(npRoot, 'analysis');
            serverPath = '/userdata/dxu/project_np/analysis';
            mountPath = fullfile('S:', serverPath);
            
            if ~isempty(npRoot) && exist(localPath, 'dir')
                p = localPath;
            elseif exist(serverPath, 'dir')
                p = serverPath;
            elseif exist(mountPath, 'dir')
                p = mountPath;
            else
                error('Cannot find analysis root');
            end
        end
        
        function cmdStr = GetJobCmd(scriptPath, varargin)
            % Get the command for running a script with a server job
            % 
            %   cmdStr = NP.Data.GetJobCmd(scriptPath)
            %   cmdStr = NP.Data.GetJobCmd(scriptPath, isRun)
            %   cmdStr = NP.Data.GetJobCmd(..., 'Queue', 'pia-batch')
            %   cmdStr = NP.Data.GetJobCmd(..., 'NumCPU', 8)
            %   cmdStr = NP.Data.GetJobCmd(..., 'RAM', [])
            %   cmdStr = NP.Data.GetJobCmd(..., 'Version', 'R2022a')
            % 
            
            p = inputParser;
            p.addOptional('isRun', false, @islogical);
            p.addParameter('Queue', 'pia-batch.q', @(x) ismember(x, {'pia-batch.q', 'skull-batch.q', 'spirit-batch'}));
            p.addParameter('NumCPU', 8, @(x) x > 0);
            p.addParameter('RAM', [], @(x) x > 0);
            p.addParameter('Version', 'R2022a', @(x) ismember(x, {'R2019a', 'R2022a'}));
            p.parse(varargin{:});
            isRun = p.Results.isRun;
            qName = p.Results.Queue;
            nCPU = p.Results.NumCPU;
            ram = p.Results.RAM;
            mlVer = p.Results.Version;
            
            if isempty(ram)
                ram = nCPU * 16;
            end
            
            scriptPath = strrep(scriptPath, '\', '/');
            scriptPath = string(scriptPath);
            [folder, name] = fileparts(scriptPath);
            logFile = strjoin([folder, name+".txt"], '/');
            
            cmdStr = sprintf("submit_job -q %s -c %d -m %d -o %s -n %s -x /data_store2/MATLAB/%s/bin/matlab %s", ...
                qName, nCPU, ram, logFile, name, mlVer, scriptPath);
            
            if isRun
                system(cmdStr);
            end
        end
        
        function SubmitJob(cmdStr)
            % 
            system(cmdStr);
        end
    end
end
