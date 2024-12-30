classdef CE
    
    methods(Static)
        function [ceArray, filePaths] = LoadSession(varargin)
            % Load NP.CodingExplorer objects from files
            % 
            %   [ceArray, filePaths] = LoadSession()
            %   [ceArray, filePaths] = LoadSession(filePaths)
            %   [ceArray, filePaths] = LoadSession(..., 'UserFunc', @(se) se)
            % 
            
            % Handle user inputs
            p = inputParser();
            p.addOptional('filePaths', '', @(x) ischar(x) || iscellstr(x) || isstring(x) || isempty(x));
            p.addParameter('UserFunc', @(x) x, @(x) isa(x, 'function_handle'));
            p.parse(varargin{:});
            filePaths = p.Results.filePaths;
            userFunc = p.Results.UserFunc;
            
            if isempty(filePaths)
                filePaths = MBrowse.Files();
            else
                filePaths = cellstr(filePaths);
            end
            if isempty(filePaths)
                ceArray = [];
                return;
            end
            
            % Preallocation
            ceArray(numel(filePaths),1) = NP.CodingExplorer();
            
            for i = 1 : numel(filePaths)
                load(filePaths{i}, 'ce');
                disp(['Loading ce ' num2str(i) ' - ' NP.SE.GetID(ce)]);
                userFunc(ce);
                ceArray(i) = ce;
            end
        end
        
    end
end

