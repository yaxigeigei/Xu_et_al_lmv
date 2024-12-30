classdef UnitCache
    
    methods(Static)
        % Unit cache
        function CreateCache(seSrc, cacheDir)
            % Extract and save data of each unit to file
            % 
            %   CreateUnitCache(sePaths)
            %   CreateUnitCache(sePaths, cacheDir)
            %   CreateUnitCache(seArray, cacheDir)
            % 
            % Input
            %   sePaths         One or more paths of saved se objects.
            %   seArray         One or more se objects.
            %   cacheDir        The directory path to save cache in. The seArray input requires cacheDir 
            %                   to be provided. If using the sePaths input, by default, cache will be 
            %                   saved in a subfolder named 'unit_cache' in the folder where the se files live. 
            %                   This default location can be overwritten if cacheDir is provided. 
            % Cache
            %   The cache file is named as "u<clusId>.mat", and the following variables are saved or 
            %   updated (if the cache file already exists).
            %   unitInfo        A row of se.userData.ksMeta.clusTb table converted to struct.
            %   tt              A 'taskTime' table where each row is an example trial of a unique stim.
            %   tv              A 'taskValue' table where each row is an example trial of a unique stim.
            %   st              Spike times saved in a #stim-element cell array where each element is 
            %                   again a #trial-element cell array.
            % 
            % See also NP.Unit.AddData2UnitCache
            
            for i = 1 : numel(seSrc)
                % Get se object
                if iscellstr(seSrc) || isstring(seSrc)
                    % From file
                    seSrc = cellstr(seSrc);
                    se = NP.SE.LoadSession(seSrc{i});
                    seDir = fileparts(seSrc{i});
                else
                    % Use input
                    se = seSrc(i);
                    assert(~isempty(cacheDir), "When using seArray as input, cacheDir must be provided to specify the location to save the cache.");
                end
%                 fprintf("\nCache unit data from %s\n", NP.SE.GetID(se));
                
                % Determine folder to save cache
                if ~exist('cacheDir', 'var') || isempty(cacheDir)
                    cacheDir = fullfile(seDir, "unit_cache");
                end
                if ~exist(cacheDir, 'dir')
                    mkdir(cacheDir);
                end
                
                % Compute sentence PETH
                senTb = NP.TaskBaseClass.SplitBySentence(se);
                [ce, senTb] = LMV.SE.ComputeSentencePETH(senTb);
                
                % Extract and save unit data
                clusTb = NP.Unit.GetClusTb(se);
                [tt, tv] = ce.GetTable('taskTime', 'taskValue');
                stSen = arrayfun(@(x) x.GetTable('spikeTime'), senTb.se, 'Uni', false);
                stUnit = cell(height(clusTb), 1);
                for u = 1 : height(clusTb)
                    stUnit{u} = cellfun(@(x) x.(u), stSen, 'Uni', false);
                end
                
                s.unitInfo = table2struct(clusTb);
                s.tt = tt;
                s.tv = tv;
                s.st = stUnit;
                NP.Unit.AddData2UnitCache(cacheDir, clusTb.clusId, s);
            end
        end
        
        function AddData2UnitCache(cacheDir, clusId, dataStruct)
            % Add data variables to unit cache
            % 
            %   AddData2UnitCache(cacheDir, clusId, dataStruct)
            % 
            % Input
            %   cacheDir        The directory path for the cache.
            %   clusId          A vector of cluster IDs.
            %   dataStruct      A struct of data to be cached. Field names determines the variable names 
            %                   used in unit cache. Two forms of field values are supported:
            %                   1) When the number of elements of a field value equals the number of units, 
            %                      each element is saved to its corresponding unit's cache.
            %                   2) When unequal, the field value is copied as a whole and saved to every 
            %                      unit caches.
            % Cache
            %   The cache file is named as "u<clusId>.mat" with the specified variables saved or updated 
            %   (if the cache file already exists).
            % 
            % See also NP.Unit.CreateUnitCache
            
            fn = fieldnames(dataStruct);
            
            for u = 1 : numel(clusId)
                % Put data for this unit in a struct
                s = struct;
                for n = 1 : numel(fn)
                    d = dataStruct.(fn{n});
                    if numel(d) == numel(clusId)
                        % Take an element from the variable
                        if iscell(d)
                            s.(fn{n}) = d{u};
                        else
                            s.(fn{n}) = d(u);
                        end
                    else
                        % Use the variable as a whole
                        s.(fn{n}) = d;
                    end
                end
                
                % Add variables to cache
                uCachePath = fullfile(cacheDir, "u"+clusId(u)+".mat");
                if exist(uCachePath, 'file')
                    fprintf("u%i: update existing unit cache\n", clusId(u));
                    save(uCachePath, '-struct', 's', '-append');
                else
                    fprintf("u%i: create unit cache\n", clusId(u));
                    save(uCachePath, '-struct', 's');
                end
            end
        end
        
        function sUnit = LoadUnitCache(unitId, varargin)
            % Plot rasters of stim aligned to speech labels
            % 
            %   sUnit = LoadUnitCache(unitId)
            %   sUnit = LoadUnitCache(unitId, dataSource)
            % 
            
            % Process inputs
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('dataSource', 'm1', @(x) ischar(x) || isstring(x) || isempty(x));
            p.parse(varargin{:});
            dataSrc = string(p.Results.dataSource);
            
            % Load cache
            if any(strcmpi(dataSrc, ["m1", "m2"]))
                srcFolder = fullfile(NP.Data.GetAnalysisRoot, "data", "se_"+lower(dataSrc), "unit_cache");
            else
                srcFolder = dataSrc;
            end
            srcFiles = fullfile(srcFolder, "u"+unitId+".mat");
            sUnit = arrayfun(@load, srcFiles, 'Uni', false);
        end
        
        function clusTb2 = AlignClusTb(clusTb1, clusTb2, isRmRddCols)
            % Sort clusTb2 in the same order as clusTb1 and remove existing columns from clusTb2
            % 
            %   clusTb2 = AlignClusTb(clusTb1, clusTb2)
            %   clusTb2 = AlignClusTb(clusTb1, clusTb2, isRmRddCols)
            % 
            
            % Remove units from clusTb2 that do not exist in clusTb1
            m = ismember(clusTb2.clusId, clusTb1.clusId);
            clusTb2(~m,:) = [];
            
            % Sort clusTb2
            [~, I] = MMath.SortLike(clusTb2.clusId, clusTb1.clusId);
            clusTb2 = clusTb2(I,:);
            
            % Remove redundant table columns
            if nargin < 3
                isRmRddCols = false;
            end
            if isRmRddCols
                m = ismember(clusTb2.Properties.VariableNames, clusTb1.Properties.VariableNames);
                clusTb2(:,m) = [];
            end
        end
        
        
        
    end
    
end