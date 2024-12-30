%% Cache epoch-sliced M1 and M2 se

for m = ["m1" "m2"]
    % Load recordings
    if m == "m0"
        seSearch = LMV.Data.FindSource([]);
    else
        seSearch = MBrowse.Dir2Table(fullfile(LMV.Data.GetAnalysisDir, 'data', 'se_'+m, '*_se_'+m+'.mat'));
    end
    sePaths = fullfile(seSearch.folder, seSearch.name);
    seEpPaths = cellstr(strrep(sePaths, '.mat', '_ep.mat'));
    
    % Reslice and clean trials
    for i = 1 : numel(sePaths)
        % Check computed
        if exist(seEpPaths{i}, 'file')
            fprintf("\nFound resliced file at %s\nTo recompute, please delete the existing file.\n", seEpPaths{i});
            continue
        else
            fprintf("\nReslicing %s\n", seSearch.name{i});
        end
        
        % Load se
        se = NP.SE.LoadSession(sePaths{i});
        
        % Slice trials to include a period before and after
        tOffset = [-1 1]; % in sec
        se = NP.SE.BleedTrials(se, tOffset);
        
        % Remove trials with bad performance
        if m ~= "m0"
            isBad = LMV.SE.IsBadTrials(se, -Inf, 0.5);
            se.RemoveEpochs(isBad);
        end
        
        % Save resliced se
        save(seEpPaths{i}, 'se', '-v7.3');
    end
end
