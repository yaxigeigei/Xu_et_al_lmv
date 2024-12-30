%% Paths

add_se_paths;

%% Select recordings

% Find all raw se files
seSearch = MBrowse.Dir2Table(fullfile(NP.Data.GetDatastoreRoot, "se", "NP*_B*_se_raw.mat"));
if height(seSearch) == 0
    disp("No se_raw file found");
    return
end

% Interactive selection of one or more recordings
[ind, isSelect] = listdlg('ListString', seSearch.name, ...
    'PromptString', 'Select one or more raw se to process', ...
    'SelectionMode', 'multiple', ...
    'ListSize', [200 400]);
if ~isSelect
    return
end

% % Manual selection
% recIds = ["NP41_B1", "NP44_B2", "NP44_B3", "NP45_B1", "NP45_B2", "NP46_B1", "NP47_B3"];
% seRawNames = recIds + "_se_raw.mat";
% ind = ismember(seSearch.name, seRawNames);

seSearch = seSearch(ind,:);

%% 

% Read spreadsheet metadata
subjectMeta = NP.Preproc.ImportMetaSheet([], "Subjects");
recMeta = NP.Preproc.ImportMetaSheet([], "Recordings");

% Iterate through recordings
for i = 1 : height(seSearch)
    % Load raw se
    disp("Loading raw se...");
    load(fullfile(seSearch.folder{i}, seSearch.name{i}), 'se');
    
    % Process metadata
    NP.SE.UpdateConventions(se); % dealing with backward compatibility
    NP.SE.AddSubjectMeta(se, subjectMeta);
    NP.SE.AddRecMeta(se, recMeta);
    NP.Unit.MergeMeta(se);
    NP.Unit.SetUniqueClusId(se);
    NP.Unit.ConvertDepth(se);
    
    % Cleaning
    if NP.SE.GetID(se) == "NP53_B1"
        NP.SE.TrimRecording(se, 'TrimWin', [92 1171]);
    else
        NP.SE.TrimRecording(se, 'TrimWin', 'sorting');
    end
    NP.SE.FixNonMonotonicTimestamps(se); % fix non-monotonic timestamps
    NP.Unit.RemoveDuplicateSpikes(se); % remove double counted spikes
    
    % Process tasks
    if ismember('taskTime', se.tableNames)
        % Slice recording to trials
        tasks = se.userData.recMeta.Task;
        tasks = strsplit(tasks, {' ', ','});
        tasks = unique(tasks, 'stable');
        NP.SE.SliceToTrials(se, 'taskNames', tasks);
        NP.SE.AddTaskValueTable(se, 'taskNames', tasks);
    end
    
    % Save processed se
    saveName = NP.SE.GetID(se) + "_se.mat";
    save(fullfile(seSearch.folder{i}, saveName), 'se', '-v7.3');
end
