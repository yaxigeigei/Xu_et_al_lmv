%% Paths

add_se_paths;

%% Select recordings

% Interactive selection of one or more recordings
ppRoot = NP.Data.GetPreprocRoot;
recSearch = MBrowse.Dir2Table(fullfile(ppRoot, "NP*_B*"));
[ind, isSelect] = listdlg('ListString', recSearch.name, ...
    'PromptString', 'Select one or more recordings to make se', ...
    'SelectionMode', 'multiple', ...
    'ListSize', [200 600]);
if ~isSelect
    return
end
recIds = string(recSearch.name(ind));

% % Manual selection
% recIds = ["NP41_B1", "NP44_B2", "NP44_B3", "NP45_B1", "NP45_B2", "NP46_B1", "NP47_B3"];

%% 

for i = 1 : numel(recIds)
    % Initialize se
    sePath = fullfile(NP.Data.GetDatastoreRoot, 'se', recIds(i) + "_se_raw.mat");
    if exist(sePath, 'file')
        % Continue from an existing se
        fprintf("\nLoad existing se from:\n%s\n", sePath);
        load(sePath, 'se');
        M = NP.SEMaker(se);
        M.isOverwrite = true;
    else
        % Initialize a new maker
        fprintf("\nCreate new se at:\n%s\n", sePath);
        M = NP.SEMaker(recIds(i));
    end
    
    % NIDQ data
    M.AddNI();
    M.ReplaceMicAudio();
    
    % LFP data
    M.AddLFP();
    
    % Sorting results
    M.AddAPMeta();
    M.AddSortingResult('SelectionPrompt', true);
    
    % Task, reading, and speech events
    if ismember("taskTime", M.se.tableNames)
        se.RemoveTable("taskTime");
    end
    M.AddTaskEvents();
    M.AddReadingEvents();
    M.AddSpeechEvents();
    
    % Audio feature extraction results
    M.AddIntensity();
    M.AddPitch();
    M.AddArtics();
    
    % Video feature extraction results
    M.AddLandmark();
    
    % Fix known issues in the recording
    M.PosthocFix();
    
    % Save SE
    se = M.se;
    save(sePath, 'se', '-v7.3');
    disp('Finished saving se');
end
