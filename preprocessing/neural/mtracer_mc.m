%% Paths

add_np_paths;

preprocRoot = NP.Data.GetPreprocRoot;

%% Recording info

if ~exist('recId', 'var')
    subjId = "NP47";
    blockId = "B3";
    recId = subjId + '_' + blockId;
end

if ~exist('ksArgs', 'var')
    warning("Using default parameters");
    ksArgs.ChannelMapFile = [];
    ksArgs.ConfigFunc = @MKilosort2.Config384;
end

%% Select recordings

% Find ap.bin files
apPattern = fullfile(preprocRoot, recId, "sglx", "*", "*", "*.ap.bin");
apSearch = MBrowse.Dir2Table(apPattern);
if height(apSearch) < 1
    error('No file matches %s', apPattern);
end
apFiles = fullfile(apSearch.folder, apSearch.name);

% Allow user to select a subset of ap.bin files
[k, isSelect] = listdlg('ListString', apFiles, ...
    'PromptString', 'Select an ap.bin file', ...
    'SelectionMode', 'multiple', ...
    'ListSize', [800 100]);
if ~isSelect
    return
end

%% Generate motion corrected files

for i = 1 : numel(k)
    % Get mc output folder
    apFile = apFiles{k(i)};
    [probeDir, apName] = fileparts(apFile);
    [runDir, probeDirName] = fileparts(probeDir);
    [sglxDir, runDirName] = fileparts(runDir);
    mcRunDirName = "mc_" + erase(runDirName, {'catgt_', 'mc_'});
    mcProbeDir = fullfile(sglxDir, mcRunDirName, probeDirName);
    
    % Find interpolant
    f2Pattern = fullfile(preprocRoot, recId, "kilosort", "*"+runDirName, probeDirName+"_map", "mc", "F2*.mat");
    f2Search = MBrowse.Dir2Table(f2Pattern);
    f2File = fullfile(f2Search.folder{1}, f2Search.name{1});
    
    % Find channel map
    ksOutDir = fileparts(f2Search.folder{1});
    cmFile = fullfile(ksOutDir, "chanMap.mat");
    
    % Motion correction
    imec = Neuropixel.ImecDataset(apFile, 'channelMap', cmFile);
    load(f2File, 'F2');
    todoList = {'ap', 'lf'};
    F2.CorrectBinary(imec, mcProbeDir, todoList);
    
    % Copy the original apFile as ap.meta and lf.meta files
    oldapName = dir(apFile);
    oldapName = oldapName.name;
    oldlfName = strrep(oldapName, '.ap.bin', '.lf.bin');
    oldapMetaName = strrep(oldapName, '.ap.bin', '.ap.meta');
    oldlfMetaName = strrep(oldapName, '.ap.bin', '.lf.meta');
    
    % Find the new ap.bin file
    newapName = dir(fullfile(mcProbeDir, '*.ap.bin'));
    newapName = newapName.name;
    newlfName = strrep(newapName, '.ap.bin', '.lf.bin');
    newapMetaName = strrep(newapName, '.ap.bin', '.ap.meta');
    newlfMetaName = strrep(newapName, '.ap.bin', '.lf.meta');

    % Copy the original ap.meta and lf.meta files to the new location (with the new name)
    copyfile(fullfile(probeDir, oldapMetaName), fullfile(mcProbeDir, newapMetaName)); % Only new name
    fprintf('Copied %s to %s\n', fullfile(probeDir, oldapMetaName), fullfile(mcProbeDir, newapMetaName));
    copyfile(fullfile(probeDir, oldlfMetaName), fullfile(mcProbeDir, newlfMetaName)); % Only new name
    fprintf('Copied %s to %s\n', fullfile(probeDir, oldlfMetaName), fullfile(mcProbeDir, newlfMetaName));

    % Rename the new ap.bin and lf.bin files to include the index number
    movefile(fullfile(mcProbeDir, newapMetaName), fullfile(mcProbeDir, oldapMetaName)); % Rename with index
    fprintf('Renamed %s to %s\n', fullfile(mcProbeDir, newapMetaName), fullfile(mcProbeDir, oldapMetaName));
    movefile(fullfile(mcProbeDir, newlfMetaName), fullfile(mcProbeDir, oldlfMetaName)); % Rename with index
    fprintf('Renamed %s to %s\n', fullfile(mcProbeDir, newlfMetaName), fullfile(mcProbeDir, oldlfMetaName));

    movefile(fullfile(mcProbeDir, newapName), fullfile(mcProbeDir, oldapName)); % Rename with index
    fprintf('Renamed %s to %s\n', fullfile(mcProbeDir, newapName), fullfile(mcProbeDir, oldapName));
    movefile(fullfile(mcProbeDir, newlfName), fullfile(mcProbeDir, oldlfName)); % Rename with index
    fprintf('Renamed %s to %s\n', fullfile(mcProbeDir, newlfName), fullfile(mcProbeDir, oldlfName));
end
