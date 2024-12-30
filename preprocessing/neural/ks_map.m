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

%% Initial run to get the spike map

for i = 1 : numel(k)
    % Get Kilosort output folder path
    apFile = apFiles{k(i)};
    [probeDir, apName] = fileparts(apFile);
    [runDir, probeDirName] = fileparts(probeDir);
    [~, runDirName] = fileparts(runDir);
    ksFolder = fullfile(preprocRoot, recId, "kilosort", runDirName, probeDirName+"_map");

    % Run Kilosort
    rez = MKilosort2.Sort(apFile, ksFolder, 'DriftMapOnly', true, ksArgs);

    % Cache a MTracer session with spike and LFP map
    mt = MTracerVM();
    mt.LoadChannelMap(fullfile(ksFolder, "chanMap.mat"));
    mt.LoadLFP(strrep(apFile, '.ap.bin', '.lf.bin'));
    mt.LoadRez(ksFolder);
    mtName = ['MTracer_session_' mt.recId '_spike-lfp-map'];
    save(fullfile(ksFolder, 'mtracer_cache', mtName), 'mt');

    % Cache spike map data for python
    cache_py_spike_map(ksFolder, rez);
end
