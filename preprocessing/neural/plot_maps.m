%% Paths

add_np_paths;

preprocRoot = NP.Data.GetPreprocRoot;

%% Recording info

if ~exist('recId', 'var')
    subjId = "NP47";
    blockId = "B3";
    recId = subjId + '_' + blockId;
end

%% Select sorts

% Find rez.mat files
rezPattern = fullfile(preprocRoot, recId, "kilosort", "**", "rez.mat");
rezSearch = MBrowse.Dir2Table(rezPattern);
if height(rezSearch) < 1
    error('No file matches %s', rezPattern);
end
rezFiles = fullfile(rezSearch.folder, rezSearch.name);

% Allow user to select a subset of rez.mat files
[k, isSelect] = listdlg('ListString', rezFiles, ...
    'PromptString', 'Select a rez.mat files', ...
    'SelectionMode', 'multiple', ...
    'ListSize', [800 100]);
if ~isSelect
    return
end

%% Load data

% Load rez.mat files
rezTb = rezSearch(k,:);
rezTb.path = rezFiles(k);
rezTb.rez = cellfun(@(x) load(x, 'rez'), rezTb.path, 'Uni', false);
rezTb.rez = cellfun(@(x) x.rez, rezTb.rez, 'Uni', false);

% Order full maps first
[~, ksDirNames] = cellfun(@(x) fileparts(string(x)), rezTb.folder);
rezTb.ksDirName = ksDirNames;
rezTb.baseName = erase(ksDirNames, "_map");
rezTb.isMapOnly = contains(ksDirNames, "map");
rezTb = sortrows(rezTb, ["baseName", "isMapOnly"], ["ascend", "descend"]);

% Find the widest time range
tWin = [Inf -Inf];
for i = 1 : height(rezTb)
    ops = rezTb.rez{i}.ops;
    tWin = [min(tWin(1), ops.tstart), max(tWin(2), ops.tstart+ops.sampsToRead)];
end
tWin = tWin / ops.fs;

%% Plot spike maps

% Initialize figure
h = min(1000, 400 * numel(k));
f = MPlot.Figure('Name', 'Spike Map', 'Position', [100 100 1200 h], 'Visible', 'on');
clf(f);

for i = 1 : height(rezTb)
    % Plot spike map
    ax = subplot(numel(k), 1, i);
    MKilosort2.PlotDriftMap(ax, rezTb.rez{i}, ~rezTb.isMapOnly(i));

    titlePattern = recId + "/.*?(?=\/rez\.mat)";
    titleStr = regexp(rezTb.path{i}, titlePattern, 'match', 'once');
    ax.Title.String = titleStr;
    ax.Title.Interpreter = 'none';

    ax.XLim = tWin;
    if i ~= numel(k)
        % ax.XLabel.String = [];
        ax.XTickLabel = [];
    end
    
    % Cache spike map data for python
    % cache_py_spike_map(ksDir, rezTb.rez{i});
end

return
%% 

% Save figure
figDir = fileparts(fileparts(rezTb.folder{1}));
timestamp = datestr(now, 'yyyymmdd_HHMMSS'); % add timestamp to the filename
saveas(f, fullfile(figDir, sprintf("%s_spike_maps_%s.png", recId, timestamp)));
