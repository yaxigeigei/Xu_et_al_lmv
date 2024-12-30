%% Movies showing mirror activity

anaDir = LMV.Data.GetAnalysisDir('movies', 'mirror_m11');

%% Load data

% Find unit info
uu = LMV.Linker.GetExampleUnitInfo("mirror_movie");

% Load M11 se
recIds = unique(cellfun(@(x) string(NP.SE.GetID(x.clusId)), uu));
sePaths = fullfile(LMV.Data.GetAnalysisDir, "linker", "se_m11", recIds+"_se.mat");
seArr = NP.SE.LoadSession(sePaths);

% Load M11 ce (sem discounted)
cePaths = fullfile(LMV.Data.GetAnalysisDir, "linker", "computed_ce", recIds+"_ce.mat");
ceArr = NP.CE.LoadSession(cePaths);

%% Construct ce objects, one for stim, one for prod, each concatenated across recordings

% Find the common set of sentences
stimIdArr = arrayfun(@(x) x.GetTable("taskValue").stimId, ceArr, 'Uni', false);
stimIdList = intersect(stimIdArr{1}, stimIdArr{2});

% % Specify the sentence of interest
% stimIdList = [ ...
%     "fsjk1_si2285", ... % the girl nodded ...
%     "mbbr0_si2315" ... % junior, what on earth ...
%     ];

% Split ce by stim and prod phase for each recording
ce2Arr = repmat(NP.CodingExplorer, [numel(ceArr), 2]);

for i = 1 : numel(ceArr)
    % Remove sentences that are not shared
    ce = ceArr(i).Duplicate;
    tv = ce.GetTable("taskValue");
    isSen = ismember(tv.stimId, stimIdList);
    ce.RemoveEpochs(~isSen);
    
    % % Reference to sentence onset times
    % ce.AlignTime('prodMatchOn', 'taskTime');
    
    % Normalize responses
    tbNames = ["resp", "sd", "sem"];
    for k = 1 : numel(tbNames)
        tb = ce.GetTable(tbNames(k));
        R = tb{:,2:end};
        R = cellfun(@(x) fillmissing(x, "nearest"), R, 'Uni', false);
        if k == 1
            [~, ~, kk] = MMath.Normalize(cell2mat(R), 'max');
            K = num2cell(repmat(kk, [size(R,1) 1]));
        end
        R = cellfun(@(x,k) x/k, R, K, 'Uni', false);
        tb{:,2:end} = R;
        ce.SetTable(tbNames(k), tb);
    end
    
    % Split stim and prod epochs
    ceTb = ce.SplitConditions("phase", "taskValue");
    ceTb = flip(ceTb);
    ce2 = ceTb.se;
    
    % Apply a consistent order of epochs
    for j = 1 : numel(ce2)
        tv = ce2(j).GetTable("taskValue");
        [~, I] = MMath.SortLike(tv.stimId, stimIdList);
        ce2(j).SortEpochs(I);
    end
    
    % Resample unit responses to the same length
    if i == 1
        tEdges = cellfun(@(x) MMath.BinCenters2Edges(x), ce2(1).GetTable("resp").time, 'Uni', false);
    end
    for j = 1 : numel(ce2)
        for k = 1 : numel(tbNames)
            tb = ce2(j).ResampleTimeSeries(tbNames(k), tEdges, 'Extrapolation', 'nearest');
            ce2(j).SetTable(tbNames(k), tb);
        end
    end
    
    ce2Arr(i,:) = ce2;
end

% Concatenate units
ce2 = [ce2Arr(:,1).CatUnits, ce2Arr(:,2).CatUnits];

% Label the phase of each ce
for p = 1 : numel(ce2)
    ce2(p).userData.phaseName = ceTb.phase(p);
end

%% Add phase info and spike times to ce objects

exampleClusIds = cellfun(@(x) x.clusId, uu);

for p = 1 : numel(ce2)
    % Make clusTb for example units
    clusTb = ce2(p).clusTb;
    [~, I] = MMath.SortLike(clusTb.clusId, exampleClusIds);
    clusTb = clusTb(I,:);
    
    % Add spike times to clusTb for each unit
    for u = 1 : height(clusTb)
        for r = 1 : numel(seArr)
            % Check if the unit is in this recording
            isUnit = clusTb.clusId(u) == NP.Unit.GetClusTb(seArr(r)).clusId;
            if ~any(isUnit)
                continue
            end
            
            % Extract spike times
            [tv, st] = seArr(r).GetTable("taskValue", "spikeTime");
            for s = 1 : numel(stimIdList)
                isSen = tv.stimId==stimIdList(s) & tv.phase==ce2(p).userData.phaseName;
                clusTb.spikeTime{u,s} = st{isSen, isUnit};
            end
        end
    end
    
    % Save data to ce
    ce2(p).userData.expClusTb = clusTb;
end

%% Generate overlay videos

% Isolate ce objects
ceStim = ce2(1);
ceProd = ce2(2);
ceUni = ceStim.Duplicate;
ceUni.userData.phaseName = "uni";

% Set up MP
if ~exist('mp', 'var')
    mp = MPlotter();
end
mp.RemovePlot();

figNum = 1;
rowDist = [2 9];
colDist = 1;

spInd = MPlot.FindSubplotInd(rowDist, colDist, 1, 1);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @cla, '', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.LinkerTrialPhone(a, b), 'ceUni', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.LinkerTrialEvents(a, b), 'ceUni', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.SentenceTitle(a, b), 'ceUni', 'epoch');
mp.AddPlot(figNum, spStr);

spInd = MPlot.FindSubplotInd(rowDist, colDist, 2, 1);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @(x,y) cla(x), '', 'epoch');
mp.AddPlot(figNum, spStr, @LMV.MP.LinkerResponseStack, 'ceStim', 'epoch');
mp.AddPlot(figNum, spStr, @LMV.MP.LinkerResponseStack, 'ceProd', 'epoch');
mp.AddPlot(figNum, spStr);

mp.RefreshAll();

f = figure(1);
% f.Position(1:2) = [0 50];
f.Position(3:4) = [800 600];
MPlot.Paperize(f, 'FontSize', 8);

mp.timeLimits = [0 3];

%% 

% Select sentences
senInd = [1 2 6 8:11];

[tt, tv] = ce2(1).GetTable("taskTime", "taskValue");

for i = senInd
    % Get time window
    tWin = [tt.prodMatchOn(i) tt.prodMatchOff(i)] + [-0.3 0.3];
    
    % Set epoch and time
    mp.timeLimits = tWin;
    mp.epoch = i;
    mp.time = tWin(1);
    
    MPlot.Paperize(f, 'FontSize', 8);
    
    % Make movie
    vidFile = fullfile(anaDir, sprintf("overlay_sent%i_%s.mp4", i, tv.stimId(i)));
    fps = 60;
    if ~isfile(vidFile)
        frames = mp.MakeVideo(f, 1/fps, 'FrameRate', fps, 'FilePath', vidFile);
    end
    
    % Save speaker audio
    stimFile = sprintf("overlay_sent%i_%s_stim.wav", i, tv.stimId(i));
    k = find(seArr(1).GetTable("taskValue").trialNum == tv.tempTrialNum(i), 1);
    if ~isfile(stimFile)
        NP.Audio.WriteTrial(seArr(1), k, tWin, 'Channel', "speaker1", 'FolderPath', anaDir, 'FileName', stimFile);
    end
end

%% Make response modulated tone

% Specify tone frequencies
F = [261.63 329.63]; % Oct4: middle C, E
F = F.*[3 2]';

% Find example units
[~, I] = MMath.SortLike(ce2(p).clusTb.clusId, exampleClusIds);

% Compute spike rate modulated tone
A = cell(ce2(1).numEpochs, numel(ce2));
for p = 1 : numel(ce2)
    [T, R] = ce2(p).GetArray('resp', [], I+1, 'DimCat', 0);
    [~, E] = ce2(p).GetArray('sem', [], I+1, 'DimCat', 0);
    R = cellfun(@(x,y) max(x-y,0), R, E, 'Uni', false);
    [T, A(:,p)] = cellfun(@(t,r) LMV.Resp.ComputeModulatedTone(t,r,F(:,p)), T, R, 'Uni', false);
end
A = cellfun(@(x,y) x+y, A(:,1), A(:,2), 'Uni', false);

% Add results to ni table
ni = ce2(1).GetTable("ni");
ni.mod = cellfun(@interp1, T, A, ni.time, 'Uni', false);
ce2(1).SetTable("ni", ni);

% Save spike modulated tones
[tt, tv] = ce2(1).GetTable("taskTime", "taskValue");
for i = senInd
    % Get time window
    tWin = [tt.prodMatchOn(i) tt.prodMatchOff(i)] + [-0.3 0.3];
    
    modFile = sprintf("overlay_sent%i_%s.wav", i, tv.stimId(i));
    % if ~isfile(modFile)
        NP.Audio.WriteTrial(ce2(1), i, tWin, 'Channel', "mod", 'FolderPath', anaDir, 'FileName', modFile);
    % end
end

%% Make stim plots

if ~exist('mp', 'var')
    mp = MPlotter();
end
mp.RemovePlot();

figNum = 1;
rowDist = [2 9];
colDist = 1;

spInd = MPlot.FindSubplotInd(rowDist, colDist, 1, 1);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @cla, '', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.LinkerTrialPhone(a, b), 'ceStim', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.LinkerTrialEvents(a, b), 'ceStim', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.SentenceTitle(a, b), 'ceStim', 'epoch');
mp.AddPlot(figNum, spStr);

spInd = MPlot.FindSubplotInd(rowDist, colDist, 2, 1);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @(x,y) cla(x), '', 'epoch');
mp.AddPlot(figNum, spStr, @LMV.MP.LinkerResponseStack, 'ceStim', 'epoch');
mp.AddPlot(figNum, spStr);

mp.RefreshAll();

f = figure(1);
f.Position(3:4) = [800 600];
MPlot.Paperize(f, 'FontSize', 8);

mp.timeLimits = [0 3];

%% 

[tt, tv] = ceStim.GetTable("taskTime", "taskValue");

for i = senInd
    % Get time window
    tWin = [tt.prodMatchOn(i) tt.prodMatchOff(i)] + [-0.3 0.3];
    
    % Set epoch and time
    mp.timeLimits = tWin;
    mp.epoch = i;
    mp.time = NaN; % hide time marker
    
    MPlot.Paperize(f, 'FontSize', 8);
    
    % Save figure
    imgFile = fullfile(anaDir, sprintf("stim_sent%i_%s.pdf ", i, tv.stimId(i)));
    if ~isfile(imgFile)
        exportgraphics(mp.plotTable.figureObj{1}, imgFile, ContentType="vector", BackgroundColor="none");
    end
end

%% Make prod plots

if ~exist('mp', 'var')
    mp = MPlotter();
end
mp.RemovePlot();

figNum = 1;
rowDist = [2 9];
colDist = 1;

spInd = MPlot.FindSubplotInd(rowDist, colDist, 1, 1);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @cla, '', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.LinkerTrialPhone(a, b), 'ceProd', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.LinkerTrialEvents(a, b), 'ceProd', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.SentenceTitle(a, b), 'ceProd', 'epoch');
mp.AddPlot(figNum, spStr);

spInd = MPlot.FindSubplotInd(rowDist, colDist, 2, 1);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @(x,y) cla(x), '', 'epoch');
mp.AddPlot(figNum, spStr, @LMV.MP.LinkerResponseStack, 'ceProd', 'epoch');
mp.AddPlot(figNum, spStr);

mp.RefreshAll();

f = figure(1);
% f.Position(1:2) = [0 50];
f.Position(3:4) = [800 600];
MPlot.Paperize(f, 'FontSize', 8);

mp.timeLimits = [0 3];

%% 

[tt, tv] = ceProd.GetTable("taskTime", "taskValue");

for i = senInd
    % Get time window
    tWin = [tt.prodMatchOn(i) tt.prodMatchOff(i)] + [-0.3 0.3];
    
    % Set epoch and time
    mp.timeLimits = tWin;
    mp.epoch = i;
    mp.time = NaN; % hide time marker
    
    MPlot.Paperize(f, 'FontSize', 8);
    
    % Save figure
    imgFile = fullfile(anaDir, sprintf("prod_sent%i_%s.pdf ", i, tv.stimId(i)));
    if ~isfile(imgFile)
        % print(mp.plotTable.figureObj{1}, '-vector', imgFile, '-dpdf');
        exportgraphics(mp.plotTable.figureObj{1}, imgFile, ContentType="vector", BackgroundColor="none");
    end
end
