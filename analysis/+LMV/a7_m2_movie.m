%% Movies showing clusters of neuronal responses tiling various task phases

anaDir = LMV.Data.GetAnalysisDir('movies', 'signrank_m2_resp');

%% Load data

% Load PETHs and UMAP embeddings
cachePath = fullfile(LMV.Data.GetAnalysisDir, 'embedding', 'umap', 'computed_umap_all.mat');
load(cachePath);
ce.clusTb.embedCoords = -ce.clusTb.embedCoords; % flip axes for better orientation

% Load responsiveness
rTest = LMV.Resp.LoadPhaseResponseTest();
sTb = NP.Unit.AlignClusTb(ce.clusTb, rTest.clusTb, true);
ce.clusTb = [ce.clusTb sTb];

% Keep responsive mPrCG units
isRm = ce.clusTb.region ~= "mPrCG" | ce.clusTb.tId1 == "none";
ceSub = ce.Duplicate;
ceSub.RemoveUnits(isRm);

%% Add data to ce

% Baseline subtraction
respTb = ceSub.GetTable('resp');
for i = 1 : height(respTb)
    isBaseline = respTb.time{i} < 0 & respTb.time{i} > -LMV.Param.respBaselineDur;
    for j = 2 : width(respTb)
        r = respTb.(j){i};
        r = r - mean(r(isBaseline));
        respTb.(j){i} = r;
    end
end
ceSub.SetTable('resp', respTb);

% Smoothing (save to spikeRate table)
fs = 1 / ceSub.userData.rsOps.rsBinSize;
respTb = ceSub.GetTable('resp');
for i = 1 : height(respTb)
    for j = 2 : width(respTb)
        r = respTb.(j){i};
        r = MNeuro.Filter1(r, fs, 'gaussian', 0.1); % 0.1 sec SD Gaussian kernel
        r = MMath.Normalize(r, 'max')*1.3;
        respTb.(j){i} = r;
    end
end
ceSub.SetTable('spikeRate', respTb, 'timeSeries');

% Add task phase events
phases = struct;
phases.atten = {'cue1On', 'stimOn'};
phases.delay = {'stimOff', 'cue3On'};
phases.init = {'cue3On', 'prodOn'};
NP.TaskBaseClass.AddEventObjects(ceSub, phases, 'taskTime');

%% Make movies with example unit responses

if ~exist('mp', 'var')
    mp = MPlotter();
end
mp.RemovePlot();

figNum = 1;
rowDist = [1 9];
colDist = [2 3];

spInd = MPlot.FindSubplotInd(rowDist, colDist, 2:numel(rowDist), 1);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @LMV.MP.EmbeddedResp, 'ceSub', 'time');

spInd = MPlot.FindSubplotInd(rowDist, colDist, 1, 2);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @cla, '', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.TrialEvents(a, b), 'ceSub', 'epoch');
mp.AddPlot(figNum, spStr);

spInd = MPlot.FindSubplotInd(rowDist, colDist, 2, 2);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @LMV.MP.PethHeatmap, 'ceSub', 'epoch');
mp.AddPlot(figNum, spStr);

mp.RefreshAll();
f = figure(figNum);
f.Position(1:2) = [50 400];
f.Position(3:4) = [1920 720];

%% 

% Set time window
tt = ce.GetTable('taskTime');
tWin = [tt.trialOn(i) tt.matchOff(i)] + [-1 1]*0.5;
mp.timeLimits = tWin;

% Set epoch and time
mp.epoch = 1;
mp.time = tWin(1);

MPlot.Paperize(f);

% Make movie
fps = 60;
vidFile = fullfile(anaDir, sprintf("m2_resp_signrank_%ifps.mp4", fps));
frames = mp.MakeVideo(f, 1/fps, 'FilePath', vidFile, 'FrameRate', fps);
