%% Movies showing clusters of neuronal responses tiling various task phases

movDir = LMV.Data.GetAnalysisDir('movies', 'bridge_play');

%% Specify example unit

expId = 1;

switch expId
    case 1
        expClusId = 410100441;
        srFile = "sr_NP41_B1_2023-01-05_20-36-10_final.mat";
        datPath = "D:\preproc\sorting_workspace\NP41_B1_g0_imec0\temp_wh.dat";
        senInd = [11 5];
        trialInd = {3:6, 3:7};
        % tRanges = {[0 8.2], [0 7.5]};
end

%% Load data

% Load se
recId = NP.SE.GetID(expClusId);
srcTb = LMV.Data.FindSource(recId);
se = NP.SE.LoadSession(srcTb.path);

% Remove other units
isUnit = expClusId == NP.Unit.GetClusTb(se).clusId;
NP.Unit.RemoveUnits(se, ~isUnit);

% Preprocess se
se = LMV.SE.Transform(se, "enrich");

% Load sorting result object (for extracting spiking audio)
load(fullfile(LMV.Data.GetAnalysisDir, "misc", "data", srFile), "sr");
sr.mdat.Filename = datPath;

%% Make sentence table

% Generate performance metrics by finding morph times
NP.SE.SetMorphTimes(se);

% Split epochs by sentences
senTb = LMV.SE.SplitBySentence(se); % poor trials are removed before spliting

% Sort trials by the duration of delay
for i = 1 : height(senTb)
    tt = senTb.se(i).GetTable("taskTime");
    delayDur = tt.prodMatchOn - tt.stimMatchOff;
    [~, I] = sort(delayDur, "ascend");
    senTb.se(i).SortEpochs(I);
end

%% Select sentence and trials

k = 1;

seSen = senTb.se(senInd(k)).Duplicate;
seSen.userData.expClusTb = NP.Unit.GetClusTb(seSen);

isRm = true(seSen.numEpochs, 1);
isRm(trialInd{k}) = false;
seSen.RemoveEpochs(isRm);

%% 

import LMV.Linker.Mov

if ~exist("mp", 'var')
    mp = MPlotter();
end

mp.RemovePlot();

figNum = 1;
rowDist = [1 2 3];
colDist = 1;

spInd = MPlot.FindSubplotInd(rowDist, colDist, 1, 1);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @cla, '', 'epoch');
% mp.AddPlot(figNum, spStr, @(a,b) Mov.TrialPhone(a, b), 'seSen', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) Mov.TrialEvents(a, b), 'seSen', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) Mov.SentenceTitle(a, b), 'seSen', 'epoch');
mp.AddPlot(figNum, spStr);

spInd = MPlot.FindSubplotInd(rowDist, colDist, 2, 1);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @cla, '', 'epoch');
mp.AddPlot(figNum, spStr, @Mov.TrialMelSpectrogram, 'seSen', 'epoch');
mp.AddPlot(figNum, spStr);

spInd = MPlot.FindSubplotInd(rowDist, colDist, 3, 1);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @cla, '', 'epoch');
mp.AddPlot(figNum, spStr, @Mov.DelayBoundary, 'seSen', 'epoch');
mp.AddPlot(figNum, spStr, @Mov.RealtimeRaster, 'seSen', 'epoch');
mp.AddPlot(figNum, spStr, @Mov.UpdateRaster, 'seSen', 'time');
mp.AddPlot(figNum, spStr);

mp.RefreshAll();
f = figure(mp.plotTable.figureObj{1});
% f.Position(1:2) = [0 50];
f.Position(3:4) = [1280 360];
MPlot.Paperize(f, 'FontSize', 8);

mp.timeLimits = [0 seSen.GetTable("taskTime").prodMatchOff(end)+0.3];

%% 

LMV.Linker.Mov.MakeBridgePlayVideos(movDir, mp, seSen);

%% Make speech audio

LMV.Linker.Mov.MakeBridgePlaySpeechAudio(movDir, seSen);

%% Make spike audio

LMV.Linker.Mov.MakeBridgePlaySpikeAudio(movDir, seSen, se, sr);
