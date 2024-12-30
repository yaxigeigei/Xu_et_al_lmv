%% Movies showing clusters of neuronal responses tiling various task phases

anaDir = LMV.Data.GetAnalysisDir('movies', 'signrank_m1_resp');

%% Load the seed session

% Load and preprocess se
recId = "NP41_B1";
srcTb = LMV.Data.FindSource(recId);
se = NP.SE.LoadSession(srcTb.path);
ops = NP.Param.Enrich;
ops.isSpkRate = true;
ops.isMel = true;
ops.isPitch = true;
ops.isArtic = true;
NP.SE.Enrich(se, ops);

% Remove non-responsive units
rTestSeed = LMV.Resp.LoadPhaseResponseTest(recId);
NP.Unit.RemoveUnits(se, ~any(rTestSeed.sigTb{:,:}, 2));

% Load sorting result object (for extracting spiking audio)
srPath = "C:\chang_lab\project_np\preproc\sorting_workspace\NP41_B1_g0_imec0\mtracer_cache\sr_NP41_B1_2023-01-05_20-36-10_final.mat";
if exist(srPath, "file")
    load(srPath, "sr");
end

%% Load cached results

% M1 sentence PETHs
% This will be reverse morphed to example trials in se
load(fullfile(LMV.Data.GetAnalysisDir, "data", "ce_m1_sentence-avg.mat"), 'ce');

% Remove non-responsive units
rTest = LMV.Resp.LoadPhaseResponseTest();
[~, I] = MMath.SortLike(rTest.clusTb.clusId, ce.clusTb.clusId);
ce.clusTb.tId1 = rTest.clusTb.tId1(I);
ce.RemoveUnits(ce.clusTb.tId1=="none");

% UMAP embedding and M2 session PETHs
% These PETHs are only used to find normalization parameters
sEmbed = load(fullfile(LMV.Data.GetAnalysisDir, 'embedding', 'umap', 'computed_umap_all.mat'));
[~, I] = MMath.SortLike(sEmbed.ce.clusTb.clusId, ce.clusTb.clusId);
ce.clusTb.embedCoords = -sEmbed.ce.clusTb.embedCoords(I,:); % flip axes for better orientation

%% Find template trials from the example recording

% Find morph times to the best matched trials
NP.SE.SetMorphTimes(se, [], 'best');

% Morph se and make sentence table
seMorph = NP.SE.MorphSession(se);
NP.Unit.AddSimSpikeTimeTable(seMorph);
[ceSeed, senTb] = LMV.SE.ComputeSentencePETH(seMorph);

% Sort rows of senTb to be consistent with ce epochs
[~, I] = MMath.SortLike(senTb.stimId, ce.GetTable('taskValue').stimId);
senTb = senTb(I,:);

% Put template trials to a new se (sentences should be in the same order as those in ce)
clear seTemp
for i = height(senTb) : -1 : 1
    m = senTb.trialNum{i} == senTb.tempTrialNum(i);
    seTemp(i) = senTb.se(i).Split({m});
end
seTemp = Merge(seTemp);
seTemp.userData = seTemp.userData(1);

%% Morph responses to template trials

% Concatenate task tables from template se and ce
[tt1, tv1] = seTemp.GetTable('taskTime', 'taskValue');
[tt2, tv2] = ce.GetTable('taskTime', 'taskValue');
tv2.alignScore(:) = -Inf; % a trick to make these trials not be templates

vn = intersect(tt1.Properties.VariableNames, tt2.Properties.VariableNames);
tt = [tt1(:,vn); tt2(:,vn)];

vn = intersect(tv1.Properties.VariableNames, tv2.Properties.VariableNames);
tv = [tv1(:,vn); tv2(:,vn)];

% Morph ce onto template trials
[tt, tv] = LMV.SE.FindMorphTimes(tt, tv, [], 'best');
ce.SetTable('taskTime', tt(height(tt1)+1:end, :));
ceMorph = NP.SE.MorphEpochs(ce);

% Add task phase events
phases = struct;
phases.atten = {'cue1On', 'stimOn'};
phases.delay = {'stimOff', 'cue3On'};
phases.init = {'cue3On', 'prodOn'};
NP.TaskBaseClass.AddEventObjects(ceMorph, phases, 'taskTime');

%% Adjust embedding responses

% Get baseline responses and normalization parameters from M2 session PETHs
[T, R] = sEmbed.ce.GetArray('resp');
rBase = mean(R(T>-LMV.Param.respBaselineDur & T<0, :), 1, 'omitmissing');
R = R - rBase;
[~, c, k] = MMath.Normalize(R, 'maxsoft', 1);

% Baseline subtraction, normalization, smoothing
resp = ceMorph.GetTable('resp');
fs = 1 / ceMorph.userData.rsOps.rsBinSize;
for i = 1 : height(resp)
    for j = 2 : width(resp)
        r = resp.(j){i};
        r = r - rBase(j-1);
        r = (r - c(j-1)) ./ k(j-1);
        r = MNeuro.Filter1(r, fs, 'gaussian', 0.1); % 0.1 sec SD Gaussian kernel
        resp.(j){i} = r;
    end
end
ceMorph.SetTable('spikeRate', resp, 'timeSeries');

%% Add embedding coordinates and responsiveness to seedClusTb for marking example units

allClusTb = ceMorph.clusTb;
seedClusTb = ceSeed.clusTb;
[~, I] = MMath.SortLike(allClusTb.clusId, seedClusTb.clusId, false);
varNames = {'tId1', 'embedCoords'};
for i = 1 : numel(varNames)
    vn = varNames{i};
    seedClusTb.(vn) = allClusTb.(vn)(I,:);
end

%% 

% Add selected example units to each se in senTb
cid = 4101e5 + [420 454 138 272 171 463 245 224]; % alt examples: 224 48 249(too sparse) 405(like 420)
uInd = NP.Unit.ClusId2Ind(cid, seedClusTb);
for i = 1 : height(senTb)
    senTb.se(i).userData.expClusTb = seedClusTb(uInd,:);
end

% Choose unit to make spike audios
cidSpk = repmat([454 245], height(senTb), 1);

%% Make movies with example unit responses

if ~exist('mp', 'var')
    mp = MPlotter();
end
mp.RemovePlot();

figNum = 1;
rowDist = [2 9];
colDist = [2 3];

spInd = MPlot.FindSubplotInd(rowDist, colDist, 2:numel(rowDist), 1);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @LMV.MP.EmbeddedResp, 'ceMorph', 'time');
mp.AddPlot(figNum, spStr, @LMV.MP.LabelEmbeddedUnits, 'senTb', 'epoch');

spInd = MPlot.FindSubplotInd(rowDist, colDist, 1, 2);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @cla, '', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.TrialPhone(a, b), 'ceMorph', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.TrialEvents(a, b), 'ceMorph', 'epoch');
mp.AddPlot(figNum, spStr, @(a,b) LMV.MP.SentenceTitle(a, b), 'ceMorph', 'epoch');
mp.AddPlot(figNum, spStr);

spInd = MPlot.FindSubplotInd(rowDist, colDist, 2, 2);
spStr = MPlotter.SubplotArgs2Str(spInd{:});
mp.AddPlot(figNum, spStr, @LMV.MP.RasterStack, 'senTb', 'epoch');
mp.AddPlot(figNum, spStr);

mp.RefreshAll();
f = figure(1);
f.Position(1:2) = [0 50];
f.Position(3:4) = [1920 800];
MPlot.Paperize(f, 'FontSize', 8);

%% 

recId = NP.SE.GetID(seTemp);
[tt, tv] = seTemp.GetTable('taskTime', 'taskValue');

rt = se.GetReferenceTime;

% Find the indices of well repeated sentences
[~, I] = sort(senTb.numTrial, 'descend');
% senInd = I(1:4);
senInd = find(senTb.numTrial >= 4);

%% 

for i = senInd(1:end)'
    % Get time window
    tWin = [tt.trialOn(i) tt.matchOff(i)] + [-1 1]*0.3;
    
    % Set epoch and time
    mp.timeLimits = tWin;
    mp.epoch = i;
    mp.time = tWin(1);
    
    MPlot.Paperize(f, 'FontSize', 8);
    
    % Make movie
    vidFile = fullfile(anaDir, sprintf("sent%i_%s_trial%i.mp4", i, recId, tv.trialNum(i)));
    fps = 60;
    frames = mp.MakeVideo(f, 1/fps, 'FrameRate', fps, 'FilePath', vidFile);
    
    % Save mic audio
    micFile = sprintf("sent%i_%s_trial%i.wav", i, recId, tv.trialNum(i));
    NP.Audio.WriteTrial(seTemp, i, tWin, 'FolderPath', anaDir, 'FileName', micFile);
end

%% Save spike audio

if ~exist("sr", "var")
    return
end

for i = senInd(1:end)'
    % Get time window
    tWin = [tt.trialOn(i) tt.matchOff(i)] + [-1 1]*0.3;
    
    for j = size(cidSpk,2)
        spkFile = fullfile(anaDir, sprintf("sent%i_%s_trial%i_u%i.wav", i, recId, tv.trialNum(i), cidSpk(i,j)));
        isTrial = se.GetTable('taskValue').trialNum == tv.trialNum(i);
        sr.WriteSpikeAudio(spkFile, cidSpk(i,j), tWin+rt(isTrial));
    end
end
