%% 

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'rsa', 'from_matt');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end

%% Load data

% load(fullfile(anaDir, "NP04_B2_preproc.mat"), 'out');
% load(fullfile(anaDir, "NP04_full_trf.mat"), 'trf');
load(fullfile(anaDir, "NP04_uniqueR2.mat"), 'mdl_cmp');
load(fullfile(anaDir, "NP04_B2-names_full_model-rep0.mat"), 'famNames', 'featNames');
load(fullfile(anaDir, "sent_ord.mat"), 'sent_ord');
featArr = readNPY(fullfile(anaDir, "NP04_B2-feature_tc_full_model-rep0.npy"));
respArr = readNPY(fullfile(anaDir, "NP04_B2-psth_tc_full_model-rep0.npy"));

famNames = string(famNames);
featNames = string(featNames);

%% Convert single timepoint feature to interval

% featArr = readNPY(fullfile(anaDir, "NP04_B2-feature_tc_full_model-rep0.npy"));

phRange = 21:36;
phArr = featArr(phRange,:);
indPh = find(sum(phArr));
for i = 1 : size(phArr,1)
    indPhEach = find(phArr(i,:));
    for j = 1 : numel(indPhEach)
        k = indPhEach(j);
        itvl = min(indPh(indPh > k) - k);
        if isempty(itvl) || itvl > 15
            phArr(i,k:k+15-1) = phArr(i,k);
        else
            phArr(i,k:k+itvl-1) = phArr(i,k);
        end
    end
end     
featArr(phRange,:) = phArr;

%% Select units

respArr = respArr(mdl_cmp.sig_units, :);

%% 

t = (1:size(featArr,2))'*1e-2;
featCell = mat2cell(featArr', numel(t), ones(size(featNames)));
respCell = mat2cell(respArr', numel(t), ones(size(mdl_cmp.sig_units)));

fTb = MSessionExplorer.MakeTimeSeriesTable(t, featCell, 'VariableNames', featNames);
rTb = MSessionExplorer.MakeTimeSeriesTable(t, respCell);

ce = NP.CodingExplorer;
ce.userData.expInfo.subjectId = 'NP04';
ce.userData.expInfo.blockId = 'B2';
ce.SetTable('feat', fTb, 'timeSeries', 0);
ce.SetTable('resp', rTb, 'timeSeries', 0);

tSlice = t(logical(fTb.sent_onset{1}));
ce.SliceSession(tSlice, 'absolute');

% 
stimId = string(sent_ord)';
[~, ia, ic] = unique(stimId);
stimDur = [diff(tSlice)-1.5; 10];
stimDur = round(stimDur(ia), 2);
stimDur = stimDur(ic);

aTb = table;
aTb.stimId = stimId;
aTb.duration = stimDur;
ce.SetTable('attr', aTb, 'eventValues');

% Slice sentences
fTb = ce.SliceTimeSeries('feat', [zeros(ce.numEpochs,1) aTb.duration] + [-1 1]*1e-3);
rTb = ce.SliceTimeSeries('resp', [zeros(ce.numEpochs,1) aTb.duration] + [-1 1]*1e-3 + 0.12); % compensate time lag in responses
ce.SetTable('feat', fTb);
ce.SetTable('resp', rTb);

% 
ceTb = ce.SplitConditions('stimId', 'attr');
ceTb(ceTb.numEpochs==1, :) = [];

%% 

for i = 1 : height(ceTb)
    [ceTb.time{i}, ceTb.mFeat{i}] = ceTb.se(i).GetArray('feat', 'DimCat', 3, 'DimAverage', 3);
    [~, ceTb.mResp{i}] = ceTb.se(i).GetArray('resp', 'DimCat', 3, 'DimAverage', 3);
    ceTb.stimIdx{i} = repmat(i, size(ceTb.time{i}));
end

stimIdx = cat(1, ceTb.stimIdx{:});
t = cat(1, ceTb.time{:});
mF = cat(1, ceTb.mFeat{:});
mR = cat(1, ceTb.mResp{:});

%% Inspect input

% Raw input
% featArr = readNPY(fullfile(anaDir, "NP04_B2-feature_tc_full_model-rep0.npy"));
% iSentOn = find(featArr(1,:));
% a = iSentOn(101);
% F = featArr(:,a:a+250);
% F = MMath.Normalize(F', 'minmax')';

% Sentence average
F = mF';
F = ceTb.mFeat{1}';
F = MMath.Normalize(F', 'minmax')';

% Plot
MPlot.Figure(345); clf
imagesc(F);
% colormap gray
ax = gca;
ax.YTick = 1 : numel(featNames);
ax.YTickLabel = featNames;
ax.TickLabelInterpreter = 'none';

%% 

famList = unique(famNames, 'stable');
iSelect = [2 3 4];
famSelect = famList(iSelect);

ss = cell(size(iSelect));

for i = 1 : numel(iSelect)
    % Get feats belonging to the feature family of interest
    k = iSelect(i);
    isFeat = famNames == famList(k);
    F = mF(:,isFeat);
    
    % Get units that have sig unique-r for this feat family
    ur = mdl_cmp.uniqueR2(:,k);
    isResp = ur > median(abs(ur))*2;
    R = mR(:,isResp);
    
    % Compute RDMs and similarity
    s = struct;
    s.t = t;
    s.stimIdx = stimIdx;
    s.repNames = ["Neural", famList(k)];
    s.input = {R, F};
    s.nBoot = 100;
    
    s = NP.RSA.ComputeRDMs(s);
    s = NP.RSA.ComputeSimilarity(s);
    
    ss{i} = s;
end

ss = cat(1, ss{:});

%% Plot RDMs

f = MPlot.Figure(160); clf
tl = tiledlayout(numel(ss), 2);
tl.Padding = 'compact';
for i = 1 : numel(ss)
    s = ss(i);
    
    ax = nexttile;
    NP.RSA.PlotRDM(s.RDM{1}, s.stimIdx);
    ax.Title.String = s.repNames(1);
    
    ax = nexttile;
    NP.RSA.PlotRDM(s.RDM{2}, s.stimIdx);
    ax.Title.String = s.repNames(2) + sprintf(", r = %.2f (p = %.2f)", s.PSM(2), s.pval(2));
end

