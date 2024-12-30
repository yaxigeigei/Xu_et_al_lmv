%% 

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'rsa');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end

%% Load data

% Sentence averaged (M1) data
ceSearch = MBrowse.Dir2Table(fullfile(NP.Data.GetAnalysisRoot, 'data', 'ce_m1_sentence-avg', '*_ce_m1.mat'));
ceArray = NP.CE.LoadSession(fullfile(ceSearch.folder, ceSearch.name));

for i = 1 : numel(ceArray)
    ce = ceArray(i);
    
    % Remove sentences that have too few repeats
    nRep = ce.GetColumn('taskValue', 'numTrial');
    ce.RemoveEpochs(nRep < 3);
end

% Load unit responsiveness table
% sTest = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'zeta', "computed_zeta_clusTb.mat"));
sTest = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'ttest', "computed_ttest_clusTb.mat"));

%% Standalone RSA

mn = 'prod';
ce = ceArray(4);
disp(NP.SE.GetID(ce));

% Select responsive units
isRec = sTest.clusTb.recId == NP.SE.GetID(ce);
pval = sTest.clusTb{isRec,mn};
uMask = any(pval*size(pval,2) < 0.05, 2);
% disp(ce.clusTb.clusId(uMask));
% uMask = true(size(pval,1), 1); % select all units

% Configure analysis
ops = ce.userData.rsOps;

featFam = struct;
featFam.Full = [{'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', 'nasals', 'approximants', 'labial', 'velar', 'coronal', 'dental'}, ...
    ops.featVars.artic, ops.featVars.inten, ops.featVars.pitch];
% featFam.Phone = [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.consonants];
% featFam.MannerPlace = {'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', 'nasals', 'approximants', 'labial', 'velar', 'coronal', 'glottal', 'dental'};
% featFam.Artic = ops.featVars.artic;
% featFam.Acous = [ops.featVars.inten ops.featVars.pitch];
% switch mn
%     case 'stim'
%         featFam.Mel = {'speaker1'};
%     case 'prod'
%         featFam.Mel = {'mic'};
% end

ops.featFam = featFam;
ops.unitInd = find(uMask);
ops.maskName = mn;
switch mn
    case 'stim'
        ops.dtNeural = -0.15;
    case 'prod'
        ops.dtNeural = 0.15;
end
ops.downsample = 2;
ops.nBoot = 0; % not computing pval

ops.selection.nIter = 1e5;
ops.selection.uFrac = 0.25;
ops.selection.fFrac = 0.25;

sRez = NP.RSA.PrepareInput(ce, ops);

sRez = NP.RSA.RunSelection(sRez);

%% 

f = MPlot.Figure(4578); clf
tl = tiledlayout(1,2);

tl.Padding = 'compact';
NP.RSA.PlotCredits(sRez);

%% 

su = sRez.selection.simUnit;
sf = sRez.selection.simFeat;
unitId = sRez.ce.clusTb.clusId(isoutlier(su));
unitId = [410100245; unitId];

f = MPlot.Figure(5479); clf
LMV.Overview.SentencesFromCache(unitId, "m1", 'TaskPhase', sRez.maskName);











