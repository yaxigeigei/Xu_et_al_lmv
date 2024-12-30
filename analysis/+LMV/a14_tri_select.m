%% Single-unit classification of sentences at each task phase

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'multiphase');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end

srcTb = NP.Data.FindSource('lmv');

%% Load data

% Task phase responsiveness
sResp = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'signrank', "computed_signrank_clusTb.mat"));
sResp.sigTb = NP.Resp.GetSigTable(sResp.clusTb);

% Sentence selectivity
selectDir = fullfile(NP.Data.GetAnalysisRoot, "sent_resp", "sen4");
sSelect = load(fullfile(selectDir, "computed_test_clusTb.mat"));
[~, I] = MMath.SortLike(sSelect.clusTb.clusId, sResp.clusTb.clusId);
sSelect.clusTb = sSelect.clusTb(I,:);

% Classifiers
mdlDir = fullfile(selectDir, "computed_mdls");
sCla.clusTb = NP.Fit.LoadModels([], mdlDir);
sCla.clusTb = cat(1, sCla.clusTb{:});
[~, I] = MMath.SortLike(sCla.clusTb.clusId, sResp.clusTb.clusId);
sCla.clusTb = sCla.clusTb(I,:);

%% 

uid = [410100463 410100441 410100450]; % 440300575

recId = "NP41_B1";
se = NP.SE.LoadSession(fullfile(NP.Data.GetAnalysisRoot, "data", "se_m1", recId+"_se_m1.mat"));
NP.Unit.SetUniqueClusId(se);
NP.Pitch.EnrichPitchTable(se);

%% 

senTb = NP.TaskBaseClass.SplitBySentence(se);
senTb = sortrows(senTb, 'numTrial', 'descend');

%% 

f = MPlot.Figure(123); clf
LMV.Plot.Overview12(senTb, 'UnitId', uid, 'TaskPhase', 'lmv');

% f = MPlot.Figure(456); clf
% LMV.Plot.Overview12(senTb, 'UnitId', uid, 'TaskPhase', 'lmv', 'Features', {'phone', 'acous', 'fujisaki'});

