%% Plot activity heatmaps (interactive) and task phase responsiveness

anaDir = LMV.Data.GetAnalysisDir('pop_dynamics');

% ceName = "ce_m2_ex3_sentence-avg";
% figDir = LMV.Data.GetAnalysisDir('pop_dynamics', ceName+"_heatmaps");

%% Load data

% Load M2 PETHs
cachePath = fullfile(anaDir, ceName+".mat");
load(cachePath, 'ce');
clusTb = ce.clusTb;

% Load signrank responsiveness
sSR = LMV.Resp.LoadPhaseResponseTest();
sTb = NP.Unit.AlignClusTb(clusTb, sSR.clusTb, true);
sTb = [clusTb sTb];
ce.clusTb = sTb;

%% Plot heatmaps with signrank test results

f = MPlot.Figure(16);
f.WindowState = "maximized";

regions = ["mPrCG", "vPrCG", "IFG", "STG"];

for r = 1 : numel(regions)
    isRegion = strcmp(regions(r), ce.clusTb.region);
    if ~any(isRegion)
        isRegion = true(size(ce.clusTb.region));
    end
    isResp = ce.clusTb.tId1 ~= "none";
    ceReg = ce.Duplicate;
    ceReg.RemoveUnits(~isRegion | ~isResp);
    ceSen = ceReg.Split(ones(ceReg.numEpochs, 1));
    
    clf(f);
    tl = tiledlayout(2,7);
    tl.Padding = 'compact';
    
    for i = 1 : ceReg.numEpochs
        ax = nexttile;
        NP.UnitPlot.PethHeatmap(ceSen(i), [], 'SortVars', {'tId1', 'tId2', 'tId3', 'tR3'}, 'SortOrder', {'ascend', 'ascend', 'ascend', 'descend'}, ...
            'Ribbon', false, 'GroupColor', LMV.Param.GetTaskPhaseColors(sSR.phaseNames));
        ax.Title.String = sprintf("%s neurons (n = %i)", regions(r), ceReg.numResp);
        ax.YLabel.String = ceSen(i).GetTable("taskValue").stimText(1);
    end
    
    exportgraphics(f, fullfile(anaDir, sprintf("peth_heatmaps_%s.png", regions(r))));
end

%% Plot heatmaps by site

regions = ["mPrCG", "vPrCG", "IFG", "STG"];

for r = 1 : numel(regions)
    % Find recordings in this region
    isRegion = strcmp(regions(r), ce.clusTb.region);
    recIds = unique(ce.clusTb.recId(isRegion));
    
    % Make figure
    f = MPlot.Figure(36); clf
    % f.WindowState = 'maximized';
    tl = tiledlayout(2, 4);
    tl.Padding = 'compact';
    
    for i = 1 : numel(recIds)
        isRec = ce.clusTb.recId == recIds(i);
        isResp = ce.clusTb.tId1 ~= "none";
        ceReg = ce.Duplicate;
        ceReg.RemoveUnits(~isRec | ~isResp);
        
        ax = nexttile;
        NP.UnitPlot.PethHeatmap(ceReg, [], 'SortVars', {'tId1', 'tId2', 'tId3', 'tR3'}, 'SortOrder', {'ascend', 'ascend', 'ascend', 'descend'}, ...
            'GroupColor', LMV.Param.GetTaskPhaseColors(sSR.phaseNames));
        ax.Title.String = sprintf("%s neurons, %s (n = %i)", regions(r), recIds(i), ceReg.numResp);
        ax.Title.Interpreter = 'none';
        ax.YLabel.FontWeight = 'bold';
    end
    return
    
    figName = sprintf("peth_heatmaps_%s.png", regions(r));
    exportgraphics(f, fullfile(anaDir, "signrank", figName));
end

%% Compare different test results

regions = ["mPrCG", "vPrCG", "IFG", "STG"];

f = MPlot.Figure(225); clf
f.WindowState = 'maximized';
for r = 1 : numel(regions)
    NP.UnitPlot.PethRespTests(ce, 'signrank', sTb, 'ttest', tTb, 'Zeta', zTb, 'Region', regions(r));
    figPath = fullfile(anaDir, "benchmark", "peth_tests_heatmaps_" + regions(r) + ".png");
    MPlot.Paperize(f);
    exportgraphics(f, figPath);
end

return
%% Examine weights of specific units
%   See 2023/12/1 notes and random.pptx

% uTbs = NP.TRF.LoadModels('phase');
% uTb = cat(1, uTbs{:});

uIds = [430100570, 430100578, 430100210, 430100403, 430100027, 430100582, 430100551, 500300446];

f = MPlot.Figure(1); clf
for i = 1 : numel(uIds)
    ax = nexttile;
    uIdx = find(uTb.clusId==uIds(i), 1);
    NP.Resp.PlotWeights(gca, uTb, uIdx, 'phase');
end

f = MPlot.Figure(2); clf
LMV.Overview.SentencesFromCache(uIds, 'DataSource', 'm1', 'StimIdList', LMV.Param.stimIdList4([1 2 4]));
