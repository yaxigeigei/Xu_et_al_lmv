%% Plot activity heatmaps (interactive) and task phase responsiveness

anaDir = LMV.Data.GetAnalysisDir('phase_resp');

%% Load data

% Load M2 PETHs
cachePath = fullfile(LMV.Data.GetAnalysisDir, "data", "ce_m2_session-avg.mat");
load(cachePath, 'ce');
clusTb = ce.clusTb;
clusTb.probeId = NP.Unit.GetProbeId(clusTb.clusId);

% Add NMF cluster ID
regions = ["mPrCG", "vPrCG", "IFG", "STG"];
for r = 1 : numel(regions)
    sNMF = load(fullfile(anaDir, 'nmf', "computed_nmfc_"+regions(r)+".mat"));
    if sNMF.region == "mPrCG"
        cid = sNMF.ce.clusTb.nmfcId;
        cid(cid==1) = 6;
        cid = cid - 1;
        sNMF.ce.clusTb.nmfcId = cid;
    end
    [~, I] = MMath.SortLike(clusTb.clusId, sNMF.ce.clusTb.clusId);
    clusTb.nmfcId(I) = sNMF.ce.clusTb.nmfcId;
    clusTb.nmfcW(I) = sNMF.ce.clusTb.nmfcW;
end
ce.clusTb = clusTb;

% Load signrank responsiveness
sSR = LMV.Resp.LoadPhaseResponseTest();
sTb = NP.Unit.AlignClusTb(clusTb, sSR.clusTb, true);

%% Plot heatmaps of responsive units ordered by NMF clusters

regions = ["mPrCG", "vPrCG", "IFG", "STG"];

f = MPlot.Figure(16); clf
tl = tiledlayout(1, numel(regions));
tl.Padding = 'compact';

for r = 1 : numel(regions)
    isRegion = strcmp(regions(r), ce.clusTb.region);
    if ~any(isRegion)
        isRegion = true(size(ce.clusTb.region));
    end
    isResp = sTb.tId1 ~= "none";
    ceSub = ce.Duplicate;
    ceSub.RemoveUnits(~isRegion | ~isResp);
    
    ax = nexttile;
    NP.UnitPlot.PethHeatmap(ceSub, [], 'SortVars', {'nmfcId', 'nmfcW'}, 'SortOrder', {'ascend', 'descend'});
    ax.Title.String = sprintf("%s neurons (n = %i)", regions(r), ceSub.numResp);
    ax.YLabel.FontWeight = 'bold';
end

MPlot.Paperize(f, 'ColumnsWide', 2.2, 'ColumnsHigh', .7);
% exportgraphics(f, fullfile(anaDir, "nmf", "peth_heatmaps.png"));
% exportgraphics(f, fullfile(anaDir, "nmf", "peth_heatmaps.pdf"));

%% Plot heatmaps by site

regions = ["mPrCG", "vPrCG", "IFG", "STG"];

for r = 1 : numel(regions)
    % Select responsive units
    isResp = sTb.tId1 ~= "none";
    
    % Find recordings in this region
    region = regions(r);
    isRegion = strcmp(region, ce.clusTb.region);
    recIds = unique(ce.clusTb.recId(isRegion));
    
    % Make figure
    f = MPlot.Figure(36); clf
    tl = tiledlayout(3, 4);
    tl.Padding = 'compact';
    
    for i = 1 : numel(recIds)
        % Find probes used in this recording
        isRec = ce.clusTb.recId == recIds(i);
        prbIds = unique(ce.clusTb.probeId(isRec));
        
        for p = 1 : numel(prbIds)
            isPrb = ce.clusTb.probeId == prbIds(p);
            ceSub = ce.Duplicate;
            ceSub.RemoveUnits(~isRec | ~isPrb |  ~isResp);
            
            ax = nexttile;
            NP.UnitPlot.PethHeatmap(ceSub, [], 'SortVars', {'nmfcId', 'nmfcW'}, 'SortOrder', {'ascend', 'descend'});
            ax.XLabel.String = [];
            ax.Title.String = sprintf("%s neurons, %s #%i (n = %i)", region, recIds(i), prbIds(p), ceSub.numResp);
            ax.Title.Interpreter = 'none';
            ax.YLabel.FontWeight = 'bold';
        end
    end
    
    MPlot.Paperize(f, 'ColumnsWide', 2.2, 'ColumnsHigh', 1.5);
    exportgraphics(f, fullfile(anaDir, "nmf", sprintf("peth_heatmaps_%s.png", region)));
    exportgraphics(f, fullfile(anaDir, "nmf", sprintf("peth_heatmaps_%s.pdf", region)));
end

%% Load additional test results

% Load t-test responsiveness
sTT = load(fullfile(anaDir, "ttest", "computed_test_clusTb.mat"));
tTb = NP.Unit.AlignClusTb(clusTb, sTT.clusTb, false);

% Load Zeta responsiveness
sZeta = load(fullfile(anaDir, 'zeta', 'computed_test_clusTb.mat'));
zTb = NP.Unit.AlignClusTb(clusTb, sZeta.clusTb, false);

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
%% Plot heatmaps with signrank test results

regions = ["mPrCG", "vPrCG", "IFG", "STG"];

f = MPlot.Figure(16); clf
tl = tiledlayout(1, numel(regions));
tl.Padding = 'compact';

for r = 1 : numel(regions)
    isRegion = strcmp(regions(r), ce.clusTb.region);
    if ~any(isRegion)
        isRegion = true(size(ce.clusTb.region));
    end
    isResp = ce.clusTb.tId1 ~= "none";
    ceSub = ce.Duplicate;
    ceSub.RemoveUnits(~isRegion | ~isResp);
    
    ax = nexttile;
    NP.UnitPlot.PethHeatmap(ceSub, [], 'SortVars', {'tId1', 'tId2', 'tId3', 'tR3'}, 'SortOrder', {'ascend', 'ascend', 'ascend', 'descend'}, ...
        'GroupColor', LMV.Param.GetTaskPhaseColors(sSR.phaseNames));
    ax.Title.String = sprintf("%s neurons (n = %i)", regions(r), ceSub.numResp);
    ax.YLabel.FontWeight = 'bold';
end

MPlot.Paperize(f, 'ColumnsWide', 2.2, 'ColumnsHigh', .7);
exportgraphics(f, fullfile(anaDir, "signrank", "peth_heatmaps.png"));
exportgraphics(f, fullfile(anaDir, "signrank", "peth_heatmaps.pdf"));

%% Plot heatmaps by site

regions = ["mPrCG", "vPrCG", "IFG", "STG"];

for r = 1 : numel(regions)
    % Select responsive units
    isResp = ce.clusTb.tId1 ~= "none";
    
    % Find recordings in this region
    region = regions(r);
    isRegion = strcmp(region, ce.clusTb.region);
    recIds = unique(ce.clusTb.recId(isRegion));
    
    % Make figure
    f = MPlot.Figure(36); clf
    tl = tiledlayout(3, 4);
    tl.Padding = 'compact';
    
    for i = 1 : numel(recIds)
        % Find probes used in this recording
        isRec = ce.clusTb.recId == recIds(i);
        prbIds = unique(ce.clusTb.probeId(isRec));
        
        for p = 1 : numel(prbIds)
            isPrb = ce.clusTb.probeId == prbIds(p);
            ceSub = ce.Duplicate;
            ceSub.RemoveUnits(~isRec | ~isPrb |  ~isResp);
            
            ax = nexttile;
            NP.UnitPlot.PethHeatmap(ceSub, [], 'SortVars', {'tId1', 'tId2', 'tId3', 'tR3'}, 'SortOrder', {'ascend', 'ascend', 'ascend', 'descend'}, ...
                'GroupColor', LMV.Param.GetTaskPhaseColors(sSR.phaseNames));
            % ax.XLabel.String = [];
            ax.Title.String = sprintf("%s neurons, %s #%i (n = %i)", region, recIds(i), prbIds(p), ceSub.numResp);
            ax.Title.Interpreter = 'none';
            ax.YLabel.FontWeight = 'bold';
        end
    end
    
    MPlot.Paperize(f, 'ColumnsWide', 2.2, 'ColumnsHigh', 1.5); %return
    exportgraphics(f, fullfile(anaDir, "signrank", sprintf("peth_heatmaps_%s.png", region)));
    exportgraphics(f, fullfile(anaDir, "signrank", sprintf("peth_heatmaps_%s.pdf", region)));
end

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
