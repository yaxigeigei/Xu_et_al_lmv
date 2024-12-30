%% NMF clustering on grand average PETHs

anaDir = LMV.Data.GetAnalysisDir('phase_resp', 'nmf');

%% Load computed nmfc

regions = ["mPrCG", "vPrCG", "IFG", "STG"];
sBatch = cell(size(regions));
for r = 1 : numel(regions)
    regionCachePath = fullfile(anaDir, "computed_nmfc_"+regions(r)+".mat");
    sBatch{r} = load(regionCachePath);
end
sBatch = cat(1, sBatch{:});

%% Find optimal number of clusters

for r = 1 : numel(sBatch)
    % Find optimal number of clusters
    sEval = sBatch(r).sEval;
    sEval.optiK = median(sEval.bootOptiK);
    sBatch(r).sEval = sEval;
    
    % Add cluster info to unit table
    ce = sBatch(r).ce;
    sOpti = sBatch(r).sClus(sEval.kList==sEval.optiK);
    ce.clusTb.nmfcId = sOpti.clustId;
    ce.clusTb.nmfcW = max(sOpti.W, [], 1)';
end

%% Gap statistics

f = MPlot.Figure(46533); clf

for r = 1 : numel(regions)
    s = sBatch(r).sEval;
    % [m, sd] = MMath.MeanStats(s.bootScore, 2);
    
    ax = nexttile;
    histogram(s.bootOptiK)
    ax.XTick = s.kList;
    xlabel('Optimal # of clusters')
    ylabel('Count of iterations')
    title(sprintf("%s, median = %i", sBatch(r).region, s.optiK));
    MPlot.Axes(ax);
end

MPlot.Paperize(f, 1, .6);
exportgraphics(f, fullfile(anaDir, "nmf_boot_k.png"));

%% Bases at different cluster numbers

k2Plot = 1 : 10;

for r = 1 : numel(regions)
    f = MPlot.Figure(46534); clf
    LMV.NMF.PlotBases(sBatch(r), k2Plot);
    MPlot.Paperize(f, 2.2, .8);
    exportgraphics(f, fullfile(anaDir, "nmf_bases_"+regions(r)+".png"));
end

%% Plot cluster means

for r = 1 : numel(regions)
    f = MPlot.Figure(225); clf
    LMV.NMF.PlotClusterMean(sBatch(r).ce);
    title(sBatch(r).region);
    MPlot.Paperize(f, 'ColumnsWide', 0.5, 'AspectRatio', 2.5);
    exportgraphics(f, fullfile(anaDir, "nmf_cluster_means_"+regions(r)+".png"));
end

return
%% Plot heatmaps

f = MPlot.Figure(15); clf
tl = tiledlayout(1, numel(regions));
tl.Padding = 'compact';

for r = 1 : numel(regions)
    s = sBatch(r);
    ce = s.ce;
    sEval = s.sEval;
    
    ceSub = ce.Duplicate;
    if s.region == "mPrCG"
        cid = ceSub.clusTb.nmfcId;
        cid(cid==1) = 6;
        cid = cid - 1;
        ceSub.clusTb.nmfcId = cid;
    end
    % isRegion = strcmp(regions(r), ce.clusTb.region);
    % if ~any(isRegion)
    %     isRegion = true(size(ce.clusTb.region));
    % end
    % ceSub = ce.Duplicate;
    % ceSub.RemoveUnits(~isRegion);
    
    ax = nexttile;
    NP.UnitPlot.PethHeatmap(ceSub, 1:sEval.optiK, 'SortVars', {'nmfcId', 'nmfcW'}, 'SortOrder', {'ascend', 'descend'});
%     title(sprintf("%s (n = %i)", regions(r), ceSub.numResp), 'Interpreter', 'none');
    ax.YLabel.String = sprintf("%s neurons (n = %i)", regions(r), ceSub.numResp);
    ax.YLabel.FontWeight = 'bold';
end

MPlot.Paperize(f, 'ColumnsWide', 2.8, 'ColumnsHigh', 1);
exportgraphics(f, fullfile(anaDir, "peth_heatmaps.png"));


