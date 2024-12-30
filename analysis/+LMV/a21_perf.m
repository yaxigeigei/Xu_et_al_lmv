%% Population sentence classification

anaDir = LMV.Data.GetAnalysisDir("coding", "pop_sen4_ecoc", "rx_detrend");

%% Load data

% 
regions = ["all", "mPrCG"];
groups = ["bridge", "mirror"];

nRegions = numel(regions);
nGroups = numel(groups);
rxCell = cell(nGroups, nRegions);
for i = 1 : nGroups
    for j = 1 : nRegions
        rxCell{i,j} = load(fullfile(anaDir, sprintf("rx_%s_%s.mat", regions(j), groups(i))));
    end
end

% Load ce for plotting
load(fullfile(LMV.Data.GetAnalysisDir, "data", "ce_m2_ex3_sentence-avg.mat"), 'ce');

%%

rxCell = LMV.Decode.ComputeXTStats(rxCell);

%%

f = MPlot.Figure(53202); clf
tl = tiledlayout("flow");
tl.Padding = "compact";

for i = 1 : numel(rxCell)
    ax = nexttile;
    LMV.Decode.PlotAccuracyHeatmap(rxCell{i}, ce);
end

MPlot.Paperize(f, 1, .8);
exportgraphics(f, fullfile(anaDir, "rx_heatmap.png"));
exportgraphics(f, fullfile(anaDir, "rx_heatmap.pdf"), ContentType="vector");

%% 

f = MPlot.Figure(53201); clf
tl = tiledlayout("vertical");
tl.Padding = "compact";

for i = 1 : numel(rxCell)
    ax = nexttile;
    LMV.Decode.PlotAccuracyTimeseries(rxCell{i}, ce, "bestTrain");
end

MPlot.Paperize(f, .6, .8);
exportgraphics(f, fullfile(anaDir, "rx_best_train.png"));
exportgraphics(f, fullfile(anaDir, "rx_best_train.pdf"), ContentType="vector");

return
%% 

f = MPlot.Figure(53203); clf
tl = tiledlayout("vertical");
tl.Padding = "compact";

for i = 1 : numel(rxCell)
    ax = nexttile;
    LMV.Decode.PlotAccuracyTimeseries(rxCell{i}, ce, "diag");
end

MPlot.Paperize(f, .6, .8);
% exportgraphics(f, fullfile(anaDir, "r_diag.png"));