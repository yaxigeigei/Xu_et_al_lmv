%% Plot example feature timeseries and neuronal activity used in fitting encoding models

anaDir = LMV.Data.GetAnalysisDir('coding', 'stRF');
recId = 'NP41_B1';
srcTb = LMV.Data.FindSource(recId);

%% Load data

se = NP.SE.LoadSession(srcTb.path);
load(fullfile(anaDir, "ce", recId+"_ce.mat"));

%% Plot example timeseries

[tt, tv] = se.GetTable('taskTime', 'taskValue');
rt = se.GetReferenceTime('taskTime');

isStim = startsWith(tv.stimText, "we've got");
indStim = find(isStim);
k = indStim(2);
% tWinEp = [tt.prodOn{k} tt.prodOff{k}];
tWinEp = tt.prodOn{k} + [0 1.01];
tWin = rt(k) + tWinEp;

uIdx = NP.Unit.ClusId2Ind(410100245, se);

f = MPlot.Figure(1077); clf
tl = tiledlayout("vertical");
tl.Padding = "compact";

ax = nexttile;
NP.TaskBaseClass.PlotTGEHier(ax, se, k, tWinEp, 'prod');
ax.XTick = [];
ax.XLabel.String = [];
ax.Visible = 'off';

ax = nexttile;
NP.UnitPlot.RasterStack(ax, se, k, tWinEp, uIdx, 'PETH', true);
ax.XTick = [];
ax.XLabel.String = [];

ax = nexttile;
LMV.TRF.PlotInputs(ax, ce, "mic", tWin);
% ax.Title.String = tv.prodText{k};

ax = nexttile;
LMV.TRF.PlotInputs(ax, ce, LMV.TRF.GetFeatureSet('artic3'), tWin, 'Normalization', 'zscore');

ax = nexttile;
LMV.TRF.PlotInputs(ax, ce, LMV.TRF.GetFeatureSet('phone'), tWin);

ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';
ax.XLabel.String = "Time (s)";
MPlot.Paperize(f, .55, 1); %return
exportgraphics(f, fullfile(anaDir, "example_feature_timeseries.png"));
exportgraphics(f, fullfile(anaDir, "example_feature_timeseries.pdf"), ContentType="vector");
