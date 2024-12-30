%% PETH heatmap of units from a single mPrCG site

figDir = "C:\chang_lab\project_np\misc\20240425_simon_keynote";

%% Load cached results

% NMF clustering results and the grand average PETHs
% The grand average is only used to find normalization parameters
cachePath = fullfile(LMV.Data.GetAnalysisDir, 'data', "ce_m2_grand-avg.mat");
load(cachePath);

%% Keep units of interest

% Load responsiveness
rTest = LMV.Resp.LoadPhaseResponseTest();
sTb = NP.Unit.AlignClusTb(ce.clusTb, rTest.clusTb, true);
ce.clusTb = [ce.clusTb sTb];

% Keep responsive mPrCG units
isRm = ce.clusTb.recId ~= "NP41_B1" | ce.clusTb.tId1 == "none";
ceSub = ce.Duplicate;
ceSub.RemoveUnits(isRm);

%% A static plot for the keynote slides

f = MPlot.Figure(26); clf
f.Position(3:4) = [700 550];
nRow = 10;
tl = tiledlayout(nRow, 1);
tl.Padding = 'compact';

ax = nexttile;
t = ceSub.GetArray('resp', 1);
tWin = t([1 end])';
colNames = {'cue1', 'stim', 'prod', 'cue3'};
cc = LMV.Param.GetTaskPhaseColors(colNames);
NP.TaskBaseClass.PlotEventWindows(ax, ceSub, 1, tWin, colNames, 'Colors', cc, 'Text', false);
% ax.Title.String = sprintf("%s neurons (n = %i)", "mPrCG", ceSub.numResp);
ax.XLabel.String = [];
ax.XTickLabel = [];
ax.YLabel.String = [];
ax.YTick = [];
axis(ax, 'off')

ax = nexttile([nRow-1 1]);
NP.UnitPlot.PethHeatmap(ceSub, [], 'SortVars', {'tId1', 'tId2', 'tId3', 'tR3'}, 'SortOrder', {'ascend', 'ascend', 'ascend', 'descend'}, ...
    'GroupColor', LMV.Param.GetTaskPhaseColors(rTest.phaseNames), 'Ribbon', false, 'Events', false);
ax.YLabel.FontWeight = 'bold';
ax.YLabel.String = sprintf("%s neurons (n = %i)", "mPrCG", ceSub.numResp);
ax.XLabel.String = "Aligned time (s)";
