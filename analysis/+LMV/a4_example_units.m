%% Single-unit classification of sentences at each task phase

anaDir = LMV.Data.GetAnalysisDir("units", "examples");
srcTb = LMV.Data.FindSource([]);

%% Plot rasters of example units

% Specify units and sentences
uu = LMV.Fig.GetExampleInfo("Fig1G");

% Convert stimText to stimId
stimDict = LMV.Param.stimDict;
stimIdList = string(fieldnames(stimDict));
stimTextList = structfun(@(x) x, stimDict);
ind = arrayfun(@(x) find(startsWith(stimTextList, x), 1), uu.startText);
uu.stimId = stimIdList(ind);

disp(uu);

f = MPlot.Figure(120020); clf
nUnit = height(uu);
tl = tiledlayout(nUnit+1, 1);
tl.Padding = 'compact';

ax = nexttile;
uData = NP.Unit.LoadUnitCache(uu.clusId(1), 'DataSource', 'm2');
se = LMV.SE.UnitCache2SE(uData);
[tt, tv] = se.GetTable('taskTime', 'taskValue');
tWin = [tt.cue1On(1) tt.matchOff(1)] + [-1 1]*.2;
NP.TaskBaseClass.PlotEventWindows(ax, se, 1, tWin, {'cue1', 'stim', 'cue3', 'prod'}, 'Colors', LMV.Param.GetTaskPhaseColors(["none", "stim", "none", "prod"]));
ax.Visible = 'off';
ax.XTick = [];
ax.XLabel.String = [];

for i = 1 : nUnit
    ax = nexttile;
    LMV.Fig.UnitRaster(ax, uu.clusId(i), 'DataSource', 'm2', 'StimIdList', uu.stimId(i));
end
MPlot.Paperize(f, 'ColumnsHigh', 0.12*(nUnit+1), 'ColumnsWide', 1.8);
exportgraphics(f, fullfile(anaDir, "fig1_unit_rasters.png"));
exportgraphics(f, fullfile(anaDir, "fig1_unit_rasters.pdf"));

return
%% 

uData = NP.Unit.LoadUnitCache(uid, 'DataSource', 'm1');
se = LMV.SE.UnitCache2SE(uData);

%% 

f = MPlot.Figure(234); clf
LMV.Overview.SessionFromCache(uid, 'DataSource', 'm2', 'StimIdList', LMV.Param.stimIdList12);
MPlot.Paperize(f, 'ColumnsHigh', 1.2, 'ColumnsWide', 2.5, 'FontSize', 5);
exportgraphics(f, fullfile(anaDir, "example_session_rasters_s12.png"));

%% 

f = MPlot.Figure(345); clf
LMV.Overview.SessionFromCache(uid, 'DataSource', 'm2', 'StimIdList', LMV.Param.stimIdList4);
MPlot.Paperize(f, 'ColumnsHigh', 0.8, 'ColumnsWide', 2.5, 'FontSize', 5);
exportgraphics(f, fullfile(anaDir, "example_session_rasters_s4.png"));

%% 

f = MPlot.Figure(123); clf
LMV.Overview.SentencesFromCache(uid, 'StimIds', LMV.Param.stimIdList4, 'Features', "phone", 'UnitLabels', "depth", 'MinTaskHight', 20);
MPlot.Paperize(f, 'ColumnsWide', 2.5, 'ColumnsHigh', 0.8, 'FontSize', 5);
exportgraphics(f, fullfile(anaDir, "example_sentence_rasters_s4.png"));

%% 
