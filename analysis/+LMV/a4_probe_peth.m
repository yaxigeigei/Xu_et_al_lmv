%% Plot example single-column PETHs

anaDir = LMV.Data.GetAnalysisDir('examples');
srcTb = LMV.Data.FindSource('probe_peth');

%% Load data

% Load M2 PETHs
cachePath = fullfile(LMV.Data.GetAnalysisDir, "data", "ce_m2_grand-avg.mat");
load(cachePath, 'ce');
clusTb = ce.clusTb;

% Add probe ID
uidChar = num2str(clusTb.clusId);
clusTb.probeId = uidChar(:,5) - '0';

% Load signrank responsiveness
rTest = LMV.Resp.LoadPhaseResponseTest();
sTb = NP.Unit.AlignClusTb(clusTb, rTest.clusTb, true);
sTb = [clusTb sTb];
ce.clusTb = sTb;

%% Plot heatmaps by site

srcTb.probeId = [0 1 0]'; % only use one probe

f = MPlot.Figure(36); clf

rowDist = [1 8];
colDist = ones(1, height(srcTb));

tl = tiledlayout(sum(rowDist), sum(colDist));
tl.Padding = 'compact';

for i = 1 : height(srcTb)
    % Prepare ce
    recId = srcTb.recId{i};
    isRec = ce.clusTb.recId == recId & ce.clusTb.probeId == srcTb.probeId(i);
    isResp = ce.clusTb.tId1 ~= "none";
    ceSub = ce.Duplicate;
    ceSub.RemoveUnits(~isRec | ~isResp);
    
    % PETH heatmap
    ntArgs = MPlot.FindTileInd(rowDist, colDist, 2, i);
    ax = nexttile(ntArgs{:});
    NP.UnitPlot.PethHeatmap(ceSub, [], 'SortVars', {'tId1', 'tId2', 'tId3', 'tR3'}, 'SortOrder', {'ascend', 'ascend', 'ascend', 'descend'}, ...
        'Ribbon', true, 'GroupColor', LMV.Param.GetTaskPhaseColors(rTest.phaseNames), 'Events', false);
    ax.Title.String = sprintf("%s neurons, %s (n = %i)", srcTb.Region{i}, recId, ceSub.numResp);
    ax.Title.Interpreter = 'none';
    % ax.Title.FontWeight = 'bold';
    
    % Task phase indicators
    tWin = ax.XLim;
    ntArgs = MPlot.FindTileInd(rowDist, colDist, 1, i);
    ax = nexttile(ntArgs{:});
    % t = ceSub.GetArray('resp', 1);
    % tWin = t([1 end])';
    colNames = {'cue1', 'stim', 'prod', 'cue3'};
    cc = LMV.Param.GetTaskPhaseColors(colNames);
    NP.TaskBaseClass.PlotEventWindows(ax, ceSub, 1, tWin, colNames, 'Colors', cc, 'Text', false);
    ax.XLabel.String = [];
    ax.XTickLabel = [];
    ax.YLabel.String = [];
    ax.YTick = [];
    axis(ax, 'off')
end

MPlot.Paperize(f, 'ColumnsWide', 2.5, 'ColumnsHigh', 0.8);
exportgraphics(f, fullfile(anaDir, "probe_peth_heatmaps.png"));


return
%% Load all recordings

srcTb = NP.Data.FindSource('lmv');
cePaths = fullfile(NP.Data.GetAnalysisRoot, 'data', 'ce_m1_sentence-avg', srcTb.recId+"_ce_m1.mat");
ceArray = NP.CE.LoadSession(cePaths);
nRec = numel(ceArray);

%% Make probe PETH plots for all recordings

egStimId = "mfxv0_si1635"; % "we've got plenty of time to think about that"

for i = 1 : nRec
    f = MPlot.Figure(690); clf
    NP.UnitPlot.ProbePETH(ceArray(i), egStimId);
%     MPlot.Paperize(f, 'ColumnsWide', 1.8, 'ColumnsHigh', 1.2);
    
    figDir = fullfile(anaDir, "probe_peth");
    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end
%     exportgraphics(f, fullfile(figDir, egStimId+"-"+NP.SE.GetRegion(ce)+"-"+NP.SE.GetID(ce)+".png"));
end

%% Load example recordings

srcTb = NP.Data.FindSource('probe_peth');
cePaths = fullfile(NP.Data.GetAnalysisRoot, 'data', 'm1_sentence-avg_ce', srcTb.recId+"_ce_m1.mat");
ceArray = NP.CE.LoadSession(cePaths);

%% Make example probe PETH plot with all heatmaps in one figure

f = MPlot.Figure(690); clf
NP.UnitPlot.ProbePETH(ceArray, egStimId);
% MPlot.Paperize(f, 'ColumnsWide', 1.8, 'ColumnsHigh', 1.2);
% exportgraphics(f, fullfile(anaDir, "example_probe_peth-"+egStimId+".png"));

return
%% Utility: Check the most repeated sentences and pick one as example

for i = 1 : nRec
    ce = ceArray(i);
    tv = ce.GetTable('taskValue');
    disp(NP.SE.GetRegion(ce));
    disp(head(sortrows(tv, 'numTrial', 'descend'), 4));
end
