%% Anatomical and electrophysiological correlations with firing sparsity

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'sparsity');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end
srcTb = LMV.Data.FindSource('lmv');

%% Load data

uid = [410100245 460100005 450100286 450100323 410100450];
uData = NP.Unit.LoadUnitCache(uid, 'DataSource', 'm1');
se = LMV.SE.UnitCache2SE(uData);

se.userData.expInfo.subjectId = 'NP00';
se.userData.expInfo.blockId = 'B0';
se.userData.recMeta.Region = '';
se.userData.recMeta.Subregion = '';

%% Compute PETH

senTb = se.SplitConditions(["stimId" "stimText"], 'taskValue');

for s = 1 : height(senTb)
    fprintf("\n%s\n", senTb.stimText(s));
    senTb.ce(s) = LMV.SE.ComputeSessionPETH(senTb.se(s));
end

%% Plotting

stimId = LMV.Param.stimIdList4(3);
[~, k] = MMath.SortLike(senTb.stimId, stimId);

nUnits = numel(uid);

for i = 1 : nUnits
    f = MPlot.Figure(34580); clf
    
    panels = {'phone', 'rasters', @LMV.PlotSE.Scalogram}';
    panelArgs = {{}, {i, 'PETH', true}, {i}}';
    rowDist = [1 1 2];
    
    axArray = NP.TaskBaseClass.PlotSE(senTb.se(k), panels, rowDist, 'PanelArgs', panelArgs);
    for j = 1 : numel(axArray)-1
        axArray(j).XTickLabel = [];
    end
    
    MPlot.Paperize(f, 'ColumnsWide', 1.2, 'ColumnsHigh', .6);
    exportgraphics(f, fullfile(anaDir, sprintf("%s_u%i_scalogram.png", stimId, uid(i))));
end


return
%% Plot figure

[~, k] = MMath.SortLike(senTb.stimId, LMV.Param.stimIdList4(3));

uInd = 1:3;
nUnits = numel(uInd);

f = MPlot.Figure(34580); clf

panels = {'phone', 'rasters', @LMV.UnitPlot.PETH, @LMV.UnitPlot.Scalogram}';
panels = repmat(panels, [1 nUnits]);

panelArgs = num2cell(num2cell(uInd));
panelArgs = repmat(panelArgs, [size(panels,1) 1]);
panelArgs(1,:) = {{}};

rowDist = [1 1 1 2];

NP.TaskBaseClass.PlotSE(senTb.se(k), panels, rowDist, 'PanelArgs', panelArgs);
