%% Examine relashionships between functional, physiological, and anatomical properties of the units

anaDir = LMV.Data.GetAnalysisDir('fpa');

%% Load data

% Load physiological information of units
sPhys = load(fullfile(NP.Data.GetAnalysisRoot, 'units', 'computed_phys_clusTb.mat'));

% Load responsiveness
rTest = LMV.Resp.LoadPhaseResponseTest();

% Load M2 session PETHs with UMAP embedding coordinates
regions = ["all" LMV.Param.regions];
cePaths = fullfile(LMV.Data.GetAnalysisDir, "embedding", "umap", "computed_umap_"+regions+".mat")';
ceArray = NP.CE.LoadSession(cePaths);

% Add attributes to each of the clusTb
for i = 1 : numel(ceArray)
    clusTb = ceArray(i).clusTb;
    
    clusTb = [clusTb NP.Unit.AlignClusTb(clusTb, sPhys.clusTb, true)];
    
    [~, I] = MMath.SortLike(rTest.clusTb.clusId, clusTb.clusId);
    clusTb.tId1 = rTest.clusTb.tId1(I);
    clusTb = [clusTb rTest.sigTb(I,:)];
    
    ceArray(i).clusTb = clusTb;
end

%% Plot interactive embedding maps

colorBy = {'tId1', 'depth', 'region', 'recording', 'waveform', 'isi'};

for r = 1 : numel(regions)
    f = MPlot.Figure(2987); clf
    tl = tiledlayout(2,3);
    tl.Padding = 'compact';
    for i = 1 : numel(colorBy)
        NP.Embed.PlotScatter(ceArray(r), colorBy{i});
    end
    MPlot.Paperize(f, 'ColumnsWide', 2.3, 'ColumnsHigh', 1.2);
    exportgraphics(f, fullfile(anaDir, sprintf("umap_%s.png", regions(r))));
end

%% Compute relative depth distributions

phaseNames = rTest.phaseNames;
dData = cell(numel(phaseNames), numel(regions));

for r = 1 : numel(regions)
    clusTb = ceArray(r).clusTb;
    isResp = clusTb.tId1~="none";
    
    for i = 1 : numel(phaseNames)
        fprintf("Compute depth distribution of %s responsive units in %s\n", phaseNames(i), regions(r));
        isPhase = clusTb.tId1==phaseNames(i);
        s = LMV.FPA.ComputeRelativeDist(clusTb.depth(isPhase), clusTb.depth(isResp), Name="depth");
        s.region = regions(r);
        s.phaseName = phaseNames(i);
        dData{i,r} = s;
    end
end

%% 

f = MPlot.Figure(17504); clf
tl = tiledlayout(1, numel(regions));
tl.Padding = "compact";

for r = 1 : numel(regions)
    ax = nexttile;
    for i = 1 : numel(phaseNames)
        s = dData{i,r};
        c = LMV.Param.GetTaskPhaseColors(s.phaseName);
        if i == 1
            plot(s.yRef, s.x/1e3, Color='k', LineStyle='-', LineWidth=2); hold on
        end
        plot(s.yTest, s.x/1e3, Color=c, LineStyle='-', LineWidth=1);
    end
    
    pStr = cellfun(@(x) sprintf("p = %.2f", x.pval), dData(:,r));
    legend(ax.Children(end-1:-1:1), pStr, Location="southwest");
    
    ax.YDir = 'reverse';
    ax.Title.String = s.region;
    ax.XLabel.String = "Probability";
    ax.YLabel.String = "Depth (mm)";
end

%% 

f = MPlot.Figure(17505); clf
tl = tiledlayout(1, numel(regions));
tl.Padding = "compact";

for r = 1 : numel(regions)
    ax = nexttile;
    for i = 1 : numel(phaseNames)
        s = dData{i,r};
        c = LMV.Param.GetTaskPhaseColors(s.phaseName);
        plot(s.yRatio, s.x/1e3, Color=c, LineStyle='-', LineWidth=1); hold on
    end
    plot([1 1]', s.x([1 end])'/1e3, Color=[0 0 0 .3]);
    
    pStr = cellfun(@(x) sprintf("p = %.2f", x.pval), dData(:,r));
    legend(flip(ax.Children(2:end)), pStr, Location="southwest");
    
    ax.YDir = 'reverse';
    ax.Title.String = s.region;
    ax.XLabel.String = "Probability";
    ax.YLabel.String = "Depth (mm)";
end

return
%% Compute conditional probability distributions

phaseNames = rTest.phaseNames;
xVars = {phaseNames, phaseNames, phaseNames};
xNames = xVars;
xNames(1:3) = {'tId1'};

yVars = {'depth', 'wfId', 'recId', 'depth'};
yNames = yVars;

pData = cell(numel(xVars), numel(regions));

for i = 1 : numel(xVars)
    for r = 1 : numel(regions)
        fprintf("Compute P(%s|%s) for %s\n", xNames{i}, yNames{i}, regions(r));
        pData{i,r} = LMV.FPA.ComputeDist(ceArray(r).clusTb, xVars{i}, yVars{i}, 'nBoot', 1000);
    end
end

%% Plot distributions

for i = 1 : numel(xVars)
    f = MPlot.Figure(3251); clf
    tl = tiledlayout('flow');
    tl.Padding = 'loose';
    for r = 1 : numel(regions)
        ax = nexttile;
        LMV.FPA.PlotDist(pData{i,r}, 'Style', "heatmap");
        ax.Title.String = sprintf("%s (n = %i)", regions(r), ceArray(r).numResp);
        ax.XTickLabelRotation = 0;
    end
    MPlot.Paperize(f, 'ColumnsWide', 2.5, 'ColumnsHigh', .35);
    exportgraphics(f, fullfile(anaDir, "p_"+xNames{i}+"_given_"+yVars{i}+".png"));
end
