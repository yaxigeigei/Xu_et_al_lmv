%% 

anaDir = LMV.Data.GetAnalysisDir("pop_dynamics", "sca");
srcTb = LMV.Data.FindSource([]);

%% Load data

load(fullfile(LMV.Data.GetAnalysisDir, "data", "ce_m2_ex3_sentence-avg.mat"), 'ce');

regions = ["all" LMV.Param.regions];
ssSCA = arrayfun(@(x) load(fullfile(anaDir, "computed_sca", "sca_"+x+".mat")), regions);
ssPCA = arrayfun(@(x) load(fullfile(anaDir, "computed_pca", "pca_"+x+".mat")), regions);

%% Organize data into tables and add metrics of model fit

scaTbs = arrayfun(@LMV.SCA.MakeDecompTable, ssSCA, 'Uni', false);
pcaTbs = arrayfun(@LMV.SCA.MakeDecompTable, ssPCA, 'Uni', false);

%% Compute VE from null distributions

nComp = 20;
nInputDims = arrayfun(@(x) size(x.X,2), ssSCA);
nullVE = arrayfun(@(x) LMV.SCA.ComputeNullVE(x, nComp), nInputDims, 'Uni', false);

%% Plot PCA library

pcaDir = LMV.Data.GetAnalysisDir("pop_dynamics", "sca", "lib_pca");

f = MPlot.Figure(84774);
f.WindowState = "maximized";

for i = 1 : numel(pcaTbs)
    tb = pcaTbs{i};
    for j = height(tb)
        clf(f);
        ve = tb.varExplained{j} ./ nullVE{i}(1:tb.nComp(j));
        LMV.SCA.PlotLatent(ce, tb.Z{j}, ve);
        % return
        exportgraphics(f, fullfile(pcaDir, sprintf("PCA_%s_nc%i.png", regions(i), tb.nComp(j))));
    end
end

%% Quantify the percent VE by PCs

f = MPlot.Figure(84773); clf
ax = nexttile;

filePath = fullfile(anaDir, "pca_ve_stats.txt");
fileID = fopen(filePath, "w");

for i = 2 : numel(pcaTbs)
    tb = pcaTbs{i};
    x = 1 : 10;
    y = tb.varExplained{1}(x) / tb.varExplained{1}(1);
    plot(x, y, 'o-', 'Color', LMV.Param.GetRegionColors(regions(i))); hold on
    
    fprintf(fileID, "%s\tPC1/PC12 = %.2f\tPC2/PC1 = %.2f\n", regions(i), y(1)/sum(y), y(2));
end
legend(regions(2:end));
ax.XLim = [0.5 x(end)+.5];
ax.YLim = [0 1];
ax.XLabel.String = "Component #";
ax.YLabel.String = "Norm. VE";
ax.Title.String = "PCA norm. VE";
MPlot.Axes(ax);
MPlot.Paperize(f, .7, .5);
exportgraphics(f, fullfile(anaDir, "pca_norm_ve.png"));

fclose(fileID);

%% Plot SCA library

scaDir = LMV.Data.GetAnalysisDir("pop_dynamics", "sca", "lib_sca");

f = MPlot.Figure(84774);
f.WindowState = "maximized";

for i = 1 : numel(scaTbs)
    tb = scaTbs{i};
    for j = 1 : height(tb)
        clf(f);
        ve = tb.varExplained{j} ./ nullVE{i}(1:tb.nComp(j));
        LMV.SCA.PlotLatent(ce, tb.Z{j}, ve);
        exportgraphics(f, fullfile(scaDir, sprintf("SCA_%s_nc%i.png", regions(i), tb.nComp(j))));
    end
end

%% Plot basis orthogonality

scaDir = LMV.Data.GetAnalysisDir("pop_dynamics", "sca", "lib_sca");

f = MPlot.Figure(84794);
f.WindowState = "maximized";

for i = 1 : numel(scaTbs)
    tb = scaTbs{i};
    clf(f);
    for j = 1 : height(tb)
        M = tb.V{j} * tb.V{j}';
        ax = nexttile;
        h = heatmap(M);
        h.Title = sprintf("Pairwise cosine, nc = %i", tb.nComp(j));
    end
    exportgraphics(f, fullfile(scaDir, sprintf("Ortho_%s.png", regions(i))));
end

%% Plot total variance explained as a function of the number of componennts

f = MPlot.Figure(84785); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";

for i = 2 : numel(scaTbs)
    tb = scaTbs{i};
    VE1 = cellfun(@(x) sum(x), tb.varExplained);
    
    VE2 = cumsum(pcaTbs{i}.varExplained{1})';
    VE2 = VE2(6:end);
    
    VE = [VE1 VE2];
    
    tb.rVE = VE1 ./ VE2;
    disp(regions(i));
    disp(tb(:,["nComp", "rVE"]));
    
    ax = nexttile;
    plot(ax, tb.nComp, VE);
    legend(["SCA", "PCA"], "Location", "southeast", "Box", "off");
    ax.XTick = tb.nComp(1:2:end);
    ax.XTickLabelRotation = 0;
    ax.XLim = tb.nComp([1 end])' + [-1 1];
    ax.YLim(1) = 0;
    ax.Title.String = sprintf("Total VE, %s", regions(i));
    ax.XLabel.String = "# of components";
    ax.YLabel.String = "%";
    MPlot.Axes(ax);
end

MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', .3);
exportgraphics(f, fullfile(anaDir, "total_VE.png"));
exportgraphics(f, fullfile(anaDir, "total_VE.pdf"));

return
%% Plot relative variance explained as a function of the number of componennts

f = MPlot.Figure(84784); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";

for i = 2 : numel(scaTbs)
    tb = scaTbs{i};
    maxComp = max(tb.nComp);
    VE = cellfun(@(x) [x NaN(1,maxComp-numel(x))], tb.varExplained, 'Uni', false);
    VE = cat(1, VE{:});
    rVE = VE ./ nullVE{i}(1:maxComp);
    cc = flip(sky(maxComp+1));
    
    ax = nexttile;
    hh = plot(ax, tb.nComp, rVE); hold on
    for j = 1 : numel(hh)
        hh(j).Color = cc(j,:);
    end
    plot(ax, tb.nComp, mean(rVE,2,"omitmissing"), 'k', 'LineWidth', 2);
    
    ax.XTick = tb.nComp(1:2:end);
    ax.XTickLabelRotation = 0;
    ax.XLim = tb.nComp([1 end])' + [-1 1];
    ax.YLim(1) = 0;
    ax.Title.String = sprintf("SCA rel. VE, %s", regions(i));
    ax.XLabel.String = "# of components";
    ax.YLabel.String = "VE / null VE";
    MPlot.Axes(ax);
end

MPlot.Paperize(f, 'ColumnsWide', 1.6, 'ColumnsHigh', .3)
% exportgraphics(f, fullfile(anaDir, "sca_rel_VE.png"));

%% Check results with explicit orthogonality constraint

sSCAo = load(fullfile(anaDir, "computed_sca_orth", "sca_all.mat"));
scaOrthTb = LMV.SCA.MakeDecompTable(sSCAo);

%% 

scaDir = LMV.Data.GetAnalysisDir("sca", "lib_sca_orth");

f = MPlot.Figure(84774);
f.WindowState = "maximized";

tb = scaOrthTb;
for j = 1 : height(tb)
    clf(f);
    ve = tb.varExplained{j} / (tb.nComp(j)/size(tb.U{j},1));
    LMV.SCA.PlotLatent(ce, tb.Z{j}, ve);
    exportgraphics(f, fullfile(scaDir, sprintf("SCA_all_nc%i.png", tb.nComp(j))));
end
