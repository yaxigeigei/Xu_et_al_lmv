%% Hierarchical clustering of linker kernels

mdlName = LMV.Linker.currentModel;
anaDir = LMV.Data.GetAnalysisDir("linker", mdlName);
srcTb = LMV.Data.FindSource([]);

%% Load data

clusTbs = LMV.Linker.LM.LoadModels(mdlName, srcTb.recId);
clusTbCat = cat(1, clusTbs{:});

%% Select units with significant models

clusTb = clusTbCat;

noMdl = cellfun(@isempty, clusTb.linker);
clusTb(noMdl,:) = [];

isNS = cellfun(@(x) x.null.r2Pval, clusTb.linker) > 0.05;
clusTb(isNS,:) = [];

%% Hierarchical clustering

% Get kernel weights
mdls = clusTb.linker;
Beta = cellfun(@(x) x.Beta, mdls, 'Uni', false);
Beta = cat(2, Beta{:});

% Compute hierarchical clustering
D = pdist(Beta', "cosine");
tree = linkage(D, "average");
leafOrder = optimalleaforder(tree, D, 'Criteria', 'adjacent');
sHC = struct;
sHC.tree = tree;
sHC.leafOrder = leafOrder;
clusTb = clusTb(sHC.leafOrder,:);

%% Inspect and annotate groups on dendrogram

% Annotate groups
groups = ["mirror", "bridge", "feedback", "other"];
clusTb.hcGroup = repmat("other", height(clusTb), 1);
clusTb.hcGroup(87:111) = "mirror";
clusTb.hcGroup(112:127) = "bridge";
clusTb.hcGroup(31:44) = "feedback";
clusTb.hcGroup = categorical(clusTb.hcGroup, groups, 'Ordinal', true);

% Save clusTb
clusTbPath = fullfile(anaDir, "linker_clusTb.mat");
if ~isfile(clusTbPath)
    save(fullfile(anaDir, "linker_clusTb.mat"), 'clusTb');
else
    disp("linker_clusTb.mat already exists. Not overwriting.");
    return
end

% Mark units
cidMark = {520200153, 410100463, 350200795};
cidMark = {460100003, 410100463, 500300446};
ccMark = [0 1 0; 1 0 0; 0 0 1];
cidMark = {460100003, 520200153, 410100155, 450100286, 430100161, 430100549};
ccMark = brewermap(numel(cidMark), 'Set1');

% Plot dendrogram
f = MPlot.Figure(112); clf
tl = tiledlayout("flow");
tl.Padding = "compact";
LMV.Linker.LM.PlotDendrogram(sHC, clusTb, 'region', 'Ribbon', 'hcGroup');
% LMV.Linker.LM.PlotDendrogram(sHC, clusTb, 'region', 'Ribbon', 'hcGroup', 'MarkUnit', cidMark, 'MarkColor', ccMark);
MPlot.Paperize(f, .5, .3);
exportgraphics(f, fullfile(anaDir, "hc_dendrogram.png"));
exportgraphics(f, fullfile(anaDir, "hc_dendrogram.pdf"));

%% Plot group kernel overlays

f = MPlot.Figure(116); clf
tl = tiledlayout(3, 1);
tl.Padding = "compact";
for g = groups(1:3)
    ax = nexttile;
    LMV.Linker.LM.PlotKernelOverlay(clusTb.linker(clusTb.hcGroup==g));
    ax.Title.String = [upper(g{1}(1)) g{1}(2:end) ' kernels'];
end
MPlot.Paperize(f, .3, .8);
exportgraphics(f, fullfile(anaDir, "hc_group_kernels.png"));
exportgraphics(f, fullfile(anaDir, "hc_group_kernels.pdf"));

%% Plot exmaple kernels

f = MPlot.Figure(115); clf
tl = tiledlayout(3, 2);
tl.Padding = "compact";
for c = [cidMark{:}]
    ax = nexttile;
    LMV.Linker.LM.FigExampleKernel(clusTb.linker{clusTb.clusId==c});
end
MPlot.Paperize(f, .7, .8);
exportgraphics(f, fullfile(anaDir, "example_kernels.png"));
exportgraphics(f, fullfile(anaDir, "example_kernels.pdf"));

%% Identify clusters with similar cutoffs

% Convert leaf indices to original order (before leaf reordering)
linkerInd = { ...
    find(clusTb.hcGroup == "mirror"), ...
    find(clusTb.hcGroup == "bridge"), ...
    find(clusTb.hcGroup == "feedback")};
linkerOrigInd = cellfun(@(x) leafOrder(x), linkerInd, 'Uni', false);

% Find the cutoff height for a given cluster
function h = FindCutoff(tree, nodeInd)
    m = size(tree, 1) + 1;
    nodeInd = nodeInd(:);
    I = [];
    while true
        M = ismember(tree(:,1:2), nodeInd);
        lastI = I;
        I = find(sum(M,2) == 2);
        if isempty(I)
            h = tree(lastI,3);
            break
        end
        isLinked = ismember(nodeInd, tree(I,:));
        nodeInd = [nodeInd(~isLinked); m+I];
    end
end

% Use the mean height to cut clusters
cutoffs = cellfun(@(x) FindCutoff(tree, x), linkerOrigInd);
meanCutoff = mean(cutoffs);

fprintf('Group cutoffs: Mirror=%.3f, Bridge=%.3f, Feedback=%.3f\n', cutoffs(1), cutoffs(2), cutoffs(3));
fprintf('Using mean cutoff = %.3f\n', meanCutoff);

hcInd = cluster(tree, 'Cutoff', meanCutoff, "criterion", "distance");
clusTbAuto = clusTb;
clusTbAuto.hcGroup = hcInd(leafOrder);

% Plot dendrogram
f = MPlot.Figure(113); clf
tl = tiledlayout("flow");
tl.Padding = "compact";
LMV.Linker.LM.PlotDendrogram(sHC, clusTbAuto, 'region', 'Ribbon', 'hcGroup');
MPlot.Paperize(f, 1, .3);
exportgraphics(f, fullfile(anaDir, "hc_dendrogram_auto.png"));
exportgraphics(f, fullfile(anaDir, "hc_dendrogram_auto.pdf"));

%% Plot auto-clustered group kernel overlays

% Get unique cluster IDs and their counts
uniqueClusters = unique(clusTbAuto.hcGroup);
clusterCounts = arrayfun(@(x) sum(clusTbAuto.hcGroup == x), uniqueClusters);

% Filter out small clusters
validClusters = uniqueClusters(clusterCounts > 5);
qt = prctile(clusterCounts(clusterCounts <= 5), [25 50 75]);
fprintf("Small clusters have median size of %i (IQR %i-%i) units\n", qt(2), qt(1), qt(3));

% Create figure for auto-clustered kernels
f = MPlot.Figure(117); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";

% Plot each valid cluster
for i = 1 : length(validClusters)
    ax = nexttile;
    clusterId = validClusters(i);
    LMV.Linker.LM.PlotKernelOverlay(clusTbAuto.linker(clusTbAuto.hcGroup == clusterId));
    ax.Title.String = sprintf('Cluster %d (n=%d)', clusterId, clusterCounts(clusterId));
    ax.YLim = [-2.5 4];
    if i > 1
        ax.YLabel.String = [];
        ax.YTickLabel = [];
    end
end

MPlot.Paperize(f, 1.5, .3);
exportgraphics(f, fullfile(anaDir, "hc_group_kernels_auto.png"));
exportgraphics(f, fullfile(anaDir, "hc_group_kernels_auto.pdf"));

return
%% Examine clustering with selected units

% Find units to mark
recIds = unique(clusTb.recId);
cidMark = cell(1,3);
cidOut = cell(1,3);

cid = LMV.Linker.GetBridgeClusId(recIds);
isPlotted = ismember(cid, clusTb.clusId);
fprintf("\nMarked %i/%i bridge units.\n", sum(isPlotted), numel(cid));
disp("Units not plotted:");
disp(cid(~isPlotted));
cidMark{1} = cid;
cidOut{1} = cid(~isPlotted);

cid = LMV.Linker.GetMirrorClusId(recIds);
isPlotted = ismember(cid, clusTb.clusId);
fprintf("\nMarked %i/%i mirror units.\n", sum(isPlotted), numel(cid));
disp("Units not plotted:");
disp(cid(~isPlotted));
cidMark{2} = cid;
cidOut{2} = cid(~isPlotted);

cid = LMV.Linker.GetFeedbackClusId(recIds);
isPlotted = ismember(cid, clusTb.clusId);
fprintf("\nMarked %i/%i feedback units.\n", sum(isPlotted), numel(cid));
disp("Units not plotted:");
disp(cid(~isPlotted));
cidMark{3} = cid;
cidOut{3} = cid(~isPlotted);

ccMark = [1 0 0; 0 1 0; 0 0 1];

%% 

f = MPlot.Figure(99); clf
LMV.Overview.SessionFromCache([cidOut{:}], 'DataSource', 'm2');

%% Make plot

f = MPlot.Figure(112); clf
f.Position(3:4) = [1600 200];
tl = tiledlayout("flow");
tl.Padding = "compact";
LMV.Linker.LM.PlotDendrogram(sHC, clusTb, 'region', 'MarkUnit', cidMark, 'MarkColor', ccMark);
