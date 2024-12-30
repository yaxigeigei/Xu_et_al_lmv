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

%% Plot dendrogram with labels

% Annotate clusters
groups = ["mirror", "bridge", "feedback", "other"];
switch mdlName
    case "smooth_lm"
        clusTb.hcGroup = repmat("other", height(clusTb), 1);
        clusTb.hcGroup(87:111) = "mirror";
        clusTb.hcGroup(112:127) = "bridge";
        clusTb.hcGroup(31:44) = "feedback";
        clusTb.hcGroup = categorical(clusTb.hcGroup, groups, 'Ordinal', true);
    otherwise
        hcInd = cluster(tree, 'Cutoff', 1.85, 'Depth', 4);
        clusTb.hcGroup = hcInd(leafOrder);
end
save(fullfile(anaDir, "linker_clusTb.mat"), 'clusTb');

% Mark clusters
cidMark = {520200153, 410100463, 350200795};
cidMark = {460100003, 410100463, 500300446};
ccMark = [0 1 0; 1 0 0; 0 0 1];
cidMark = {460100003, 520200153, 410100155, 450100286, 430100161, 430100549};
ccMark = brewermap(numel(cidMark), 'Set1');

% Make plot
f = MPlot.Figure(112); clf
tl = tiledlayout("flow");
tl.Padding = "compact";
LMV.Linker.LM.PlotDendrogram(sHC, clusTb, 'region', 'Ribbon', 'hcGroup');
% LMV.Linker.LM.PlotDendrogram(sHC, clusTb, 'region', 'Ribbon', 'hcGroup', 'MarkUnit', cidMark, 'MarkColor', ccMark);
% MPlot.Paperize(f, .5, .3);
% exportgraphics(f, fullfile(anaDir, "hc_dendrogram.png"));
% exportgraphics(f, fullfile(anaDir, "hc_dendrogram.pdf"));

%% Plot group kernel overlay

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

%% Quantify kernel peaks

statsFile = fullfile(anaDir, "kernel_peak_stats.txt");
fileID = fopen(statsFile, "w");

% Peak stats within each group
pk = struct;
for g = ["mirror", "feedback"]
    m = clusTb.hcGroup==g;
    [~, pk.(g)] = cellfun(@(x) findpeaks(x.Beta, x.dt, NPeaks=1, SortStr="descend"), clusTb.linker(m));
    p = signrank(pk.(g));
    fprintf(fileID, "%s kernel peak time stats:\n", g);
    fprintf(fileID, "Mean ± SD: %.2f ± %.2f. Median (IQR): %.2f (%.2f-%.2f)\n", mean(pk.(g)), std(pk.(g)), median(pk.(g)), prctile(pk.(g),25), prctile(pk.(g),75));
    fprintf(fileID, "Different from 0s with p-val = %e\n\n", p);
end

% Peak difference between mirror and feedback
p = ranksum(pk.mirror, pk.feedback);
fprintf(fileID, "mirror and feedback kernel peak time is different at p-val = %e\n", p);

fclose (fileID);

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
