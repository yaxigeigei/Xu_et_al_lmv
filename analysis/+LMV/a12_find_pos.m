%% Distributions of linking positions

mdlName = LMV.Linker.currentModel;
anaDir = fullfile(LMV.Data.GetAnalysisDir, "linker", mdlName);
srcTb = LMV.Data.FindSource([]);

%% Load cached data

% Profiles
cacheDir = LMV.Data.GetAnalysisDir("linker", mdlName, "computed_profile");
sPf = cell(size(srcTb.recId));
for recIdx = 1 : height(srcTb)
    recId = NP.SE.GetID(srcTb.name{recIdx});
    disp(recId);
    cacheFile = fullfile(cacheDir, sprintf("%s.mat", recId));
    sPf{recIdx} = load(cacheFile, 'ce', 'recId');
end

% Linker clusTb with M2 sentence PETHs
load(fullfile(LMV.Data.GetAnalysisDir, "linker", mdlName, "linker_clusTb_peth.mat"), "clusTb");

%% Update linked positions

for recIdx = 1 : numel(sPf)
    s = sPf{recIdx};
    LMV.Linker.LM.ComputeLinkingScores(s.ce);
    s.posTb = LMV.Linker.LM.FindLinkedPositions(s.ce);
    sPf{recIdx} = s;
end

%% Test if linked positions are sentence-selective

posTbs = cell(size(sPf));
for j = 1 : numel(sPf)
    disp(sPf{j}.recId);
    posTb = LMV.Linker.LM.IsPosSelective(sPf{j}.posTb, clusTb);
    posTbs{j} = posTb(posTb.isSelective,:);
end
posTb = cat(1, posTbs{:});

%% Cache results

save(fullfile(anaDir, "linked_positions.mat"), "posTb");

%% Organize posTb by linker types

typeTb = table;
typeTb.types = LMV.Linker.types';
typeTb.clusId = arrayfun(@(x) clusTb.clusId(clusTb.hcGroup==x), typeTb.types, 'Uni', false);
for i = 1 : height(typeTb)
    isLinker = ismember(posTb.clusId, typeTb.clusId{i});
    typeTb.posTb{i} = posTb(isLinker & posTb.score > LMV.Linker.scoreTh,:);
end

%% Plot linking positions for individual units

f = MPlot.Figure(310413); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";
for i = 1 : height(typeTb)
    ax = nexttile;
    LMV.Linker.LM.PlotLinkPositions(ax, typeTb.posTb{i});
    ax.XLim(1) = 0;
    title(sprintf("%s positions", typeTb.types(i)));
end
MPlot.Paperize(f, 2.2, .4);
exportgraphics(f, fullfile(anaDir, "linked_positions.png"));
print(f, '-vector', fullfile(anaDir, "linked_positions.pdf"), '-dpdf');

%% Plot time histograms of relative linking positions

f = MPlot.Figure(310424); clf
tl = tiledlayout("vertical");
tl.Padding = "compact";
for i = 1 : height(typeTb)
    ax = nexttile;
    LMV.Linker.LM.PlotLinkPositionHist(ax, typeTb.posTb{i});
    title(sprintf("%s positions", typeTb.types(i)));
end
MPlot.Paperize(f, .4, .7);
exportgraphics(f, fullfile(anaDir, "linked_pos_histograms.png"));
exportgraphics(f, fullfile(anaDir, "linked_pos_histograms.pdf"));

fileID = fopen(fullfile(anaDir, "uniform_dist_test.txt"), "w");
for i = 1 : height(typeTb)
    p = LMV.Linker.LM.TestUniformDist(typeTb.posTb{i});
    fprintf(fileID, "Test %s linking positions is uniform, p = %.2e\n", typeTb.types(i), p);
end
fclose(fileID);

%% Plot histograms showing the number of linking positions per linked sentence

f = MPlot.Figure(310433); clf
tl = tiledlayout(height(typeTb), 2);
tl.Padding = "compact";

for i = 1 : height(typeTb)
    ax = nexttile;
    LMV.Linker.LM.PlotNumLinkHist(ax, typeTb.posTb{i}, 'position');
    title(sprintf("%s", typeTb.types(i)));
    
    ax = nexttile;
    LMV.Linker.LM.PlotNumLinkHist(ax, typeTb.posTb{i}, 'sentence');
    title(sprintf("%s", typeTb.types(i)));
end

MPlot.Paperize(f, 1, .7);
exportgraphics(f, fullfile(anaDir, "link_num_histograms.png"));

return
%% Plot unit number heatmap of linking positions

stimIdList = LMV.Param.stimIdList14;

f = MPlot.Figure(310425); clf
tl = tiledlayout("vertical");
tl.Padding = "compact";

posTb = typeTb.posTb{1};
N = cell(size(stimIdList));
senLb = string();
for i = 1 : numel(stimIdList)
    isSen = posTb.stimId == stimIdList(i);
    sen = posTb.stim(find(isSen,1));
    wrd = Cut(sen);
    ph = Cut(wrd);
    senLb(i) = wrd(1:2).GetAllParentLabel + " ...";
    tEdges = [double(ph); ph(end).T.tmax];
    N{i} = histcounts(posTb.t0(isSen), tEdges);
end

L = cellfun(@numel, N);
[~, I] = sort(L);
maxL = max(L);
M = NaN(numel(stimIdList), maxL);
for i = 1 : numel(stimIdList)
    M(i, 1:numel(N{i})) = N{i};
end

ntArgs = MPlot.FindTileInd([1 3], 1, 1, 1);
ax = nexttile(ntArgs{:});
b = bar(mean(M, 'omitmissing'), 'histc');
b.EdgeColor = "none";
b.FaceColor = [0 0 0]+.5;
ax.XTick = [];
MPlot.Axes(ax);
ax.XLim = [1 maxL];

ntArgs = MPlot.FindTileInd([1 3], 1, 2, 1);
ax = nexttile(ntArgs{:});
hm = heatmap(1:maxL, senLb(I), M(I,:));
hm.MissingDataColor = 'w';
hm.GridVisible = "off";
hm.XLabel = "Position mirrored (in phoneme index)";
MPlot.Paperize(f, 1.2, 0.6);
