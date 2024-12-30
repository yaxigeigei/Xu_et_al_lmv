%% Plot SCA dynamics

anaDir = LMV.Data.GetAnalysisDir("pop_dynamics", "sca");
srcTb = LMV.Data.FindSource([]);

%% Load data

% Load SCA results
load(fullfile(anaDir, "regTb.mat"), "regTb");
region = "mPrCG";
nComp = 12;
regTb = regTb(regTb.cond==region & regTb.nComp==nComp, :);

% Unit info
clusTb = py.pandas.read_pickle(fullfile(anaDir, "df", "clus.pkl"));
clusTb = table(clusTb);
isRegion = clusTb.region == region;
isResp = any(clusTb{:,4:end-1}, 2);
clusTb = clusTb(isRegion & isResp, :);

% Sentence-averaged responses
load(fullfile(LMV.Data.GetAnalysisDir, "data", "ce_m2_ex3_sentence-avg.mat"), 'ce');

% %% Load SCA data
% 
% % Source data
% load(fullfile(LMV.Data.GetAnalysisDir, "pop_dynamics", "ce_m2_ex3_sentence-avg.mat"), 'ce');
% 
% % Unit info
% clusTb = py.pandas.read_pickle(fullfile(anaDir, "df", "clus.pkl"));
% clusTb = table(clusTb);
% isRegion = clusTb.region == region;
% isResp = any(clusTb{:,4:end-1}, 2);
% clusTb = clusTb(isRegion & isResp, :);
% 
% % Load SCA output
% scaTbs = LMV.SCA.LoadResults(fullfile(anaDir, "computed_sca", "sca_"+region+".mat"), nComp);
% regTb = table;
% regTb.region = region;
% regTb = [regTb vertcat(scaTbs{:})];
% regTb = LMV.SCA.EnrichResultTable(regTb);
% 

%% Find preferred sentences for SCA components

% Sort components
[~, veInd] = sort(regTb.relVE{1}, 'descend');
switch region
    case 'mPrCG'
        veInd = veInd([1:3 5 4 6:end]);
    otherwise
end
Z = regTb.Z{1}(:,veInd);
V = regTb.V{1}(veInd,:);
CI = regTb.ZnullCI{1}(:,veInd,:);

[indAct, stimText] = LMV.SCA.FindActivatedSentences2(Z, CI, ce);

%% Load linker data

mdlName = "smooth_lm";
lkTb = LMV.Linker.LoadClusTb(mdlName);
lkTb(lkTb.hcGroup=="other",:) = [];

% load(fullfile(LMV.Data.GetAnalysisDir, "linker", mdlName, "linked_positions.mat"), "posTb");

%% Find best linked sentence for each linker unit

for i = 1 : height(lkTb)
    % pTb = posTb(posTb.clusId==lkTb.clusId(i), :);
    pTb = lkTb.posTb{i};
    if height(pTb) && any(pTb.isOut)
        lkTb.stimText{i} = pTb.stim(1:min(3,end)).GetParentLabel;
    else
        lkTb.stimText{i} = string(missing);
    end
end

%% 

f = MPlot.Figure(448); clf
tl = tiledlayout(3,4);
tl.Padding = "compact";

types = LMV.Linker.types;
% types = "bridge";
nType = numel(types);
vAll = cell(nComp,1);
vLkAll = cell(nComp, nType);

zEdges = linspace(-3, 8, 20);
hArgs = {'Normalization', 'probability', 'EdgeColor', 'none'};
% hArgs = {'Normalization', 'cdf', 'DisplayStyle', 'stairs'};

for i = 5 : nComp
    % 
    v = V(i,:)';
    v = zscore(v);
    vAll{i} = v;
    
    ax = nexttile;
    histogram(v, zEdges, hArgs{:}, 'FaceColor', 'k', 'FaceAlpha', 0.3); hold on
    
    for j = 1 : nType
        isType = lkTb.hcGroup==types(j);
        isSen = cellfun(@(x) any(ismember(x, stimText{i})), lkTb.stimText);
        isLk = ismember(clusTb.clusId, lkTb.clusId(isType & isSen));
        vLk = v(isLk);
        vLkAll{i,j} = vLk;
        
        if j ~= 2
            continue
        end
        
        if isempty(vLk)
            pval = NaN;
        else
            [~, pval] = kstest2(v, vLk);
        end
        histogram(vLk, zEdges, hArgs{:}, 'FaceColor', lines(1));
        MPlot.Axes(ax);
        title(ax, sprintf("Comp %i, %s, p=%.3f", i, types(j), pval));
    end
end

v = cat(1, vAll{:});

for j = 1 : nType
    ax = nexttile;
    
    histogram(v, zEdges, hArgs{:}, 'FaceColor', 'k', 'FaceAlpha', 0.3); hold on
    
    vLk = cat(1, vLkAll{:,j});
    histogram(vLk, zEdges, hArgs{:}, 'FaceColor', lines(1));
    
    if isempty(vLk)
        pval = NaN;
    else
        [~, pval] = kstest2(v, vLk);
    end
    
    xlabel("Z-scored loading");
    ylabel("Probability");
    MPlot.Axes(ax);
    title(ax, sprintf("All, %s, p=%.3f", types(j), pval));
    
    fprintf("\n%s, all components:\n\tmean±sd %.1f±%.1f, p=%.3f, n = %i; all, n = %i\n", types(j), mean(vLk), std(vLk), pval, numel(vLk), numel(v));
end

MPlot.Paperize(f, 1.5, .8);
% exportgraphics(f, fullfile(anaDir, "linker_loading.png"));
% exportgraphics(f, fullfile(anaDir, "linker_loading.pdf"));

return
%% Plot weights by regions

clusTb = py.pandas.read_pickle(fullfile(anaDir, "df", "clus.pkl"));
clusTb = table(clusTb);
isRegion = clusTb.region == "mPrCG";
isResp = any(clusTb{:,4:end-1}, 2);
clusTb(~(isRegion & isResp),:) = [];

types = [LMV.Linker.types "other"];
clusTb.hcGroup(:) = "other";
for i = 1 : numel(types)
    isLk = ismember(clusTb.clusId, lkTb.clusId(lkTb.hcGroup==types(i)));
    clusTb.hcGroup(isLk) = types(i);
end

f = MPlot.Figure(448); clf
tl = tiledlayout(3,4);
tl.Padding = "compact";
for i = 1 : nComp
    v = V(i,:)';
    ax = nexttile;
    boxplot(v, clusTb.hcGroup);
    % ax.XTick = 1 : numel(regions);
    % ax.XTickLabel = regions;
    ax.YLim = [-.2 .3];
    MPlot.Axes(ax);
    title(ax, "Component "+i);
end

%% Find high loading units

for regIdx = 1 %: numel(regions)
    region = region(regIdx);
    regTb = scaTbs{regIdx};
    k = find(regTb.nComp==nComp, 1);
    
    ve = regTb.varExplained{k} ./ nullVE{regIdx};
    [~, veInd] = sort(ve, 'descend');
    
    
    V = regTb.V{k}(veInd,:);
    % V = rotatefactors(V', 'Method', 'orthomax')';
    
    th = prctile(abs(V(:)), 95);
    isHigh = abs(V) > th;
    [~, uInd] = find(isHigh);
    uInd = unique(uInd);
    isUnit = MMath.Ind2Logical(uInd, size(V,2));
end

%% 

V = regTb.V{k}(veInd,:);
% V = rotatefactors(V', 'Method', 'orthomax')';

th = prctile(abs(V(:)), 95);
isHigh = abs(V) > th;
[~, uInd] = find(isHigh);
uInd = unique(uInd);
isUnit = MMath.Ind2Logical(uInd, size(V,2));

%% 

f = MPlot.Figure(815); clf
tl = tiledlayout(numel(region), 1);
tl.Padding = "compact";
for i = 1 : numel(region)
    m = clusTb.region==region(i) & isUnit;
    ax = nexttile;
    x = (1 : size(V,1))';
    plot(V(x,m), '-', 'Color', [0 0 0 .1]);
    title(region(i));
    xlabel("Components");
    xlim(x([1 end])' + [-1 1]);
    ax.Box = "off";
end

%% 

V = regTb.V{k}(veInd,:);

clusTb = py.pandas.read_pickle(fullfile(anaDir, "df", "clus.pkl"));
clusTb = table(clusTb);
clusTb(~any(clusTb{:,4:end},2),:) = [];

region = LMV.Param.regions;
recIds = unique(clusTb.recId, 'stable');
nComp = size(V,1);

f = MPlot.Figure(449); clf
tl = tiledlayout(3,4);
tl.Padding = "compact";
for i = 1 : nComp
    v = V(i,:)';
    th = prctile(v, 95);
    v(v < th) = NaN;
    
    [G, GID] = findgroups(clusTb.region);
    [m, sd] = splitapply(@(x) MMath.MeanStats(x), v, G);
    [~, I] = MMath.SortLike(GID, region);
    m = m(I);
    sd = sd(I);
    
    ax = nexttile;
    bar(m); hold on
    errorbar(m, sd);
    ax.XTick = 1 : numel(region);
    ax.XTickLabel = region;
    ax.Box = "off";
    title(ax, "Component "+i);
end

