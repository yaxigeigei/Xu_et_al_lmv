%% Analyze latent activity

anaDir = LMV.Data.GetAnalysisDir("pop_dynamics", "sca");
srcTb = LMV.Data.FindSource([]);

%% Load data

% Load SCA results
load(fullfile(anaDir, "regTb.mat"), "regTb", "clusTb");
regions = regTb.region;
nComp = regTb.nComp(1);

% Sentence-averaged responses
load(fullfile(LMV.Data.GetAnalysisDir, "data", "ce_m2_ex3_sentence-avg.mat"), 'ce');

%% Plot latent in original order

for i = 1 : height(regTb)
    f = MPlot.Figure(84774+i); clf(f);
    tl = tiledlayout(3, 4, 'Parent', f);
    tl.Padding = "compact";
    LMV.SCA.PlotLatentNoCI(ce, regTb.Z{i}, regTb.relVE{i}, tl);
    MPlot.Paperize(f, 'ColumnsWide', 1.4, 'ColumnsHigh', .6);
end

%% Plot latent in custom order

for i = 1 : height(regTb)
    ve = regTb.relVE{i};
    [~, I] = sort(ve, 'descend');
    switch regTb.cond(i)
        case 'mPrCG'
            I = I([1:3 5 4 6:end]);
    end
    
    f = MPlot.Figure(84774+i); clf(f);
    tl = tiledlayout(3, 4, 'Parent', f);
    tl.Padding = "compact";
    LMV.SCA.PlotLatent(ce, regTb.Z{i}(:,I), regTb.ZnullCI{i}(:,I,:), ve(I), tl);
    MPlot.Paperize(f, 'ColumnsWide', 1.4, 'ColumnsHigh', .6);
    exportgraphics(f, fullfile(anaDir, sprintf("sca_latent_%s.png", regTb.cond(i))));
    exportgraphics(f, fullfile(anaDir, sprintf("sca_latent_%s.pdf", regTb.cond(i))), ContentType="vector");
end

%% Plot basis orthogonality

f = MPlot.Figure(84794);

for i = 1 : height(regTb)
    cond = regTb.cond(i);
    [~, I] = sort(regTb.relVE{i}, 'descend');
    V = regTb.V{i}(I,:);
    M = V * V';
    
    clf(f);
    h = heatmap(M);
    h.Title = cond;
    h.XLabel = "Component";
    h.YLabel = "Component";
    h.FontSize = 12;
    h.Interpreter = "none";
    MPlot.Paperize(f, 'ColumnsWide', 0.4, 'ColumnsHigh', 0.35);
    exportgraphics(f, fullfile(anaDir, sprintf("sca_ortho_%s.png", cond)));
    exportgraphics(f, fullfile(anaDir, sprintf("sca_ortho_%s.pdf", cond)));
end

%% Quantify and plot sentence activations

f = MPlot.Figure(83754); clf
tl = tiledlayout("flow");
tl.Padding = "compact";
ax = nexttile;
for i = 1 : height(regTb)-1
    region = regTb.region(i);
    saTb = LMV.SCA.QuantifySentenceActivations(ce, regTb.Z{i}, regTb.ZnullCI{i});
    
    hh(i) = plot(saTb.span, saTb.relVar, '.', 'Color', LMV.Param.GetRegionColors(region)); hold on
    
    m = saTb.isAct;
    plot(saTb.span(m), saTb.relVar(m), 'o', 'Color', LMV.Param.GetRegionColors(region));
end
lgd = legend(hh(1:4), regions(1:4), Location="eastoutside");
ax.YLabel.String = "Relative variance";
ax.XLabel.String = "Span (frac. of trial)";
MPlot.Axes(ax);

MPlot.Paperize(f, .75, .3);
exportgraphics(f, fullfile(anaDir, "sca_activations.png"));
exportgraphics(f, fullfile(anaDir, "sca_activations.pdf"));

%% Test difference in activation during stim

[~, regionInd] = MMath.SortLike(regions, ["mPrCG", "vPrCG"], false);
stimCompInd = [4 6]; % mPrCG comp 4 and vPrCG comp 6

% Extract single-trial activations
actCell = cell(size(stimCompInd));
for i = regionInd'
    actCell{i} = LMV.SCA.ExtractActivations(ce, regTb.Ztr{i}(:,stimCompInd(i),:), "stim");
end

% Test significance
[~, p, ks2stat] = kstest2(actCell{1}(:), actCell{2}(:));

% Plot histograms with statistics
f = MPlot.Figure(47833); clf
histogram(actCell{1}(:), 20, 'FaceColor', LMV.Param.GetRegionColors(regions(1))); hold on
histogram(actCell{2}(:), 20, 'FaceColor', LMV.Param.GetRegionColors(regions(2)));
legend(regions(1:2));
ax = MPlot.Axes(gca);
ax.Title.String = sprintf("Mean phase activation (p = %.1e, D* = %.2f)", p, ks2stat);
ax.XLabel.String = "Mean activation (AU)";
ax.YLabel.String = "# of samples";
MPlot.Paperize(f, 1, .7);
exportgraphics(f, fullfile(anaDir, "mean_stim_activation_mPrCG-vPrCG.png"));

%% Test difference in sustained sentence-selective activations

% Match components using correlation
[~, condInd] = MMath.SortLike(regTb.cond, ["mPrCG", "mPrCG_no-bridge"], false);
M = corr(regTb.Z{condInd});
II{1} = 1:nComp;
[~, II{2}] = max(M, [], 2);

f = MPlot.Figure(7151); clf
for i = 1 : numel(condInd)
    ax = nexttile;
    imagesc(M(:,II{i}));
    ax.XLabel.String = "mPrCG";
    if i == 1
        ax.YLabel.String = "mPrCG_no-bridge";
    else
        ax.YLabel.String = "mPrCG_no-bridge (matched)";
    end
    ax.YLabel.Interpreter = "none";
    ax.Title.String = "Pairwise Pearson's r";
end
MPlot.Paperize(f, 1.2, 0.5);
exportgraphics(f, fullfile(anaDir, "match_components.png"));

% Find components with sustained sentence selective activations
[indAct, stimText] = LMV.SCA.FindActivatedSentences2(regTb.Z{condInd(1)}, regTb.ZnullCI{1}, ce);
cInd = find(~cellfun(@isempty, indAct))';
sInd = [indAct{cInd}];

% Extract activations
actCell = cell(size(condInd));
for i = 1 : numel(condInd)
    k = condInd(i);
    A = LMV.SCA.ExtractActivations(ce, regTb.Z{k}(:,II{i}));
    actCell{i} = arrayfun(@(a,b) A(a,b), sInd, cInd);
end
actCell{:}

% Test statistical significance
actRedu = (actCell{1} - actCell{2})./actCell{1};
actReduPrct = prctile(actRedu*100, [25 50 75]);
[p, ~, stats] = signrank(actCell{:});

% Make scatter plot and test
f = MPlot.Figure(7152); clf
plot([0 6], [0 6], 'Color', [0 0 0 .5]); hold on
plot(actCell{:}, 'o');
ax = MPlot.Axes(gca);
ax.XLabel.String = "mPrCG";
ax.YLabel.String = "mPrCG no-bridge";
ax.Title.String = sprintf("Reduction %.0f (%.0f-%.0f)%%, median (IQR)\nWilcoxon signed rank test p = %.1e", ...
    actReduPrct(2), actReduPrct(1), actReduPrct(3), p);
MPlot.Paperize(f, 0.5, 0.5);
exportgraphics(f, fullfile(anaDir, "mPrCG_no-bridge_latent_mag.png"));

return
%% Quantify and plot phase activations

f = MPlot.Figure(83775); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";
for i = 1 : numel(regions)
    ax = nexttile;
    region = regions(i);
    tb = compTbs{i};
    k = find(tb.nComp==nComp, 1);
    
    ve = tb.varExplained{k} ./ nullVE{i};
    [~, I] = sort(ve, 'descend');
    
    saTb = LMV.SCA.QuantifyPhaseActivations(ce, tb.Z{k}(:,I));
    
    bb = bar(saTb{end:-1:1,2:end}', 'stacked', 'w'); hold on
    
    ax.YLim(2) = 100;
    % ax.YLabel.String = "Variance";
    ax.YTickLabel = [];
    ax.XTickLabel = ["atten", "stim", "delay", "init", "prod"];
    ax.XTickLabelRotation = 90;
    ax.Title.String = region;
    MPlot.Axes(ax);
end

MPlot.Paperize(f, 1, .4);
exportgraphics(f, fullfile(anaDir, "sca_phase_activations.png"));

%% 

X = ssSCA(2).X;
sComp = ssSCA(2).decomp_12;
Z = X*sComp.U + sComp.b_u;

figure(1); clf
plot([Z(:,1) sComp.Z(:,1)])
plot([Z(:,1) sComp.Z(:,1)])

figure(2); clf
histogram(Z(:,1)-sComp.Z(:,1))
