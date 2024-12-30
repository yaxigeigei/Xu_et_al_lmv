%% Plot summary goodness-of-fit of all units and all models

anaDir = LMV.Data.GetAnalysisDir('trf');
srcTb = LMV.Data.FindSource([]);

%% Load test results

rTest = LMV.Resp.LoadPhaseResponseTest();
clusTb = rTest.clusTb;
rSig = rTest.sigTb{:,["stim" "prod" "prod"]};

sTest = LMV.Resp.LoadSentenceSelectTest();
sSig = sTest.sigTb{:,["stim" "prod" "prod"]};

%% Load data

sets = ["phone", "strf", "artic"]; % "combo1", "combo3pm", "combo3akt", 
targets = LMV.TRF.targets;

u = height(clusTb);
s = numel(sets);
t = numel(targets);

setMdls = cell(size(sets));
rrr = NaN(u, s, t);
ppp = rrr;
for i = 1 : s
    % Load models
    mdlNames = sets(i) + "_" + targets;
    mdlTbs = LMV.TRF.LoadModels(mdlNames);
    mdlTb = cat(1, mdlTbs{:});
    
    % Make unit order in mdlTb consistent with testTb
    [~, I] = MMath.SortLike(mdlTb.clusId, clusTb.clusId);
    mdlTb = mdlTb(I,:);
    
    % Get computed models
    mdls = mdlTb{:,mdlNames};
    hasMdl = ~cellfun(@isempty, mdls);
    hasNull = cellfun(@(x) isfield(x, 'null'), mdls);
    
    % Extract r
    rr = NaN(u,t);
    rr(hasMdl) = cellfun(@(x) x.r2, mdls(hasMdl));
    rr(rr<0) = 0;
    rr = sqrt(rr);
    rrr(:,i,:) = rr;
    
    % Extract pvals
    pp = NaN(u,t);
    pp(hasNull) = cellfun(@(x) x.null.r2Pval, mdls(hasNull));
    ppp(:,i,:) = pp;
    
    setMdls{i} = mdlTb;
end

%% Plot r-value CDFs

regions = unique(clusTb.region);

for k = 1 : numel(regions)
    f = MPlot.Figure(22200); clf
    tl = tiledlayout(t,s);
    tl.Padding = 'compact';
    
    isReg = clusTb.region == regions(k);
    
    for i = 1 : t
        for j = 1 : s
            ax = nexttile;
            
            % Responsive units
            r = rrr(isReg,j,i);
            m = rSig(isReg,i) > 0;
            r(~m) = NaN;
            histogram(r, 0:.001:1, 'Normalization', 'cdf', 'DisplayStyle', 'stairs'); hold on
            
            % And selective units
            r = rrr(isReg,j,i);
            m = rSig(isReg,i) > 0 & sSig(isReg,i) > 0;
            r(~m) = NaN;
            histogram(r, 0:.001:1, 'Normalization', 'cdf', 'DisplayStyle', 'stairs'); hold on
            
            % And r-value significant units
            r = rrr(isReg,j,i);
            m = rSig(isReg,i) > 0 & sSig(isReg,i) > 0 & ppp(isReg,j,i) < 0.05;
            r(~m) = NaN;
            histogram(r, 0:.001:1, 'Normalization', 'cdf', 'DisplayStyle', 'stairs'); hold on
            
            if i == 1
                ax.Title.String = sets(j);
            end
            if j == 1
                ax.YLabel.String = targets(i);
            end
            if i == t
                ax.XLabel.String = 'r-value';
            end
            if i == t && j == s
                ldg = legend(["responsive", "selective", "encoding"], 'Location', 'northwest');
            end
            ax.XLim = [0 .8];
            ax.YLim = [0 .6];
            ax.YTick = 0:.1:1;
            MPlot.Axes(ax);
        end
    end
    MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', .8);
    exportgraphics(f, fullfile(anaDir, sprintf("r_cdf_%s.png", regions(k))));
end

%% Plot boxplots

regions = LMV.Param.regions;
a = numel(regions);

tar = cellstr(targets);
tar = cellfun(@(x) x(1), tar, 'Uni', false);

f = MPlot.Figure(22300); clf
tl = tiledlayout(a,s);
tl.Padding = 'compact';

for k = 1 : a
    isReg = clusTb.region == regions(k);
    
    for j = 1 : s
        rr = rrr(isReg,j,:);
        rr = squeeze(rr);
        
        m = sSig(isReg,:) > 0;
        rr(~m) = NaN;
        
        ax = nexttile;
        boxplot(ax, rr, tar);
        ax.YLim = [-.05 .8];
        ax.Title.String = sets(j);
        MPlot.Axes(ax);
        
        if j == 1
            ax.YLabel.String = regions(k);
        end
    end
end

MPlot.Paperize(f, 'ColumnsWide', 1.2, 'ColumnsHigh', 1.2);
exportgraphics(f, fullfile(anaDir, "r_boxplot_t-resp.png"));

