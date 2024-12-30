%% Spatial distributions of linker units

mdlName = "smooth_lm";
anaDir = LMV.Data.GetAnalysisDir("linker", mdlName);
srcTb = LMV.Data.FindSource([]);

%% Load data

clusTb = LMV.Linker.LoadClusTb(mdlName);
rTest = LMV.Resp.LoadPhaseResponseTest();

% Find units responsive to both stim and prod
isResp = all(rTest.sigTb{:, ["stim", "prod"]}, 2);

% Find manually identified units
recIds = srcTb.recId;
cid = cellfun(@(x) LMV.Linker.GetSelectedClusId(x, recIds), LMV.Linker.types, 'Uni', false);
isSelect = ismember(rTest.clusTb.clusId, [cid{:}]);

% Include responsive and the identified units
rTb = rTest.clusTb(isResp | isSelect, :);

%% For each linker type, plot the number of units in each recording

types = LMV.Linker.types([2 1 3]);
[recIds, ia] = unique(rTest.clusTb.recId, 'stable');

nr = numel(recIds);
nt = numel(types);
nLinkers = zeros(nr, nt);
for i = 1 : nr
    for j = 1 : nt
        nLinkers(i,j) = sum(clusTb.recId==recIds(i) & clusTb.hcGroup==types(j));
    end
end

f = MPlot.Figure(6); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";

for i = 1 : nt
    ax = nexttile;
    p = nLinkers(:,i);
    m = p > 0;
    h = donutchart(p(m));
    h.CenterLabel = [types(i) "kernels"];
    h.CenterLabelFontSize = 12;
    h.Names = replace(recIds(m), '_', '-');
    h.LabelStyle = "namedata";
    h.ColorOrder = LMV.Param.GetRegionColors(rTest.clusTb.region(ia(m)));
end

MPlot.Paperize(f, 1.6, .3);
exportgraphics(f, fullfile(anaDir, "linker_n_units.png"));
exportgraphics(f, fullfile(anaDir, "linker_n_units.pdf"), ContentType="vector");

%% For each linker type, plot the percentage of units from each region (normalized)

types = LMV.Linker.types([2 1 3]);
regions = LMV.Param.regions;

nr = numel(regions);
nt = numel(types);
nLinkers = zeros(nr, nt);
nUnit = zeros(nr, 1);
nRec = zeros(nr, 1);

for i = 1 : nr
    for j = 1 : nt
        nLinkers(i,j) = sum(clusTb.region==regions(i) & clusTb.hcGroup==types(j));
    end
    isRegion = rTest.clusTb.region==regions(i);
    nRec(i) = numel(unique(rTest.clusTb.recId(isRegion)));
    nUnit(i) = sum(isRegion);
end

nRec(regions=="IFG") = nRec(regions=="IFG") + 2; % for two double probes
nRec(regions=="mPrCG") = nRec(regions=="mPrCG") + 1; % for one double probe

f = MPlot.Figure(7); clf
tl = tiledlayout(2,3);
tl.Padding = "compact";

% Normalize by the total number of units recorded in each region
pLinkers = nLinkers./nUnit;
pLinkers = pLinkers./sum(pLinkers);

for i = 1 : nt
    ax = nexttile;
    p = pLinkers(:,i);
    m = p > 0;
    h = donutchart(p(m));
    h.CenterLabel = [types(i) "kernels"];
    h.CenterLabelFontSize = 12;
    h.Names = regions(m);
    h.LabelStyle = "namepercent";
    h.ColorOrder = LMV.Param.GetRegionColors(regions(m));
end

% Normalize by the number of sites in each region
pLinkers = nLinkers./nRec;
pLinkers = pLinkers./sum(pLinkers);

for i = 1 : nt
    ax = nexttile;
    p = pLinkers(:,i);
    m = p > 0;
    h = donutchart(p(m));
    h.CenterLabel = [types(i) "kernels"];
    h.CenterLabelFontSize = 12;
    h.Names = regions(m);
    h.LabelStyle = "namepercent";
    h.ColorOrder = LMV.Param.GetRegionColors(regions(m));
end

MPlot.Paperize(f, 1.6, .6);
exportgraphics(f, fullfile(anaDir, "linker_region_percentages.png"));
exportgraphics(f, fullfile(anaDir, "linker_region_percentages.pdf"), ContentType="vector");

%% For each region, plot the percentage of units in each linker type

types = [LMV.Linker.types([2 1 3]) "other"];
regions = LMV.Param.regions;

cc = LMV.Param.GetLinkerColors(types);

nr = numel(regions);
nt = numel(types);
nLinkers = zeros(nr, nt);
for i = 1 : nr
    for j = 1 : nt
        nLinkers(i,j) = sum(clusTb.region==regions(i) & clusTb.hcGroup==types(j));
    end
end

f = MPlot.Figure(6); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";

for i = 1 : nr
    ax = nexttile;
    p = nLinkers(i,:);
    m = p > 0;
    h = piechart(p(m));
    h.Title = regions(i);
    h.Names = types(m);
    h.LabelStyle = "namepercent";
    h.ColorOrder = LMV.Param.GetLinkerColors(types(m));
end

MPlot.Paperize(f, 1.6, .3);
exportgraphics(f, fullfile(anaDir, "linker_comp_by_region.png"));
exportgraphics(f, fullfile(anaDir, "linker_comp_by_region.pdf"), ContentType="vector");

%% Plot raw unit depths

types = LMV.Linker.types([2 1 3]);
regions = LMV.Param.regions;

cTb = clusTb;
cTb.region = categorical(cTb.region, regions, 'Ordinal', true);
cTb = sortrows(cTb, {'region', 'recId'});

f = MPlot.Figure(16); clf
tl = tiledlayout(1, numel(types));
tl.Padding = "compact";

for i = 1 : numel(types)
    isType = cTb.hcGroup==types(i);
    recIds = unique(cTb.recId(isType), 'stable');
    
    ax = nexttile;
    for j = 1 : numel(recIds)
        isUnit = isType & cTb.recId==recIds(j);
        region = cTb.region(find(isUnit,1));
        x = repelem(j, sum(isUnit));
        plot(x, cTb.depth(isUnit), 'o', 'Color', LMV.Param.GetRegionColors(region)); hold on
    end
    ax.YDir = "reverse";
    ax.XLim = [0 j+1];
    ax.XTick = 1:j;
    ax.XTickLabel = recIds;
    ax.TickLabelInterpreter = "none";
    ax.Title.String = types(i);
    MPlot.Axes(ax);
end

MPlot.Paperize(f, 1.5, 0.6);
exportgraphics(f, fullfile(anaDir, "unit_depth_by_sites.png"));

%% Fit depth distributions

cTb = clusTb;
cTb.region = categorical(cTb.region, LMV.Param.regions, 'Ordinal', true);
cTb = sortrows(cTb, {'region', 'recId'});

types = LMV.Linker.types([2 1 3]);
regions = ["mPrCG", "STG"];

pdLink = cell(numel(types), numel(regions));
pdResp = cell(numel(types), numel(regions));
pval = NaN(size(pdLink));

x = linspace(0, 5000, 100);
dArgs = {'Kernel', 'Kernel', 'normal', 'Bandwidth', 500};

for i = 1 : numel(types)
    isType = cTb.hcGroup==types(i);
    % regions = string(unique(cTb.region(isType)));
    
    for j = 1 : numel(regions)
        % Linker distribution
        isLink = isType & cTb.region==regions(j);
        if sum(isLink) > 5
            d1 = cTb.depth(isLink);
            pd = fitdist(d1, dArgs{:});
            pdLink{i,j} = pdf(pd, x);
        else
            pdLink{i,j} = x*NaN;
        end
        
        % All responsive distribution
        isResp = rTest.clusTb.region==regions(j) & any(rTest.sigTb{:,:},2);
        d2 = rTest.clusTb.depth(isResp);
        pd = fitdist(d2, dArgs{:});
        pdResp{i,j} = pdf(pd, x);
        
        if sum(isLink) > 5
            [~, pval(i,j)] = kstest2(d1, d2);
        end
    end
end

pval

%% 

f = MPlot.Figure(17); clf
tl = tiledlayout(2, numel(types));
tl.Padding = "compact";

for j = 1 : numel(regions)
    for i = 1 : numel(types)
        ax = nexttile;
        
        yLink = pdLink{i,j};
        yLink = yLink/sum(yLink);
        cc = LMV.Param.GetRegionColors(regions(j));
        plot(yLink, x/1e3, '-', 'Color', cc); hold on
        
        yResp = pdResp{i,j};
        yResp = yResp/sum(yResp);
        plot(yResp, x/1e3, ':', 'Color', cc, 'LineWidth', 2);
        
        ax.YDir = "reverse";
        ax.XLim = [0 .035];
        ax.TickLabelInterpreter = "none";
        ax.Title.String = types(i);
        MPlot.Axes(ax);
    end
end

MPlot.Paperize(f, .8, 1);
exportgraphics(f, fullfile(anaDir, "depth_dist_abs.png"));

%% 

f = MPlot.Figure(18); clf
tl = tiledlayout(1, numel(regions));
tl.Padding = "compact";

cc = zeros(3);
ls = ["-", "--", ":"];

for j = 1 : numel(regions)
    ax = nexttile;
    
    for i = 1 : numel(types)
        yLink = pdLink{i,j};
        yLink = yLink/sum(yLink);
        
        yResp = pdResp{i,j};
        yResp = yResp/sum(yResp);
        
        plot(yLink./yResp, x/1e3, Color=cc(i,:), LineStyle=ls(i), LineWidth=1.5); hold on
    end
    
    plot([1 1]', x([1 end])'/1e3, 'Color', [0 0 0 .3]);
    
    ax.XLabel.String = "Enrichment (fold)";
    ax.YLabel.String = "Distance from surface (mm)";
    ax.YDir = "reverse";
    ax.XLim = [0 3];
    ax.TickLabelInterpreter = "none";
    ax.Title.String = regions(j);
    MPlot.Axes(ax);
end

legend(ax.Children(end:-1:2), types, Location="eastoutside");

MPlot.Paperize(f, .7, .5);
exportgraphics(f, fullfile(anaDir, "depth_dist_ratio.png"));
exportgraphics(f, fullfile(anaDir, "depth_dist_ratio.pdf"));

