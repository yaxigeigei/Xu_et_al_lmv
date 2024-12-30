%% Spatial distributions of linker units

mdlName = "smooth_lm";
anaDir = LMV.Data.GetAnalysisDir("linker", mdlName);

%% Load data

% Load linker info
clusTb = LMV.Linker.LoadClusTb(mdlName);

% Load waveform info
sPhys = load(fullfile(LMV.Data.GetAnalysisDir, "units", "computed_phys_clusTb.mat"));
[~, I] = MMath.SortLike(sPhys.clusTb.clusId, clusTb.clusId);
clusTb.wfId = sPhys.clusTb.wfId(I);

%% Plot waveforms by depth

types = [LMV.Linker.types([2 1 3]), "other"];
regions = ["mPrCG", "STG"];

cTb = clusTb;
cTb.region = categorical(cTb.region, regions, 'Ordinal', true);
cTb = sortrows(cTb, {'region', 'recId'});

f = MPlot.Figure(16); clf
tl = tiledlayout(1, numel(regions));
tl.Padding = "compact";
for r = 1 : numel(regions)
    ax = nexttile;
    k = 0;
    xLb = string([]);
    for t = 1 : numel(types)
        isType = cTb.hcGroup==types(t);
        isRegion = cTb.region==regions(r);
        isRange = cTb.depth < 5e3;
        isUnit = isType & isRegion & isRange;
        if sum(isUnit) < 3
            continue
        end
        k = k + 1;
        xLb(end+1) = types(t);
        NP.UnitPlot.WaveformAtDepth(ax, cTb(isUnit,:), k, 'TimeScaling', 0.1, 'WaveformScaling', .2);
    end
    ax.YLim = [0 5];
    ax.XLim = [0.4 k+.6];
    ax.XTick = 1:k;
    ax.XTickLabel = xLb;
    ax.YLabel.String = "Distance from surface (mm)";
    ax.TickLabelInterpreter = "none";
    ax.Title.String = regions(r);
    MPlot.Axes(ax);
end

MPlot.Paperize(f, 0.8, 0.5);
exportgraphics(f, fullfile(anaDir, "unit_waveform_at_depth.png"));
exportgraphics(f, fullfile(anaDir, "unit_waveform_at_depth.pdf"));
