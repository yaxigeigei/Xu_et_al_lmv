%% Create a library of plots for inspecting results

mdlName = "smooth_lm";
anaDir = LMV.Data.GetAnalysisDir("linker", mdlName);
libDir = LMV.Data.GetAnalysisDir("linker", mdlName, "lib_tuning");

srcTb = LMV.Data.FindSource([]);

%% Find units

% Load linker clustering result
clusTb = LMV.Linker.LoadClusTb(mdlName);

% Get clusIds of each linker group
types = LMV.Linker.types;
hcGroups = arrayfun(@(x) clusTb.clusId(clusTb.hcGroup==x)', types, 'Uni', false);

% Load unit cache
uCache = cellfun(@NP.Unit.LoadUnitCache, hcGroups, 'Uni', false);

%% Percentage of time-locked tuning

sets = ["strf", "artic3", "phone"];
targets = LMV.TRF.targets;
mdlNames = sets + "_" + targets';

r2frac = NaN(numel(mdlNames), numel(types));
for t = 1 : numel(types)
    for m = 1 : numel(mdlNames)
        mn = mdlNames(m);
        TRFs = cellfun(@(x) x.trf.(mn), uCache{t}, 'Uni', false)';
        TRFs = LMV.RF.UpdateModelR2(TRFs);
        r2 = NaN(size(TRFs));
        for u = 1 : numel(TRFs)
            if ~isempty(TRFs{u})
                r2(u) = TRFs{u}.r2;
            end
        end
        r2frac(m,t) = mean(r2>0);
    end
end
r2frac = reshape(r2frac, [size(mdlNames), numel(types)]);

f = MPlot.Figure(779); clf
tl = tiledlayout("vertical");
tl.Padding = "compact";

cc = lines(numel(sets));
mdlInd = [1 3; 1 3; 2 3];

for t = 1 : numel(targets)
    ax = nexttile;
    ind = mdlInd(t,:);
    v = squeeze(r2frac(t,ind,:))';
    bb = bar(v);
    for i = 1 : numel(bb)
        bb(i).FaceColor = cc(ind(i),:);
    end
    ax.XTickLabel = types;
    ax.YLabel.String = "Frac.";
    ax.Title.String = targets(t);
    MPlot.Axes(ax);
    lgd = legend(bb, sets(ind), Location="eastoutside");
end

MPlot.Paperize(f, 0.8, 0.7);
exportgraphics(f, fullfile(anaDir, "frac_time-locked.png"));
exportgraphics(f, fullfile(anaDir, "frac_time-locked.pdf"));

%% Extract models with matched tuning time between stim and prod

sets = "phone";
targets = ["stim", "prod"];
mdlNames = sets + "_" + targets;

mTRF = cell(size(types));
mRF = mTRF;
for l = 1 : numel(types)
    % Extract models
    nUnits = numel(uCache{l});
    TRFs = cell(nUnits, 1);
    for u = 1 : nUnits
        TRFs{u} = cellfun(@(x) uCache{l}{u}.trf.(x), mdlNames, 'Uni', false);
    end
    TRFs = cat(1, TRFs{:});
    isRm = any(cellfun(@isempty, TRFs), 2);
    TRFs(isRm,:) = [];
    
    % Check time-locking
    TRFs = LMV.RF.UpdateModelR2(TRFs);
    r2 = NaN(size(TRFs));
    for m = 1 : numel(TRFs)
        if ~isempty(TRFs{m})
            r2(m) = TRFs{m}.r2;
        end
    end
    TRFs(all(isnan(r2),2),:) = [];
    
    % Match the tuning time between stim and prod
    for u = 1 : size(TRFs,1)
        TRFs(u,:) = LMV.Linker.MatchRF(TRFs(u,:));
    end
    mTRF{l} = TRFs;
    
    RFs = TRFs;
    for u = 1 : size(TRFs,1)
        for p = 1 : size(TRFs,2)
            TRF = TRFs{u,p};
            RF = TRF.mdls{TRF.r2Idx};
            RF.feats = TRF.feats;
            RF.resps = TRF.resps;
            RFs{u,p} = RF;
        end
    end
    mRF{l} = RFs;
end

%% Visualize phoneme RF weights

nTargets = numel(targets);
nMaxUnit = 15;

for l = 1 : numel(types)
    % Initialze figure
    f = MPlot.Figure(780); clf
    f.WindowState = "maximized";
    tl = tiledlayout(nTargets, nMaxUnit);
    tl.Padding = "tight";
    
    TRFs = mTRF{l};
    for u = 1 : min(size(TRFs,1), nMaxUnit)
        for p = 1 : size(TRFs,2)
            ntArgs = MPlot.FindTileInd(nTargets, nMaxUnit, p, u);
            ax = nexttile(ntArgs{:});
            LMV.TRF.PlotWeights2(TRFs{u,p});
            if targets(p) == "prod"
                ax.XDir = "reverse";
            end
        end
    end
    exportgraphics(f, fullfile(libDir, sprintf("%s_matched_stRF.png", types(l))))
    % break
end

for  l = 1 : numel(types)
    % Initialze figure
    f = MPlot.Figure(790); clf
    f.WindowState = "maximized";
    tl = tiledlayout(nTargets, nMaxUnit);
    tl.Padding = "tight";
    
    RFs = mRF{l};
    for u = 1 : min(size(RFs,1), nMaxUnit)
        for p = 1 : size(RFs,2)
            ntArgs = MPlot.FindTileInd(nTargets, nMaxUnit, p, u);
            ax = nexttile(ntArgs{:});
            LMV.RF.PlotWeights(RFs{u,p}, 'Orientation', 'vertical');
        end
    end
    exportgraphics(f, fullfile(libDir, sprintf("%s_matched_RF.png", types(l))))
    % break
end

%% Plot heatmaps of optiaml RFs

% Initialze figure
f = MPlot.Figure(810); clf
% f.WindowState = "maximized";

nTargets = numel(targets);
nTypes = numel(types);
tl = tiledlayout(nTargets, nTypes);
tl.Padding = "compact";

for k = 1 : numel(types)
    % Construct table for plot
    mdls = mTRF{k};
    cid = cellfun(@(x) x.resps, mdls(:,1));
    cid = str2double(erase(cid, "u"));
    [~, I] = MMath.SortLike(clusTb.clusId, cid);
    stimTb = clusTb(I,:);
    stimTb.mdls = mdls(:,1);
    
    % Plot weights as a heatmap
    ntArgs = MPlot.FindTileInd(nTargets, nTypes, 1, k);
    ax = nexttile(ntArgs{:});
    uInd = LMV.RF.PlotPopWeights(ax, stimTb, 'SortBy', 'hc');
    ax.Title.String = types(k) + " " + targets(1);
    
    prodTb = stimTb;
    prodTb.mdls = mdls(:,2);
    prodTb = prodTb(uInd,:);
    
    ntArgs = MPlot.FindTileInd(nTargets, nTypes, 2, k);
    ax = nexttile(ntArgs{:});
    LMV.RF.PlotPopWeights(ax, prodTb, 'SortBy', 'none');
    ax.Title.String = types(k) + " " + targets(2);
end

MPlot.Paperize(f, 2, 1.5);
% figPath = fullfile(libDir, "matched_tuning_heatmaps.png");
% exportgraphics(f, figPath);
% print(f, strrep(figPath, ".png", ".pdf"), '-dpdf');
