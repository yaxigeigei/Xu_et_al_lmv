%% Create a library of plots for inspecting results

anaDir = LMV.Data.GetAnalysisDir('coding', 'stRF');
srcTb = LMV.Data.FindSource([]);

%% 

load(fullfile(anaDir, "segmental_clusTb.mat"));

load(fullfile(LMV.Data.GetAnalysisDir, "data", "ce_m2_session-avg.mat"), "ce");
fMel = ce.userData.melMeta.F;

%% Get units with valid models

sets = LMV.TRF.featSets;
targets = ["stim", "prod"];

nTargets = numel(targets);
nSets = numel(sets);
mdlTbs = cell(nSets, nTargets);

for i = 1 : nSets
    for j = 1 : nTargets
        tb = clusTb;
        mn = sets(i)+"_"+targets(j);
        isMdl = ~cellfun(@isempty, tb.(mn));
        mdls = tb.(mn)(isMdl);
        
        if sets(i) == "strf"
            for k = 1 : numel(mdls)
                mdls{k}.feats = fMel;
            end
        end
        
        tb = tb(isMdl, 1:end-9);
        tb.mdls = mdls;
        mdlTbs{i,j} = tb;
    end
end

%% Define example units

uIds = cell(nSets, nTargets);

for i = 1 : nSets
    for j = 1 : nTargets
        switch sets(i)+"_"+targets(j)
            case "phone_stim"
                uIds{i,j} = 440300585;
            case "phone_prod"
                uIds{i,j} = 410100245;
            case "strf_stim"
                % uIds{i,j} = mdlTbs{i,j}.clusId(endsWith(mdlTbs{i,j}.region, "PrCG")); % 450200052
                uIds{i,j} = 450200052;
            case "artic3_prod"
                uIds{i,j} = 540100167;
        end
    end
end

%% Plot optimal weights of example units

f = MPlot.Figure(710); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";

for i = 1 : numel(sets)
    for j = 1 : numel(targets)
        tb = mdlTbs{i,j};
        uid = uIds{i,j};
        if ismember(tb.mdls{1}.name, ["strf_prod", "artic3_stim"])
            continue
        end
        
        for k = 1 : numel(uid)
            isMdl = tb.clusId==uid(k);
            M = tb.mdls{isMdl};
            m = M.mdls{M.r2Idx};
            m.name = M.name;
            m.feats = M.feats;
            m.resps = M.resps;
            % m.Beta = zscore(m.Beta);
            
            ax = nexttile;
            LMV.RF.PlotWeights(m, 'Orientation', 'vertical', 'SigColor', 'k', 'Parent', ax);
            
            switch sets(i)
                case "phone"
                    ax.YTickLabel = MPlot.StaggerLabels(MLing.ARPA2IPA(ax.YTickLabel), -5);
                case "artic3"
                    ax.YTickLabel = MPlot.StaggerLabels(NP.Artic.GetLabels(ax.YTickLabel), -8);
                case "strf"
                    ax.YDir = "normal";
                    ax.YTick = [1:20:numel(m.feats) numel(m.feats)];
                    ax.YTickLabel = round(m.feats(ax.YTick)/1e3, 1);
                    ax.YLabel.String = "Frequency (Hz)";
                otherwise
                    ax.YTickLabel = feats;
            end
            MPlot.Axes(ax);
        end
    end
end

MPlot.Paperize(f, 1, 0.55);
figPath = fullfile(anaDir, "example_unit_weights.png");
exportgraphics(f, figPath);
print(f, strrep(figPath, ".png", ".pdf"), '-dpdf');

%% Plot heatmaps of optiaml RFs

% Collect tables
tb2plot = {};
u2mark = {};
for i = 1 : numel(sets)
    for j = 1 : numel(targets)
        tb = mdlTbs{i,j};
        if ~ismember(tb.mdls{1}.name, ["strf_prod", "artic3_stim"])
            tb2plot{end+1} = tb;
            u2mark{end+1} = uIds{i,j};
        end
    end
end

% Plot
f = MPlot.Figure(810); clf
% f.WindowState = "maximized";

nUnits = cellfun(@height, tb2plot)+50;
tl = tiledlayout(1, sum(nUnits));
tl.Padding = "compact";

for i = 1 : numel(tb2plot)
    ntArgs = MPlot.FindTileInd(1, nUnits, 1, i);
    ax = nexttile(ntArgs{:});
    LMV.RF.PlotPopWeights(ax, tb2plot{i}, 'SortBy', 'peak', 'MarkUnit', u2mark{i}(1));
end

MPlot.Paperize(f, 2.2, 0.7);
figPath = fullfile(anaDir, "pop_weight_heatmaps.png");
exportgraphics(f, figPath);
print(f, strrep(figPath, ".png", ".pdf"), '-dpdf');

