%% 

anaDir = LMV.Data.GetAnalysisDir('units');
srcTb = LMV.Data.FindSource([]);

%% Load QC tables

qcPaths = fullfile(anaDir, 'computed_qc', srcTb.recId+"_clusTb.mat");
uTbs = cell(size(qcPaths));
for i = 1 : numel(qcPaths)
    load(qcPaths(i), 'clusTb');
    uTbs{i} = clusTb;
end

%% Load se and compute log ISI histograms

tEdges = cumsum([repelem(1e-3, 20) repelem(2e-3, 20) repelem(4e-3, 20) repelem(8e-3, 20)]);

for i = 1 : height(srcTb)
    % Check cache
    cacheDir = fullfile(anaDir, "computed_log_isi");
    cachePath = fullfile(cacheDir, srcTb.recId{i}+"_clusTb_log-isi.mat");
    if exist(cachePath, 'file')
        fprintf("Cached log ISI histograms is found for %s.\n", srcTb.recId{i});
        load(cachePath, "clusTb");
        uTbs{i} = clusTb;
        continue
    end
    
    % Load se
    se = NP.SE.LoadSession(srcTb.path{i});
    se = se.Duplicate({'spikeTime'});
    se.SliceSession(0, 'absolute');
    
    % Compute log ISI histogram
    st = se.GetTable('spikeTime');
    clusTb = uTbs{i};
    for u = 1 : width(st)
        % Compute inter-spike intervals
        tSpk = unique(st.(u){1});
        isi = diff(tSpk);
        
        % ISI histogram
        nISI = histcounts(isi, [tEdges Inf]);
        nISI = nISI(1:end-1);
        clusTb.isiEdges{u} = tEdges;
        clusTb.isiCount{u} = nISI;
    end
    uTbs{i} = clusTb;
    
    % Cache table
    if ~exist(cacheDir, 'dir')
        mkdir(cacheDir);
    end
    save(cachePath, "clusTb");
end
uTb = cat(1, uTbs{:});

% Load waveform results
sWf = load(fullfile(anaDir, "computed_phys_clusTb.mat"), "clusTb");
% subtypes = categorical(["Pyr", "PosN", "PosW", "FS"], ["FS", "Pyr", "PosN", "PosW"], 'Ordinal', true)';
sWf.clusTb.wfId = subtypes(sWf.clusTb.wfId);
uTb = [uTb sWf.clusTb(:,["wfId" "t2p"])];

%% UMAP embedding of unit waveform

H = cat(1, uTb.isiCount{:});

rng(61);
Z = run_umap(H, 'metric', 'correlation', 'verbose', 'text', 'randomize', true, 'n_neighbors', 20);

%% Labeling the main clusters with K-means on the embedding

rng(61);
k = 4;
isiKmId = kmeans(Z, k);

subtypes = categorical(["Long2", "Burst2", "Long1", "Burst1"], ["Burst1", "Burst2", "Long1", "Long2"], 'Ordinal', true)';
uTb.isiKmId = subtypes(isiKmId);

%% Interactive plot of embedding and clustering

% With ISIH embedding
uTb.embedCoords = Z;
gVars = {'isiKmId', 'depth', 'wfId', 'region'};

for i = 1 : numel(gVars)
    f = MPlot.Figure(9700+i); clf
    NP.UnitPlot.WaveformEmbedding(uTb, gVars{i}, 'logISI');
    MPlot.Paperize(f, 'ColumnsWide', 0.8, 'Aspect', .8);
    exportgraphics(f, fullfile(anaDir, "umap_isi_by_"+gVars{i}+".png"));
end

% With waveform embedding
uTb.embedCoords = sWf.clusTb.embedCoords .* [ones(height(uTb),1) ones(height(uTb),1)];
gVars = {'isiKmId'};

for i = 1 : numel(gVars)
    f = MPlot.Figure(9800+i); clf
    NP.UnitPlot.WaveformEmbedding(uTb, gVars{i}, 'logISI');
    MPlot.Paperize(f, 'ColumnsWide', 0.8, 'Aspect', .8);
    exportgraphics(f, fullfile(anaDir, "umap_waveform_by_"+gVars{i}+".png"));
end

uTb.embedCoords = Z;

%% Plot average ISIH of each cluster

f = MPlot.Figure(9800); clf
NP.UnitPlot.MeanIsiArray(uTb, 'isiKmId');
MPlot.Paperize(f, 'ColumnsWide', .8, 'ColumnsHigh', .6);
exportgraphics(f, fullfile(anaDir, "kmeans_centroid_isi.png"));

%% Physiological-anatomical distributions

xVars = {'isiKmId', 'isiKmId'};
xNames = xVars;

yVars = {'wfId', 'depth'};
yNames = yVars;

regions = [unique(uTb.region); "all"];

for i = 1 : numel(xVars)
    f = MPlot.Figure(3251); clf
    tl = tiledlayout('flow');
    tl.Padding = 'loose';
    for r = 1 : numel(regions)
        % Select units in the region of interest
        isRegion = strcmp(regions(r), uTb.region);
        if ~any(isRegion)
            isRegion = true(size(uTb.region));
        end
        uTbSub = uTb(isRegion,:);
        
        % Plot distribution
        ax = nexttile;
        s = LMV.FPA.ComputeDist(uTbSub, xVars{i}, yVars{i}, 'nBoot', 1e3);
        LMV.FPA.PlotDist(s, 'Style', "heatmap");
        ax.CLim(2) = prctile(ax.Children(end).CData(:), 99);
        ax.Title.String = sprintf("%s (n = %i)", regions(r), sum(isRegion));
    end
    MPlot.Paperize(f, 'ColumnsWide', 1.1, 'ColumnsHigh', .9);
    exportgraphics(f, fullfile(anaDir, "p_"+xNames{i}+"_given_"+yVars{i}+".png"));
end

%% Save unit tables with physiological info

clusTb = uTb;
save(fullfile(anaDir, "computed_phys_clusTb.mat"), 'clusTb');

return
%% Sanity check if clusters were chosen correctly

m = find(t2p > 0.5);
if numel(m) > 50
    m = randsample(m, 50);
end

f = MPlot.Figure(214040); clf
panels = cell(10, 5);
panelArgs = panels;
for i = 1 : numel(m)
    panels{i} = 'waveform';
    panelArgs{i} = {i};
end
NP.UnitPlot.Array(uTb(m,:), panels', 'PanelArgs', panelArgs');


