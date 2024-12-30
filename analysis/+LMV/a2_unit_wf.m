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
uTb = cat(1, uTbs{:});

%% UMAP embedding of unit waveform

W = cellfun(@(x) x', uTb.waveformMed, 'Uni', false);
W = cellfun(@(x) x(:), W, 'Uni', false);
W = cat(2, W{:});

rng(61);
Z = run_umap(W', 'metric', 'correlation', 'verbose', 'text', 'randomize', false, 'n_neighbors', 10);
uTb.embedCoords = Z;

%% Find and label the clusters of waveform embedding

rng(61);
k = 4;
gmm = fitgmdist(Z, k);
wfId = cluster(gmm, Z);

subtypes = categorical(["Pyr", "PosN", "FS", "PosW"], ["FS", "Pyr", "PosN", "PosW"], 'Ordinal', true)';
uTb.wfId = subtypes(wfId);

%% Interactive plot of embedding and clustering

gVars = {'wfId', 'depth', 'region'};

for i = 1 : numel(gVars)
    f = MPlot.Figure(9700+i); clf
    NP.UnitPlot.WaveformEmbedding(uTb, gVars{i});
    MPlot.Paperize(f, 'ColumnsWide', .4, 'Aspect', .7, 'FontSize', 5, 'Zoom', 3);
    exportgraphics(f, fullfile(anaDir, "umap_waveform_by_"+gVars{i}+".png"));
    exportgraphics(f, fullfile(anaDir, "umap_waveform_by_"+gVars{i}+".pdf"));
end

%% Plot average waveform of each cluster

f = MPlot.Figure(9800); clf
NP.UnitPlot.MeanWaveformArray(uTb, 'wfId');
MPlot.Paperize(f, 'ColumnsWide', .4, 'Aspect', 1);
exportgraphics(f, fullfile(anaDir, "centroid_waveforms.png"));
exportgraphics(f, fullfile(anaDir, "centroid_waveforms.pdf"));

%% Report numbers of units by waveform type

fileID = fopen(fullfile(anaDir, 'waveform_stats.txt'), 'w');
for i = 1 : k
    fprintf(fileID, '%s\t%i\n', string(subtypes(i)), sum(uTb.wfId==subtypes(i)));
end
fclose(fileID);

%% Plot 

[~, ~, depthBin] = histcounts(uTb.depth, [0 2000 4000 6000]);
depthTypes = ["OOR", "0-2mm", "2-4mm", "4-6mm"];
uTb.depthType = depthTypes(depthBin+1)';

pyrTb = uTb(uTb.wfId=="Pyr",:);
pyrTb(pyrTb.depthType=="OOR",:) = [];

f = MPlot.Figure(9810); clf
NP.UnitPlot.MeanWaveformArray(pyrTb, 'depthType');
MPlot.Paperize(f, 'ColumnsWide', .4, 'Aspect', 1);
exportgraphics(f, fullfile(anaDir, "pyr_waveforms_by_depth.png"));
exportgraphics(f, fullfile(anaDir, "pyr_waveforms_by_depth.pdf"));

%% Plot peak-to-trough time of clusters 1 and 4

t2pTypes = subtypes([3 1])';
m = any(uTb.wfId == t2pTypes, 2);
t2p = NaN(size(m));

for i = 1 : height(uTb)
    % 
    if ~m(i)
        continue
    end
    
    W = uTb.waveformMed{i};
    [nChan, nTime] = size(W);
    c0 = round(nChan/2);
    W = W(c0,:);
    Wtr = W;
    
    [~, iTr] = min(Wtr);
    [~, iPk] = max(W(iTr:end));
    t2p(i) = iPk / 30;
end

uTb.t2p = t2p;

%% Distribution of trough-to-peak time 

cc = brewermap(k, 'Set1');
% t2pTypes = ["FS", "PyrU", "PyrL"];

f = MPlot.Figure(9820); clf
for i = 1 : numel(t2pTypes)
    if ismember(t2pTypes(i), subtypes)
        m = uTb.wfId == t2pTypes(i);
    elseif t2pTypes(i) == "PyrU"
        m = uTb.wfId == "Pyr" & uTb.depthType == "0-3mm";
    elseif t2pTypes(i) == "PyrL"
        m = uTb.wfId == "Pyr" & uTb.depthType == "3-6mm";
    end
    h = histogram(t2p(m), 0:0.05:1.5); hold on
    h.FaceColor = cc(i,:);
end
ax = MPlot.Axes(gca);
ax.XLabel.String = 'Trough-to-peak time (ms)';
ax.YLabel.String = '# of units';
lgd = legend(t2pTypes);

MPlot.Paperize(f, .5, .25);
exportgraphics(f, fullfile(anaDir, "trough-to-peak_histograms.png"));
exportgraphics(f, fullfile(anaDir, "trough-to-peak_histograms.pdf"));

%% Physiological-anatomical distributions

xVars = {'wfId'};
xNames = xVars;

yVars = {'depth'};
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


