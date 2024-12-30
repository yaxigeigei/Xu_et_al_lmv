%% 

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'trf');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end

%% Load computed data

% Load TRF models
setName = "combo3pm";
mdlNames = setName + "_" + NP.TRF.targets;
mdlDir = fullfile(anaDir, setName, "mdls");
mdlTbs = NP.TRF.LoadModels(mdlNames, mdlDir);
mdlTb = cat(1, mdlTbs{:});

% Load physiological information of units
sPhys = load(fullfile(NP.Data.GetAnalysisRoot, 'units', 'computed_isi_clusTb.mat'));
clusIdPhys = NP.Unit.ToUniqueClusId(sPhys.clusTb.clusId, sPhys.clusTb.recId);
[clusIdPhys, I] = MMath.SortLike(clusIdPhys, mdlTb.clusId);
physVars = {'wfId', 't2p', 'isiKmId'};
mdlTb(:,physVars) = sPhys.clusTb(I,physVars);

% Load TRF clustering results
load(fullfile(anaDir, setName, "computed_hc.mat"), 'sHCxAll');

%% Add data to clusTb

% Add interpreted cluster names
s = sHCxAll;
[s, d] = NP.TRF.NameClusters(s, setName+"-1");
clusNames = fieldnames(d);

% Sort models by optimal leaf order for consistency with clustering heatmap plots
leafOrder = optimalleaforder(s.Z, s.D, 'Criteria', 'adjacent');
fns = fieldnames(s);
for i = 1 : numel(fns)
    fn = fns{i};
    if size(s.(fn), 1) == numel(leafOrder)
        fprintf("Reorder '%s'\n", fn)
        s.(fn) = s.(fn)(leafOrder,:);
    end
end

% Match clusTb rows to s.resps
[~, I] = MMath.SortLike("u"+mdlTb.clusId, s.resps);
mdlTb2 = mdlTb(I,:);
mdlTb2.hcId = s.hcId;
mdlTb2.hcName = s.hcName;
mdlTb2.embedCoords = s.umapCoords;
mdlTb2.mdlName = s.mdlNames;
for i = 1 : numel(s.mdlNames)
    mdlTb2.mdl{i} = mdlTb2.(s.mdlNames(i)){i};
end
mdlTb2.maxBeta = s.maxBeta;

%% Compute average TRFs

regions = ["mPrCG", "STG"];
avgTbs = cell(numel(regions), numel(mdlNames));
for i = 1 : numel(regions)
    isRegion = mdlTb2.region == regions(i);
    for j = 1 : numel(mdlNames)
        isMdl = mdlTb2.mdlName == mdlNames(j);
        subTb = mdlTb2(isRegion & isMdl,:);
        avgTbs{i,j} = NP.TRF.AverageModels(subTb, mdlNames(j), 'hcName', clusNames);
        avgTbs{i,1}.(mdlNames(j)) = avgTbs{i,j}.(mdlNames(j)); % put all models in the first table
    end
end

%% Plot cluster mean TRFs

for i = 1 : numel(regions)
    f = MPlot.Figure(102000+i); clf
    NP.TRF.PlotWeightsCombo(avgTbs{i,1}, mdlNames, 'UnitsPerPage', s.nc, 'Page', 1);
    exportgraphics(f, fullfile(anaDir, setName, sprintf('avg_trf_weights_%s.png', regions(i))));
end

%% Plot FPA on embedding

regions = ["mPrCG", "STG"];
colorByList = ["hcName", "recId", "depth", "waveform"];

for i = 1 : numel(regions)
    f = MPlot.Figure(4300+i); clf
    NP.TRF.PlotEmbeddingFPA(s, mdlTb2, colorByList, regions(i), mdlNames);
    MPlot.Paperize(f, 'ColumnsWide', 1.2, 'ColumnsHigh', 0.8);
    exportgraphics(f, fullfile(anaDir, setName, sprintf("embed_%s_by_%s.png", regions(i), strjoin(colorByList, '_'))));
end

%% Compute spatial similarity

recs = unique(mdlTb2.recId);
simRec = cell(numel(recs), numel(mdlNames));
for i = 1 : numel(recs)
    isRec = mdlTb2.recId == recs(i);
    for j = 1 : numel(mdlNames)
        isMdl = mdlTb2.mdlName == mdlNames(j);
        subTb = mdlTb2(isRec & isMdl,:);
        simRec{i,j} = NP.TRF.ComputePairwiseCorr(subTb, 0);
    end
end

%% 

regions = ["mPrCG", "STG"];
recRegions = cellfun(@(x) x.region, simRec);
simReg = cell(numel(regions), numel(mdlNames));
for i = 1 : numel(regions)
    isRegion = mdlTb2.region == regions(i);
    for j = 1 : numel(mdlNames)
        % Combine pairwise data
        m = recRegions(:,j) == regions(i);
        ss = simRec(m,j);
        sCat = ss{1};
        sCat.r = cell2mat(cellfun(@(x) x.r, ss, 'Uni', false));
        sCat.d = cell2mat(cellfun(@(x) x.d, ss, 'Uni', false));
        
        % Compute similarity
        simReg{i,j} = NP.TRF.ComputeSpatialSimilarity(sCat, 100);
    end
end

%% Plot spatial similarity

recRegions = cellfun(@(x) x.region, simRec);

f = MPlot.Figure(4400); clf
tl = tiledlayout(numel(regions), numel(mdlNames));
tl.Padding = 'compact';
for i = 1 : numel(regions)
    for j = 1 : numel(mdlNames)
        ax = nexttile;
        m = recRegions(:,j) == regions(i);
        NP.TRF.PlotSpatialSimilarity(simReg{i,j}, simRec(m,j));
    end
end
MPlot.Paperize(f, 'ColumnsWide', 2, 'ColumnsHigh', 0.8);
exportgraphics(f, fullfile(anaDir, setName, "trf_spatial_corr.png"));


