%% Embedding of responsive units

anaDir = LMV.Data.GetAnalysisDir('embedding', 'umap');
figDir = LMV.Data.GetAnalysisDir('embedding', 'umap', 'resp_units');

%% Load data

% Load grand average (M2) PETHs
load(fullfile(NP.Data.GetAnalysisRoot, 'data', 'ce_m2_grand-avg.mat'), 'ce');

% Load responsiveness
rTest = LMV.Resp.LoadPhaseResponseTest();
tTb = rTest.clusTb;

% Add tId to clusTb in ce
[~, I] = MMath.SortLike(tTb.clusId, ce.clusTb.clusId);
ce.clusTb.tId1 = tTb.tId1(I);

% Remove non-responsive untis
ce.RemoveUnits(ce.clusTb.tId1=="none");

% Split ce by regions
regions = ["all", "mPrCG", "vPrCG", "IFG", "STG"];
ceArray = cell(size(regions));
for r = 1 : numel(regions)
    ceRegional = ce.Duplicate;
    if regions(r) ~= "all"
        ceRegional.RemoveUnits(ce.clusTb.region ~= regions(r));
    end
    ceArray{r} = ceRegional;
end
ceArray = cat(1, ceArray{:});

%% Compute UMAP embedding

for r = 1 : numel(regions)
    ce = ceArray(r);
    NP.Embed.BatchUMAP(ce);
    cachePath = fullfile(anaDir, "computed_umap_"+regions(r)+".mat");
    save(cachePath, 'ce');
end

