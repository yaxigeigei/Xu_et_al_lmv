%% Hierarchical clustering of 

anaDir = LMV.Data.GetAnalysisDir('trf');
srcTb = LMV.Data.FindSource([]);

%% Load computed data

% Load sentence selectivity
sTest = LMV.Resp.LoadSentenceSelectTest();
testTb = sTest.clusTb;
sSig = sTest.sigTb{:,["stim" "prod" "prod"]};

% Load TRF models
setName = "combo3pm";
mdlNames = setName + "_" + NP.TRF.targets;
mdlDir = fullfile(anaDir, setName, "mdls");
mdlTbs = NP.TRF.LoadModels(mdlNames, mdlDir);
mdlTb = cat(1, mdlTbs{:});

% Make unit order in mdlTb consistent with testTb
[~, I] = MMath.SortLike(mdlTb.clusId, testTb.clusId);
mdlTb = mdlTb(I,:);

% Exclude non-selective or insignificant models
mdls = mdlTb{:,mdlNames};
isSig = cellfun(@(x) isfield(x, 'null') && x.null.r2Pval < 0.05, mdls);
indExclude = find(~sSig | ~isSig);
for i = indExclude(:)'
    mdls{i} = [];
end
mdlTb{:,mdlNames} = mdls;
isEpt = all(cellfun(@isempty, mdls), 2);
mdlTb(isEpt,:) = [];

%% Compute hierarchical clustering

nc = 5;
hcName = "hc" + nc;

regions = unique(mdlTb.region);
% regions = "mPrCG";

cachePath = fullfile(anaDir, setName, hcName, "computed_hc.mat");

if exist(cachePath, 'file')
    % Load previously computed
    load(cachePath);
else
    % Compute clustering
    sHC = cell(numel(regions), numel(mdlNames));
    sHCxPhase = cell(numel(regions), 1);
    sHCxReg = cell(1, numel(mdlNames));
    
    for i = 1 : numel(regions)
        % Select models
        isRegion = mdlTb.region == regions(i);
        regMdlTb = mdlTb(isRegion,:);
        mdls = regMdlTb{:,mdlNames};
        regMat = repmat(regMdlTb.region, [1 numel(mdlNames)]);
        
        % Cluster TRFs within each task phase
        for j = 1 : numel(mdlNames)
            s = NP.TRF.ClusterModels(mdls(:,j));
            s.regions = regMdlTb.region(s.hasMdl);
            sHC{i,j} = s;
        end
        
        % Cluster TRFs across task phases
        s = NP.TRF.ClusterModels(mdls);
        s.regions = regMat(s.hasMdl);
        sHCxPhase{i} = s;
    end
    
    % Cluster TRFs across regions but for each task phase
    mdls = mdlTb{:,mdlNames};
    for i = 1 : numel(mdlNames)
        s = NP.TRF.ClusterModels(mdls(:,i));
        s.regions = mdlTb.region(s.hasMdl);
        sHCxReg{i} = s;
    end
    
    % Cluster TRFs across all
    mdls = mdlTb{:,mdlNames};
    regMat = repmat(mdlTb.region, [1 numel(mdlNames)]);
    s = NP.TRF.ClusterModels(mdls);
    s.regions = regMat(s.hasMdl);
    sHCxAll = s;
    
    % Cache results
    save(cachePath, 'sHC', 'sHCxPhase', 'sHCxReg', 'sHCxAll');
end

%% Plot results of clustering with each region and task phase

% for i = 1 : numel(regions)
%     f = MPlot.Figure(4400+i); clf
%     f.Position(3:4) = [900 450];
%     NP.TRF.PlotClusterCombo(sHC(i,:));
%     MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', 0.9, 'FontSize', 5);
%     exportgraphics(f, fullfile(anaDir, setName, hcName, sprintf("clus_trf_each_%s.png", regions(i))));
% end

%% Plot results of clustering across all regions and task phases

for i = 1 : numel(regions)
    f = MPlot.Figure(4400+i); clf
    NP.TRF.PlotXAllClusterCombo(sHCxAll, regions(i));
%     MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', 0.9, 'FontSize', 5);
%     exportgraphics(f, fullfile(anaDir, setName, hcName, sprintf("clus_trf_xall_%s.png", regions(i))));
end

%% Plot results of clustering across task phases

for i = 1 : numel(regions)
    f = MPlot.Figure(4400+i); clf
    NP.TRF.PlotXModelClusterCombo(sHCxPhase{i});
%     MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', 0.9, 'FontSize', 5);
%     exportgraphics(f, fullfile(anaDir, setName, hcName, sprintf("clus_trf_xphase_%s.png", regions(i))));
end

%% Plot results of clustering across regions

for i = 1 : numel(regions)
    f = MPlot.Figure(4400+i); clf
    NP.TRF.PlotXRegionClusterCombo(sHCxReg, regions(i));
%     MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', 0.9, 'FontSize', 5);
%     exportgraphics(f, fullfile(anaDir, setName, hcName, sprintf("clus_trf_xregion_%s.png", regions(i))));
end

return
%% Plot the same units across periods

for i = 1 : numel(regions)
    for j = 1 : numel(mdlNames)
        f = MPlot.Figure(4400+j); clf
        NP.TRF.PlotClusterComboX(sHC(i,:,1), j);
%         MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', 0.9, 'FontSize', 5);
        exportgraphics(f, fullfile(anaDir, setName, sprintf("clus_trf_x_%s_by_%s.png", regions(i), mdlNames(j))));
    end
end

%% 

for i = 1 : numel(regions)
    for j = 1 : numel(mdlNames)
        f = MPlot.Figure(4400+j); clf
        NP.TRF.PlotClusterComboX(sHC(i,:,2), j);
%         MPlot.Paperize(f, 'ColumnsWide', .7, 'Aspect', 1.25, 'FontSize', 5);
        exportgraphics(f, fullfile(anaDir, setName, sprintf("clus_trf_x_sig_%s_by_%s.png", regions(i), mdlNames(j))));
    end
end

