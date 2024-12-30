%% Plot summary goodness-of-fit of all units and all models

anaDir = LMV.Data.GetAnalysisDir('coding', 'stRF');
srcTb = LMV.Data.FindSource([]);

%% Cache data tables for segmental models

cachePath = fullfile(anaDir, "segmental_clusTb.mat");

if ~exist(cachePath, 'file')
    % Load models
    sets = ["phone", "strf", "artic3"];
    targets = LMV.TRF.targets;
    mdlNames = sets+"_"+targets(:);
    clusTbs = LMV.TRF.LoadModels(mdlNames, srcTb.recId);
    clusTb = cat(1, clusTbs{:});
    
    % Find segmental models and their r2 values
    C = clusTb{:,mdlNames};
    r2 = NaN(size(C));
    pval = r2;
    for i = 1 : numel(C)
        if isempty(C{i})
            continue
        end
        
        C{i} = LMV.RF.UpdateModelR2(C{i});
        
        if C{i}.r2 < 0.01 || C{i}.null.r2Pval > 0.05
            C{i} = []; % remove non-segmental or insignificant models
        else
            r2(i) = C{i}.r2;
            pval(i) = C{i}.null.r2Pval;
        end
    end
    clusTb{:,mdlNames} = C;
    
    r2Tb = array2table(r2, 'VariableNames', mdlNames, 'RowNames', "u"+clusTb.clusId);
    pvalTb = array2table(pval, 'VariableNames', mdlNames, 'RowNames', "u"+clusTb.clusId);
    
    save(cachePath, 'sets', 'targets', 'mdlNames', 'clusTb', 'r2Tb', 'pvalTb');
else
    load(cachePath);
end

%% Load test results

rTest = LMV.Resp.LoadPhaseResponseTest();
[~, I] = MMath.SortLike(rTest.clusTb.clusId, clusTb.clusId);
rSig = rTest.sigTb{I,["stim" "prod" "prod"]};

% sTest = LMV.Resp.LoadSentenceSelectTest();
% [~, I] = MMath.SortLike(sTest.clusTb.clusId, clusTb.clusId);
% sSig = sTest.sigTb{I,["stim" "prod" "prod"]};

%% Compute the fraction of explainable units with bootstrap

% Group vPrCG and mPrCG as PrCG
clusTb.region2 = clusTb.region;
clusTb.region2(ismember(clusTb.region, ["vPrCG", "mPrCG"])) = "PrCG";

regions = ["PrCG", "IFG", "STG"];
targets = LMV.TRF.targets;
sets = {["strf", "phone"], ["strf", "phone"], ["artic3", "phone"]};

frac = cell(size(targets));
rng(61);

for i = 1 : numel(targets)
    ff = NaN(numel(regions), numel(sets{i}), 3);
    for j = 1 : numel(regions)
        isRegion = clusTb.region2 == regions(j);
        for k = 1 : numel(sets{i})
            % Find the subset of interest
            isResp = rSig(:,i) > 0;
            isUnit = isRegion & isResp;
            
            % Compute the fraction
            mn = sets{i}(k)+"_"+targets(i);
            isSeg = ~isnan(r2Tb.(mn)(isUnit));
            
            % ff(j,k,1) = mean(isSeg);
            % G = findgroups(clusTb.recId(isUnit));
            % ff(j,k,2:3) = MMath.BootCI(500, @mean, isSeg, 'Groups', G);
            
            [ff(j,k,1), ~, ~, ff(j,k,2:3)] = MMath.MeanStats(isSeg, 1, 'NBoot', 1000);
        end
    end
    frac{i} = ff;
end

%% Plot fraction of units explainable

f = MPlot.Figure(6360); clf
tl = tiledlayout("vertical");
tl.Padding = "compact";

for i = 1 : numel(targets)
    ax = nexttile;
    
    ff = frac{i};
    bb = bar(ff(:,:,1)); hold on
    for j = 1 : numel(bb)
        bb(j).FaceColor = LMV.Param.GetModelColors(sets{i}(j));
        errorbar(bb(j).XEndPoints, bb(j).YData, bb(j).YData'-ff(:,j,2), ff(:,j,3)-bb(j).YData', 'k.');
    end
    
    legend(bb, replace(sets{i}, ["strf", "artic3", "phone"], ["spectral", "artic", "phone"]), 'Location', 'north');
    
    ax.XTickLabel = regions;
    ax.XTickLabelRotation = 0;
    ax.YLim = [0 .3];
    ax.YTick = 0:.1:ax.YLim(2);
    ax.Title.String = targets(i);
    ax.Title.Interpreter = "none";
    ax.YLabel.String = "Frac. of resp. units";
    MPlot.Axes(ax);
end

MPlot.Paperize(f, 'ColumnsWide', .5, 'ColumnsHigh', .8);
% exportgraphics(f, fullfile(anaDir, "stRF_frac_explainable.png"));
% exportgraphics(f, fullfile(anaDir, "stRF_frac_explainable.pdf"));

%% Plot overlay of r2 timecourse

sets = [repelem("strf",2), "artic3", repelem("phone",3)];
targets = [LMV.TRF.targets LMV.TRF.targets];
regions = LMV.Param.regions(1:2);

f = MPlot.Figure(3292); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";

for k = 1 : numel(sets)
    % Find models
    mn = sets(k)+"_"+targets(k);
    isMdl = false(height(clusTb),1);
    for i = 1 : height(clusTb)
        if isempty(clusTb.(mn){i})
            continue
        end
        
        % Disambiguate feedback vs prod by testing if the model has higher r2
        if targets(k)=="stim" || isempty(clusTb.(sets(k)+"_feedback"){i}) || isempty(clusTb.(sets(k)+"_prod"){i})
            isMdl(i) = true;
        elseif targets(k)=="prod"
            isMdl(i) = clusTb.(mn){i}.r2 > clusTb.(sets(k)+"_feedback"){i}.r2;
        elseif targets(k)=="feedback"
            isMdl(i) = clusTb.(mn){i}.r2 > clusTb.(sets(k)+"_prod"){i}.r2;
        end
    end
    isRegion = ismember(clusTb.region, regions);
    mdls = clusTb.(mn)(isMdl & isRegion);
    
    % Sort models
    [~, I] = sort(cellfun(@(x) x.r2t, mdls), 'descend');
    mdls = mdls(I);
    
    % Print the best model
    [~, I] = max(cellfun(@(x) x.r2, mdls));
    disp(mdls{I}.resps);
    
    % Make plot
    ax = nexttile;
    LMV.RF.PlotR2Overlay(mdls, 'Parent', ax);
    ax.Title.String = sprintf("%s %s", sets(k), targets(k));
end

% MPlot.Paperize(f, 'ColumnsWide', 1.6, 'ColumnsHigh', 0.7);
% exportgraphics(f, fullfile(anaDir, "stRF_r2_ladder.png"));

return
%% Plot violin scatters for inspection

targets = LMV.TRF.targets;

f = MPlot.Figure(6364); clf
tl = tiledlayout("vertical");
tl.Padding = "compact";

LMV.TRF.PlotR2Scatter("phone", targets, clusTb, r2Tb, rSig);

LMV.TRF.PlotR2Scatter("artic3", targets, clusTb, r2Tb, rSig);

MPlot.Paperize(f, 'ColumnsWide', 1, 'ColumnsHigh', 0.7);
% exportgraphics(f, fullfile(anaDir, "stRF_r2_violin.png"));

%% Plot r2 CDFs

sets = [repelem("strf",2), "artic3", repelem("phone",3)];
targets = [LMV.TRF.targets LMV.TRF.targets];
regions = LMV.Param.regions;

nr = numel(regions);
nm = numel(sets);

f = MPlot.Figure(6363); clf
tl = tiledlayout(2, 3);
tl.Padding = "compact";

for i = 1 : nm
    mn = sets(i)+"_"+targets(i);
    ax = nexttile;
    for j = 1 : nr
        % Get r2 values
        r2 = r2Tb.(mn);
        
        % Find the subset of interest
        isRegion = clusTb.region == regions(j);
        isResp = rSig(:,mod(i-1,3)+1) > 0;
        isUnit = isRegion & isResp;
        
        binEdges = 0:0.01:1;
        histogram(r2(isUnit), binEdges, 'Normalization', 'cdf', 'DisplayStyle', 'stairs', 'EdgeColor', LMV.Param.GetRegionColors(regions(j)));
        hold on
    end
    ax.XLim = [0 0.4];
    ax.YLim = [0 0.3];
    ax.Title.String = mn;
    ax.Title.Interpreter = "none";
    ax.XLabel.String = "r-squared";
    ax.YLabel.String = "Frac. of resp. units";
    MPlot.Axes(ax);
end

MPlot.Paperize(f, 'ColumnsWide', 1, 'ColumnsHigh', 0.7);
% exportgraphics(f, fullfile(anaDir, "stRF_r2_cdf.png"));

%% 

cid = [460100250 520200153 540100167 410100245 540100167];
LMV.TRF.PlotWeightsComboFromCache(cid, 'm1', "phone_"+LMV.TRF.targets);
