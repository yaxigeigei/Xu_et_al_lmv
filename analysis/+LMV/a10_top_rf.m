%% Plot summary goodness-of-fit of all units and all models

anaDir = LMV.Data.GetAnalysisDir('coding', 'stRF');
srcTb = LMV.Data.FindSource([]);

%% Load test results

rTest = LMV.Resp.LoadPhaseResponseTest();
clusTb = rTest.clusTb;
rSig = rTest.sigTb{:,["stim" "prod" "prod"]};

sTest = LMV.Resp.LoadSentenceSelectTest();
sSig = sTest.sigTb{:,["stim" "prod" "prod"]};

%% Load models and extract r values

sets = ["phone", "strf", "artic3"];
targets = LMV.TRF.targets;

nu = height(clusTb);
ns = numel(sets);
nt = numel(targets);

setMdls = cell(size(sets));
rrr = NaN(nu, ns, nt);
ppp = rrr;
for i = 1 : ns
    % Load models
    mdlNames = sets(i) + "_" + targets;
    mdlDir = fullfile(anaDir, sets(i), "mdls");
    mdlTbs = LMV.TRF.LoadModels(mdlNames);
    mdlTb = cat(1, mdlTbs{:});
    
    % Make unit order in mdlTb consistent with testTb
    [~, I] = MMath.SortLike(mdlTb.clusId, clusTb.clusId);
    mdlTb = mdlTb(I,:);
    
    % Get computed models
    mdls = mdlTb{:,mdlNames};
    hasMdl = ~cellfun(@isempty, mdls);
    hasNull = cellfun(@(x) isfield(x, 'null'), mdls);
    
    % Extract r
    rr = NaN(nu,nt);
    rr(hasMdl) = cellfun(@(x) LMV.RF.FindPeakR2(x.r2each, x.dt), mdls(hasMdl));
    rr(rr<0) = 0;
    rr = sqrt(rr);
    rrr(:,i,:) = rr;
    
    % Extract pvals
    pp = NaN(nu,nt);
    pp(hasNull) = cellfun(@(x) x.null.r2Pval, mdls(hasNull));
    ppp(:,i,:) = pp;
    
    setMdls{i} = mdlTb;
end

%% Inspect RF weights of top units

targets = ["stim", "feedback", "prod"];
regions = LMV.Param.regions;

f = MPlot.Figure(3298);

for i = 1 : numel(sets)
    mdlTb = setMdls{i};
    for j = 1 : numel(targets)
        mdlNames = sets(i) + "_" + targets;
        for k = 1 : numel(regions)
            isReg = clusTb.region == regions(k);
            isSig = rSig(:,j) > 0;
            uInd = find(isReg & isSig);
            if isempty(uInd)
                continue
            end
            
            [~, I] = sort(rrr(uInd,i,j), 'descend');
            uInd = uInd(I);
            uIndArray = NP.UnitPlot.MakeArrayID(uInd, 15);
            
            clf(f);
            LMV.TRF.PlotWeightsCombo(mdlTb(uIndArray(:,1),:), mdlNames);
            figName = sprintf("%s_trf_top_%s.png", regions(k), targets(j));
            figDir = fullfile(anaDir, sets(i));
            exportgraphics(f, fullfile(figDir, figName));
            % return
        end
    end
end
