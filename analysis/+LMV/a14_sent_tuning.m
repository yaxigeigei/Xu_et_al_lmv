%% Compute sentence tuning

anaDir = LMV.Data.GetAnalysisDir('sent_resp', 'tuning');

srcTb = LMV.Data.FindSource([]);

%% Load test data

resp = LMV.Resp.LoadPhaseResponse(srcTb.recId);
rTest = LMV.Resp.LoadPhaseResponseTest(srcTb.recId);
sTest = LMV.Resp.LoadSentenceSelectTest(srcTb.recId);

%% Compute sentence tuning

phaseNames = ["stim", "delay", "init", "prod"];
stimIdList = LMV.Param.stimIdList12;

nRec = height(srcTb);
nPhase = numel(phaseNames);
nSen = numel(stimIdList);
clusTbs = cell(nRec, nPhase);

for i = 1 : nRec
    fprintf("\nCompute sentence tunings for %s\n", srcTb.recId{i});
    for j = 1 : nPhase
        % Unpack variables
        respTb = resp(i).respTb;
        dropoutTb = resp(i).dropoutTb;
        phaseTb = resp(i).phaseTb;
        stimIdTb = resp(i).stimIdTb;
        clusTb = resp(i).clusTb;
        
        % Add phase name to clusTb
        pn = phaseNames(j);
        clusTb.tuningPhase(:) = pn;
        
        % Set dropout periods to NaN
        respTb{:,:}(dropoutTb{:,:}) = NaN;
        
        % Keep responsive and sentence selective units
        isInc = rTest(i).sigTb.(pn) & sTest(i).sigTb.(pn);
%         isInc = rTest(i).sigTb.(pn) > 0;
        respTb = respTb(:,isInc);
        clusTb = clusTb(isInc,:);
        
        % 
        T = NaN(width(respTb), nSen, 4);
        for k = 1 : nSen
            varIdx = find(stimIdTb.Properties.VariableNames == stimIdList(k));
            if isempty(varIdx)
                if j == 1
                    fprintf("%s does not have the sentence %s\n", srcTb.recId{i}, stimIdList(k));
                end
                continue
            end
            isSen = stimIdTb.(varIdx);
            dResp = respTb{isSen & phaseTb.(pn), :} - respTb{isSen & phaseTb.baseline, :};
            
            % [T(:,k,1), T(:,k,2), T(:,k,3)] = MMath.MeanStats(dResp, 1);
            [T(:,k,1), ~, T(:,k,2)] = MMath.MedianStats(dResp, 1);
            T(:,k,4) = sum(~isnan(dResp), 1);
        end
        clusTb.tuningMean = T(:,:,1);
        clusTb.tuningErr = T(:,:,3);
        clusTb.tuningN = T(:,:,4);
        clusTbs{i,j} = clusTb;
    end
end

clusTbCat = cat(1, clusTbs{:});

%% Plot tuning heatmaps

% regions = LMV.Param.regions;
regions = "mPrCG";
nReg = numel(regions);

f = MPlot.Figure(5656); clf
f.WindowState = 'maximized';
tl = tiledlayout(nReg, nPhase);
tl.Padding = 'compact';

for i = 1 : nReg
    for j = 1 : nPhase
        ax = nexttile;
        
        rn = regions(i);
        pn = phaseNames(j);
        isUnit = clusTbCat.region == rn & clusTbCat.tuningPhase == pn;
        if ~any(isUnit)
            ax.Visible = 'off';
            continue
        end
        M = clusTbCat.tuningMean(isUnit,:)';
        cid = clusTbCat.clusId(isUnit);
        
        M = MMath.Normalize(M, 'minmax');
        % M = softmax(M);
        [M, I] = LMV.NMF.SortPatterns(M, 'peak');
        cid = cid(I);
        
        imagesc(ax, M');
        ax.XTick = 1 : nSen;
        ax.XTickLabel = "s"+ax.XTick;
        ax.YTick = 1 : sum(isUnit);
        ax.YTickLabel = "u"+cid;
        ax.Title.String = sprintf("%s: %s", rn, pn);
        MPlot.Axes(ax);
    end
end

exportgraphics(f, fullfile(anaDir, sprintf("sentence_tuning_heatmaps_%s.png", strjoin(regions, "_"))));

%% Plot 12-sentence spike rasters

% regions = LMV.Param.regions;
regions = "mPrCG";
nReg = numel(regions);

for i = 1 : nReg
    for j = 1 : nPhase
        
        rn = regions(i);
        pn = phaseNames(j);
        rasterDir = fullfile(anaDir, rn+"_"+pn);
        if ~exist(rasterDir, 'dir')
            mkdir(rasterDir);
        end
        
        isUnit = clusTbCat.region == rn & clusTbCat.tuningPhase == pn;
        if ~any(isUnit)
            continue
        end
        M = clusTbCat.tuningMean(isUnit,:)';
        cid = clusTbCat.clusId(isUnit);
        
        M = MMath.Normalize(M, 'minmax');
        M = softmax(M);
        [M, I] = LMV.NMF.SortPatterns(M, 'peak');
        cid = cid(I);
        
        cidArr = NP.UnitPlot.MakeArrayID(cid, 10);
        
        for k = 1 : size(cidArr,2)
            cidNames = "u"+cidArr(:,k);
            cidNames(ismissing(cidNames)) = [];
            figPath = fullfile(rasterDir, sprintf("page%i_%s.png", k, strjoin(cidNames, "_")));
            if exist(figPath, 'file')
                continue
            end
            
            f = MPlot.Figure(6656); clf
            f.WindowState = 'maximized';
            phaseArg = pn;
            if ~ismember(phaseArg, ["stim", "prod"])
                phaseArg = "full";
            end
%             LMV.Overview.SentencesFromCache(cidArr(:,k), 'StimIdList', stimIdList, 'TaskPhase', phaseArg);
            LMV.Overview.SessionFromCache(cidArr(:,k), 'DataSource', "m2", 'StimIdList', stimIdList);
            exportgraphics(f, figPath);
        end
    end
end

return
%% Plot tuning violin-scatter plots

regions = LMV.Param.regions;
nReg = numel(regions);

f = MPlot.Figure(5657); clf
f.WindowState = 'maximized';
tl = tiledlayout(nReg, nPhase);
tl.Padding = 'compact';

for i = 1 : nReg
    for j = 1 : nPhase
        ax = nexttile;
        
        rn = regions(i);
        pn = phaseNames(j);
        isUnit = clusTbCat.region == rn & clusTbCat.tuningPhase == pn;
        if ~any(isUnit)
            ax.Visible = 'off';
            continue
        end
        M = clusTbCat.tuningMean(isUnit,:)';
        
        M = MMath.Normalize(M, 'minmax');
        M = softmax(M);
        for k = 1 : nSen
            MPlot.ViolinScatter(k, M(k,:));
            hold(ax, 'on');
        end
        
        ax.XLim = [0 nSen+1];
        ax.YLim = [.04 .2];
        ax.YLabel.String = "softmax()";
        ax.XTick = 1 : nSen;
        ax.XTickLabel = "s"+ax.XTick;
        ax.Title.String = sprintf("%s: %s", rn, pn);
        MPlot.Axes(ax);
    end
end

