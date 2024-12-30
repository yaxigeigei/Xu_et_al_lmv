%% Multi-phoneme modulation

anaDir = LMV.Data.GetAnalysisDir("peri_event");
srcTb = LMV.Data.FindSource([]);

%% Extract sequence responses

cachePath = fullfile(anaDir, "seq_resp_recTb.mat");

if isfile(cachePath)
    load(cachePath, "recTb");
else
    recTb = srcTb;
    for recIdx = 1 : height(srcTb)
        % Load cached sequence data
        recId = NP.SE.GetID(srcTb.name{recIdx});
        inputFile = fullfile(LMV.Data.GetAnalysisDir, "peri_event", "extracted_seq", recId+"_seqData.mat");
        if ~isfile(inputFile)
            fprintf("%s does not have cached sequence data.\n", recId);
            continue
        end
        s = load(inputFile, "se", "seqData");
        
        % Extract sequence responses
        targets = LMV.TRF.targets;
        recTb{recIdx,targets} = LMV.PE.ExtractSeqResp(s, targets);
    end
    save(cachePath, "recTb");
end

%% Compute modulation

regions = LMV.Param.regions;
isRegion = ismember(recTb.Region, regions(1:2));

targets = LMV.TRF.targets;
clusTbs = cell(size(targets));

for tIdx = 1 : numel(targets)
    target = targets(tIdx);
    disp(target);
    
    clusTb = cat(1, recTb.(target){isRegion});
    clusTb.mod = NaN(height(clusTb), 3);
    for u = 1 : height(clusTb)
        disp("u"+clusTb.clusId(u));
        
        % Find the best tuned seed
        stRF = clusTb.mdls{u};
        opRF = stRF.mdls{stRF.r2Idx};
        [~, I] = max(opRF.Beta);
        seed = stRF.feats(I);
        
        % Compute modulation
        rTb = clusTb.seqR{u};
        isSeed = rTb.seqStr1==seed;
        rTb = rTb(isSeed,:);
        
        mi = LMV.PE.ComputeSeqMod(rTb);
        nShuffle = 200;
        miNull = zeros(nShuffle, numel(mi));
        for n = 1 : nShuffle
            miNull(n,:) = LMV.PE.ComputeSeqMod(rTb, true);
        end
        miPval = MMath.EstimatePval(mi, miNull)';
        
        [kwPval, kwTb] = arrayfun(@(x) kruskalwallis(rTb.optiResp, rTb.(x), 'off'), "seqStr"+(1:3), 'Uni', false);
        kwPval = cat(2, kwPval{:});
        kwSize = cellfun(@(x) x{2,5}, kwTb);
        
        clusTb.seed(u) = seed;
        clusTb.modPval(u,:) = miPval;
        clusTb.mod(u,:) = mi;
        clusTb.kwPval(u,:) = kwPval;
        clusTb.kwSize(u,:) = kwSize;
    end
    clusTbs{tIdx} = clusTb;
end

%% Make plots and print stats

modFile = fullfile(anaDir, "seq_mod.txt");
fileID = fopen(modFile, "w");

f = MPlot.Figure(830); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";

for i = 1 : numel(clusTbs)
    mi = clusTbs{i}.mod;
    miPval = clusTbs{i}.modPval;
    mi(miPval>0.05) = NaN;
    mi = mi + 1;
    
    m2 = rmmissing(mi(:,2));
    [~, q] = iqr(m2);
    pval = signrank(m2);
    fprintf(fileID, "%s: 2-phone, n = %i units, median (IQR): %.1f (%.1f-%.1f), p = %.1e\n", targets(i), numel(m2), median(m2), q(1), q(2), pval);
    
    m3 = rmmissing(mi(:,3));
    [~, q] = iqr(m3);
    pval = signrank(m3);
    fprintf(fileID, "%s: 3-phone, n = %i units, median (IQR): %.1f (%.1f-%.1f), p = %.1e\n\n", targets(i), numel(m3), median(m3), q(1), q(2), pval);
    
    ax = nexttile;
    nPh = 2:3;
    boxplot(mi(:,nPh), 'Colors', 'k');
    % ax.YLim(1) = 1; % modulation depth is by definition ≥ 1
    ax.YTick = 1:2:ax.YLim(2);
    ax.XTickLabel = nPh+"-phone";
    ax.YLabel.String = "Modulation depth";
    ax.Title.String = targets(i);
    MPlot.Axes(ax);
end

fclose(fileID);

MPlot.Paperize(f, 1.2, 0.3);
exportgraphics(f, fullfile(anaDir, "seq_mod.png"));
exportgraphics(f, fullfile(anaDir, "seq_mod.pdf"));

return
%% 

regions = LMV.Param.regions;
isRegion = ismember(recTb.Region, regions(1:2));

targets = LMV.TRF.targets;
clusTbs = cell(size(targets));

for tIdx = 1 : numel(targets)
    target = targets(tIdx);
    clusTb = cat(1, recTb.(target){isRegion});
    for i = 1 : height(clusTb)
        % Compute modulation
        rTb = clusTb.seqR{i};
        return
        % M = cellfun(@(x) std(x), rTb{:,:}, 'Uni', false);
        M = cellfun(@(x) max(x)-min(x), rTb{:,:}, 'Uni', false);
        % M(cellfun(@isempty, M)) = {NaN};
        M = cell2mat(M);
        M = M./cat(1, rTb.(1){:});
        M(isinf(M)) = NaN;
        mTb = array2table(M, "RowNames", rTb.Properties.RowNames, "VariableNames", rTb.Properties.VariableNames);
        
        % Find the best tuned seed
        stRF = clusTb.mdls{i};
        opRF = stRF.mdls{stRF.r2Idx};
        [~, I] = max(opRF.Beta);
        seed = stRF.feats(I);
        clusTb.topSeed(i) = seed;
        clusTb.topMod(i,:) = mTb{seed,:};
    end
    clusTbs{tIdx} = clusTb;
end

%% 

f = MPlot.Figure(830); clf
tl = tiledlayout("vertical");
tl.Padding = "compact";
for i = 1 : numel(clusTbs)
    
    mi = clusTbs{i}.mod;
    p = clusTbs{i}.modPval;
    mi(p>0.05) = NaN;
    mi = mi+1;
    
    ax = nexttile;
    plot(mi', 'Color', [0 0 0 .3]);
    ax.XLim = [.5 3.5];
    ax.XTick = 1:3;
    ax.XTickLabel = ax.XTick;
    ax.YScale = 'log';
    ax.YTick = [1 2 4 8];
    ax.YLabel.String = "Rel. mod. depth";
    ax.Title.String = targets(i);
    MPlot.Axes(ax);
end
MPlot.Paperize(f, 0.4, 0.8);

