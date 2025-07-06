%% Create a library of plots for inspecting results

mdlName = LMV.Linker.currentModel;
srcTb = LMV.Data.FindSource([]);

%% Prepare data for plots

cacheDir = LMV.Data.GetAnalysisDir("linker", mdlName, "computed_profile");

sPf = cell(size(srcTb.recId));

for recIdx = 1 : height(srcTb)
    % Get recording ID
    recId = NP.SE.GetID(srcTb.name{recIdx});
    disp(recId);
    cacheFile = fullfile(cacheDir, sprintf("%s.mat", recId));
    if exist(cacheFile, 'file')
        sPf{recIdx} = load(cacheFile);
        continue
    end
    
    
    % Load sequence data (only for plotting)
    seqData = LMV.Linker.LoadSeqData(recId);
    
    % Load ce
    cePath = fullfile(LMV.Data.GetAnalysisDir("linker"), "computed_peth_sem-discounted", recId+"_ce.mat");
    load(cePath, 'ce');
    
    % Load RF models
    rfTb = LMV.TRF.LoadModels("phone_stim", recId);
    [~, I] = MMath.SortLike(rfTb.clusId, ce.clusTb.clusId, false);
    ce.clusTb.phone_stim = rfTb.phone_stim(I);
    
    % Load linker models
    lkTb = LMV.Linker.LM.LoadModels("smooth_lm", recId);
    [~, I] = MMath.SortLike(lkTb.clusId, ce.clusTb.clusId, false);
    ce.clusTb.linker = lkTb.linker(I);
    
    % Only keep units with seq data (i.e. responsive during both stim and prod)
    ceClusId = ce.clusTb.clusId;
    seqClusId = NP.Unit.GetClusTb(seqData.se).clusId;
    [~, uInd] = MMath.SortLike(ceClusId, seqClusId);
    isRm = true(size(ceClusId));
    isRm(uInd) = false;
    ce.RemoveUnits(isRm);
    
    
    % Align to stim onsets
    ce.AlignTime('prodMatchOn', 'taskTime');
    
    % Find linking positions
    LMV.Linker.LM.ComputeLinkingScores(ce);
    
    
    s = struct;
    s.recId = string(recId);
    s.ce = ce;
    s.triTb = seqData.triTb;
    sPf{recIdx} = s;
    
    % Cache results
    save(cacheFile, '-struct', 's');
end

return
%% Add prediction and scores to ce

cacheDir = LMV.Data.GetAnalysisDir("linker", mdlName, "computed_link");

for i = 1 : height(srcTb)
    % Load ce
    recId = srcTb.recId(i);
    cePath = fullfile(LMV.Data.GetAnalysisDir("linker"), "computed_peth_sem-discounted", recId+"_ce.mat");
    load(cePath, 'ce');
    
    % Load RF models (for completeness)
    rfTb = LMV.TRF.LoadModels("phone_stim", recId);
    [~, I] = MMath.SortLike(rfTb.clusId, ce.clusTb.clusId, false);
    ce.clusTb.phone_stim = rfTb.phone_stim(I);
    
    % Load linker models
    lkTb = LMV.Linker.LM.LoadModels(mdlName, recId);
    [~, I] = MMath.SortLike(lkTb.clusId, ce.clusTb.clusId, false);
    ce.clusTb.linker = lkTb.linker(I);
    
    % Remove units without linker model (i.e. those not responsive to either stim or prod)
    isRm = cellfun(@isempty, ce.clusTb.linker);
    ce.RemoveUnits(isRm);
    
    % Align to sentence onsets
    ce.AlignTime('prodMatchOn', 'taskTime');
    
    % Compute scores
    LMV.Linker.LM.ComputeLinkingScores(ce);
    
    % Find linking positions
    posTb = LMV.Linker.LM.FindScorePeaks(ce);
    
    % Cache results
    outPath = fullfile(cacheDir, )
    save(cacheFile, 'ce');
end

