%% Fit linear models that predict prod activity from temporal stim activity for each unit

mdlName = LMV.Linker.currentModel;
mdlDir = LMV.Data.GetAnalysisDir("linker", mdlName, "mdls");
srcTb = LMV.Data.FindSource([]);

%% 

for recIdx = 1 : height(srcTb)
    % Get recording ID
    recId = NP.SE.GetID(srcTb.name{recIdx});
    disp(recId);
    
    % Check models
    mdlFile = sprintf("%s_clusTb.mat", recId);
    mdlPath = fullfile(mdlDir, mdlFile);
    if exist(mdlPath, 'file')
        fprintf("\nSkip computed models from\n%s\n", mdlPath);
        continue
    end
    
    % Load ce
    cePath = fullfile(LMV.Data.GetAnalysisDir("linker"), "computed_peth_sem-discounted", recId+"_ce.mat");
    load(cePath, 'ce');
    
    % Find units responsive to both stim and prod
    rTest = LMV.Resp.LoadPhaseResponseTest(recId);
    [~, I] = MMath.SortLike(rTest.clusTb.clusId, ce.clusTb.clusId);
    isResp = all(rTest.sigTb{I, ["stim", "prod"]}, 2);
    
    % Find manually identified units
    cid = arrayfun(@(x) LMV.Linker.GetSelectedClusId(x, recId), LMV.Linker.types, 'Uni', false);
    cid = [cid{:}];
    isSelect = ismember(ce.clusTb.clusId, cid);
    
    % Fit models
    uInd = find(isResp | isSelect);
    tic
    clusTb = LMV.Linker.LM.FitEncoding(ce, uInd);
    toc
    
    % Save models
    save(mdlPath, 'clusTb');
end

return
%% 

mdlDir = LMV.Data.GetAnalysisDir("linker", "smooth_lm_trials", "mdls");
srcTb = LMV.Data.FindSource([]);

%% 

for recIdx = 1 : height(srcTb)
    % Get recording ID
    recId = NP.SE.GetID(srcTb.name{recIdx});
    disp(recId);
    
    % Check models
    mdlFile = sprintf("%s_clusTb.mat", recId);
    mdlPath = fullfile(mdlDir, mdlFile);
    if exist(mdlPath, 'file')
        fprintf("\nSkip computed models from\n%s\n", mdlPath);
        continue
    end
    
    
    % Load se
    sePath = fullfile(LMV.Data.GetAnalysisDir, "linker", "se_m11", recId+"_se.mat");
    se = NP.SE.LoadSession(sePath);
    
    % % Slice epochs to include a period before and after
    % tOffset = [-1 1]*0.5; % in sec
    % se = NP.SE.BleedTrials(se, tOffset);
    % NP.Unit.SetUniqueClusId(se); % make sure cluster IDs are unique
    
    % Remove pairs of epochs with bad performance
    tt = se.GetTable('taskTime');
    isBad = LMV.SE.IsBadTrials(se, -Inf, 0.5);
    isBad = isBad | isnan(tt.prodMatchOff);
    isBad = isBad | circshift(isBad, -1); % get the pairs
    se.RemoveEpochs(isBad);
    
    
    % Resample spike rates
    tt = se.GetTable('taskTime');
    ops = NP.Param.Resample();
    ops.rsWin = round([tt.trialOn tt.prodMatchOff] + [0 .5], 6); % round off eps introduced during morphing
    ops.rsBinSize = 0.01;
    tEdges = arrayfun(@(a,b) a : ops.rsBinSize : b, ops.rsWin(:,1), ops.rsWin(:,2), 'Uni', false);
    respTb = se.ResampleTimeSeries("spikeRate", tEdges);
    
    % Make ce
    ce = NP.CodingExplorer(se);
    ce.SetTable("resp", respTb, "timeSeries", se.GetReferenceTime);
    
    
    % Find units responsive to both stim and prod
    rTest = LMV.Resp.LoadPhaseResponseTest(recId);
    [~, I] = MMath.SortLike(rTest.clusTb.clusId, ce.clusTb.clusId);
    isResp = all(rTest.sigTb{I, ["stim", "prod"]}, 2);
    
    % Find manually identified units
    cid = arrayfun(@(x) LMV.Linker.GetSelectedClusId(x, recId), LMV.Linker.types, 'Uni', false);
    cid = [cid{:}];
    isSelect = ismember(ce.clusTb.clusId, cid);
    
    
    % Fit models
    uInd = find(isResp | isSelect);
    tic
    clusTb = LMV.Linker.LM.FitEncoding(ce, uInd);
    toc
    
    % Save models
    save(mdlPath, 'clusTb');
end



