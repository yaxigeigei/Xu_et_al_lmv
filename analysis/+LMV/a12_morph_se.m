%% Transform each se such that the prod times are aligned to stim times

cacheDir = LMV.Data.GetAnalysisDir("linker", 'se_m11');
srcTb = LMV.Data.FindSource([]);

%% Make morphed se

for k = 1 : height(srcTb)
    % Check computed
    m11Path = fullfile(cacheDir, srcTb.name{k});
    if exist(m11Path, 'file')
        fprintf("\nSkip cached se\n%s\n", m11Path);
        continue
    end
    
    % Load and enrich se
    se = NP.SE.LoadSession(srcTb.path{k}, 'UserFunc', @(x) x.RemoveTable('LFP', 'niTime'));
    ops = NP.Param.Enrich;
    ops.isSpkRate = true;
    ops.spkBinSize = 0.01;
    ops.spkKerSize = 0.015;
    ops.isSpkSpan = true;
    ops.isMel = true;
    se = LMV.SE.Transform(se, 'enrich', ops);
    
    % Morph prod to stim
    se = LMV.Linker.MorphSession(se);
    
    % Save morphed se
    save(m11Path, 'se', '-v7.3');
end

%% Make seTb

seTbDir = LMV.Data.GetAnalysisDir("linker", "computed_seTb");

for k = 1 : height(srcTb)
    % Check computed
    seTbPath = fullfile(seTbDir, strrep(srcTb.name{k}, '_se.mat', '_seTb.mat'));
    if exist(seTbPath, 'file')
        fprintf("\nSkip computed seTb\n%s\n", seTbPath);
        continue
    end
    
    
    % Load stim-prod aligned se
    se = NP.SE.LoadSession(fullfile(anaDir, "se_m11", srcTb.name{k}));
    
    % Slice epochs to include a period before and after
    tOffset = [-1 1]*0.5; % in sec
    se = NP.SE.BleedTrials(se, tOffset);
    NP.Unit.SetUniqueClusId(se); % make sure cluster IDs are unique
    
    % Remove pairs of epochs with bad performance
    tt = se.GetTable('taskTime');
    isBad = LMV.SE.IsBadTrials(se, -Inf, 0.5);
    isBad = isBad | isnan(tt.prodMatchOff);
    isBad = isBad | circshift(isBad, -1); % get the pairs
    se.RemoveEpochs(isBad);
    
    
    % Compute cross-correlation between stim and prod Mel sprectrogram
    ops = NP.Param.Resample();
    ops.rsBinSize = 0.01;
    ops.rsShifts = -0.5 : ops.rsBinSize : 0.1; % amount to shift prod sampling window
    se.userData.xrOps = ops;
    xrMel = LMV.Linker.ComputeMelXC(se, ops);
    se.userData.xrMel = xrMel;
    
    % Make seTb
    seTb = NP.SE.SplitConditions(se, 'conditionVars', {'taskName', 'phase', 'stimText'});
    for i = 2 : height(seTb)
        seTb.se(i).userData = [];
    end
    
    % Save seTb
    seTbPath = fullfile(seTbDir, strrep(srcTb.name{k}, '_se.mat', '_seTb.mat'));
    save(seTbPath, 'seTb', '-v7.3');
end

