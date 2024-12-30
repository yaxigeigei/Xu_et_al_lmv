%% Extract and cache peri-phoneme sequences

cacheDir = LMV.Data.GetAnalysisDir("linker", "extracted_seq");
srcTb = LMV.Data.FindSource([]);

%% Run extraction for each recording

for i = 1 : height(srcTb)
    % Check computed
    recId = string(srcTb.recId{i});
    cachePath = fullfile(cacheDir, sprintf("%s_seqData.mat", recId));
    if exist(cachePath, 'file')
        fprintf("\nNot overwrite existing cache for %s\n", recId);
        continue
    end
    
    
    % Load se
    sePath = fullfile(LMV.Data.GetAnalysisDir, "linker", "se_m11", srcTb.name{i});
    se11 = NP.SE.LoadSession(sePath);
    clusTb = NP.Unit.GetClusTb(se11);
    
    % Find units that are responsive to both stim and prod
    sTest = LMV.Resp.LoadPhaseResponseTest(recId);
    [~, I] = MMath.SortLike(sTest.clusTb.clusId, clusTb.clusId);
    isResp = all(sTest.sigTb{I, ["stim", "prod"]}, 2);
    
    % Get manually identified units
    cid = arrayfun(@(x) LMV.Linker.GetSelectedClusId(x, recId), LMV.Linker.types, 'Uni', false);
    cid = [cid{:}];
    isSelect = ismember(clusTb.clusId, cid);
    
    % Remove other units
    isUnit = isResp | isSelect;
    if ~any(isUnit)
        continue
    end
    NP.Unit.RemoveUnits(se11, ~isUnit);
    
    % Offload bulky data
    se11.RemoveTable('ni');
    
    
    % Reslice to trials
    [tt, tv] = se11.GetTable("taskTime", "taskValue");
    rt = se11.GetReferenceTime("taskTime");
    tSlice = [0; tt.cue1On(3:2:end)] + rt(1:2:end);
    seTr = se11.Duplicate;
    seTr.SliceSession(tSlice, 'absolute');
    seTr.SetTable("taskValue", tv(2:2:end,:), 'eventValues');
    LMV.SE.StandardizeDataType(seTr);
    
    % Remove trials with bad performance
    tt = seTr.GetTable('taskTime');
    isBad = LMV.SE.IsBadTrials(seTr, -Inf, 0.5);
    seTr.RemoveEpochs(isBad);
    
    % Exclude the incomplete last trial of NP54_B1
    if recId == "NP54_B1"
        seTr.RemoveEpochs(seTr.GetTable("taskValue").trialNum == 72);
    end
    
    
    % Extract peri-phone data
    [featTb, senTb, se] = LMV.PE.ExtractPeriPhoneData(seTr);
    
    % Identify a standard set of phonemes
    isStim = featTb.source=="stim";
    [N, C] = histcounts(categorical(featTb.name(isStim)));
    [N, I] = sort(N, 'descend');
    seeds = string(C(I))';
    
    % Extract sequence data
    phases = ["stim", "prod"];
    sources = ["stim", "prod"];
    positions = { ...
        {[-1 0 1]}, ...
        {[-1 0 1]}, ...
        };
    levels = "word";
    seqData = cell2table(cell(numel(levels), numel(phases)), 'VariableNames', phases, 'RowNames', levels);
    for j = 1 : numel(phases)
        for k = 1 : numel(levels)
            seqData.(phases(j)){k} = LMV.PE.ComputeSeq(featTb, senTb, ...
                'Source', sources(j), 'SeedFeatures', seeds, 'Level', levels(k), 'Positions', positions{j}, 'IncludePartial', true);
        end
    end
    
    
    % Save cache
    save(cachePath, 'se', 'seeds', 'seqData', "-v7.3");
end
