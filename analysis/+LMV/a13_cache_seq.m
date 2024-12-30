%% Extract and cache peri-phoneme sequences

cacheDir = LMV.Data.GetAnalysisDir("peri_event", "extracted_seq");
srcTb = LMV.Data.FindSource([]);

%% Run extraction for each recording

for i = 1 : height(srcTb)
    % Check computed
    recId = NP.SE.GetID(srcTb.name{i});
    cachePath = fullfile(cacheDir, sprintf("%s_seqData.mat", recId));
    if exist(cachePath, 'file')
        fprintf("\nNot overwrite existing cache for %s\n", recId);
        continue
    end
    
    
    % Load se
    sePath = srcTb.path{i};
    se = NP.SE.LoadSession(sePath, 'UserFunc', @(x) x.RemoveTable('LFP', 'niTime')); % also offload bulky data
    
    % Exclude non-responsive units
    sTest = LMV.Resp.LoadPhaseResponseTest(recId);
    isUnit = sTest.sigTb.stim | sTest.sigTb.prod;
    if ~any(isUnit)
        continue
    end
    NP.Unit.RemoveUnits(se, ~isUnit);
    
    % Add features
    ops = NP.Param.Enrich;
    ops.isFiltSpeech = true;
    ops.isMel = true;
    ops.isPitch = true;
    ops.isArtic = true;
    se = LMV.SE.Transform(se, 'enrich', ops);
    
    % Exclude the incomplete last trial of NP54_B1
    if recId == "NP54_B1"
        se.RemoveEpochs(72);
    end
    
    
    % Extract peri-phone data
    [featTb, senTb, se] = LMV.PE.ExtractPeriPhoneData(se);
    
    % Identify a standard set of phonemes
    isStim = featTb.source=="stim";
    [N, C] = histcounts(categorical(featTb.name(isStim)));
    [N, I] = sort(N, 'descend');
    feats = string(C(I));
    
    % Extract sequence data
    phases = ["stim", "feedback", "prod"];
    sources = ["stim", "prod", "prod"];
    positions = { ...
        {0, [0 -1], [0 -1 -2]}, ...
        {0, [0 -1], [0 -1 -2]}, ...
        {0, [0 1], [0 1 2]} ...
        };
    levels = ["phone", "syll", "word"];
    seqData = cell2table(cell(numel(levels), numel(phases)), 'VariableNames', phases, 'RowNames', levels);
    for j = 1 : numel(phases)
        for k = 1 : numel(levels)
            seqData.(phases(j)){k} = LMV.PE.ComputeSeq(featTb, senTb, ...
                'Source', sources(j), 'SeedFeatures', feats, 'Level', levels(k), 'Positions', positions{j});
        end
    end
    
    
    % Save cache
    save(cachePath, 'se', 'feats', 'seqData', "-v7.3");
end
