%% Export data and metadata for contrastive learning analysis in Python

srcTb = LMV.Data.FindSource([]);

seDir = LMV.Data.GetAnalysisDir("data", "se_m11");
pklDir = LMV.Data.GetAnalysisDir("data", "pkl_m11");

% pyenv("Version", "C:\Users\many\miniconda3\envs\lmv\python.exe")

%% Exports unit responses and task intervals for construction of se objects in Python

for i = 1 : height(srcTb)
    % Load session
    recId = srcTb.recId(i);
    sePath = fullfile(seDir, recId+"_se.mat");
    se = NP.SE.LoadSession(sePath);
    
    % Reslice se to a single epoch
    se1 = se.Duplicate({'spikeRate'});
    se1.SliceSession(0, 'absolute');
    
    % Make spike rate table
    srTb = se1.GetTable("spikeRate");
    srCell = srTb{:,:};
    srTb = table(srCell{:}, 'VariableNames', srTb.Properties.VariableNames);
    
    % Make matched stim and prod interval table
    [ttTb, tvTb] = se.GetTable("taskTime", "taskValue");
    rt = se.GetReferenceTime();
    isStim = tvTb.phase == "stim";
    isProd = tvTb.phase == "prod";

    stimTb = table( ...
        ttTb.stimOn(isStim) + rt(isStim), ...
        ttTb.stimOff(isStim) + rt(isStim), ...
        tvTb.stimId(isStim), ...
        tvTb.stimText(isStim), ...
        'VariableNames', ["start_time", "stop_time", "stim_id", "transcript"]);

    prodTb = table(...
        ttTb.stimOn(isStim) + rt(isProd), ... % use the same stimOn and stimOff but for prod phase by adding rt(isProd)
        ttTb.stimOff(isStim) + rt(isProd), ...
        tvTb.stimId(isProd), ...
        tvTb.prodText(isProd), ...
        tvTb.alignScore(isProd), ...
        tvTb.tReact(isProd), ...
        'VariableNames', ["start_time", "stop_time", "stim_id", "transcript", "align_score", "reaction_time"]);

    isMatched = all(~isnan([stimTb.start_time, stimTb.stop_time, prodTb.start_time, prodTb.stop_time]), 2);
    stimTb = stimTb(isMatched,:);
    prodTb = prodTb(isMatched,:);

    % Save tables as dataframes
    outDir = fullfile(pklDir, recId);
    if ~isfolder(outDir)
        mkdir(outDir);
    end
    py.pandas.DataFrame(srTb).to_pickle(fullfile(outDir, "spike_rate.pkl.gz"));
    py.pandas.DataFrame(stimTb).to_pickle(fullfile(outDir, "matched_stim.pkl"));
    py.pandas.DataFrame(prodTb).to_pickle(fullfile(outDir, "matched_prod.pkl"));
end

%% Export unit info

% Load linker info
tbLink = LMV.Linker.LoadClusTb();

for i = 1 : height(srcTb)
    % Load task phase responsiveness test results
    recId = srcTb.recId{i};
    sResp = LMV.Resp.LoadPhaseResponseTest(recId);
    clusTb = [sResp.clusTb(:,["clusId", "recId", "subjectId", "region"]) sResp.sigTb];
    
    % Add linker info
    clusTb.linkerGroup(:) = "none";
    for j = 1 : height(tbLink)
        m = tbLink.clusId(j) == clusTb.clusId;
        if ~any(m)
            continue
        end
        clusTb.linkerGroup(m) = tbLink.hcGroup(j);
    end
    
    % Save clusTb as csv file
    csvPath = fullfile(pklDir, recId, "units.csv");
    writetable(clusTb, csvPath);
end


