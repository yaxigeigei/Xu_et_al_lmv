%% Single-unit classification of sentences at each task phase

% Analysis folder
anaDir = fullfile(NP.Data.GetAnalysisRoot, 'multiphase', 'extracted_feat');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end

% Source data
srcTb = NP.Data.FindSource('lmv');

%% 

for i = 1 : height(srcTb)
    %% Load and standardize se
    
    % Check file
    cePath = fullfile(anaDir, strrep(srcTb.name{i}, 'se', 'ce'));
    if exist(cePath, 'file')
        fprintf("Not overwriting the existing ce file at\n%s\n\n", srcTb.path{i});
        continue
    end
    
    % Load recording
    se = NP.SE.LoadSession(srcTb.path{i});
    
    %% Add derived features
    
    % Standard enrichment
    ops = NP.Param.Enrich;
    ops.isMel = true;
    ops.isPitch = true; % this also adds 'trial' table, no need to run separately
    ops.isArtic = true;
    se = LMV.SE.Transform(se, 'enrich', ops); % this only keeps LMV trials
    
    % Add word level features
    LMV.Word.AddWordEvents(se, "lmv_syllabification_lookup.csv");
    LMV.Word.AddWordTimeseriesTable(se);
    
    %% Resampling
    
    % Initialize option struct from the enrichment option struct
    ops = NP.Param.Resample(ops);
    
    % Specify the grouping event
    evtNames = {'stim', 'prod'};
    sEvents = struct('tableName', 'taskTime', 'varName', evtNames, 'readout', 0);
    ops.eventVar = sEvents;
    
    % Specify the features to resample
    cFeats = { ...
        'inten',    {'env', 'peakEnv', 'peakRate'}, @sum; ...
        'pitch',    {'rF0', 'voicing', 'phrase', 'accent'}, @max; ...
        'pitch',    {'voicing'}, @sum; ...
        'word',     {'nSyll'}, @sum; ...
        'word',     {'blick', 'wordFreq'}, @mean; ...
        };
    sFeats = struct('tableName', cFeats(:,1), 'varNames', cFeats(:,2), 'readout', cFeats(:,3));
    ops.rsVars = sFeats;
    
    % Specify the responses to resample
    sResp = struct;
    sResp.tableName = 'spikeTime';
    sResp.varNames = se.GetTable('spikeTime').Properties.VariableNames;
    sResp.readout = 'rate';
    dt = [0.12 -0.12];
    
    % Resampling
    ce = NP.CodingExplorer;
    ce.userData = se.userData;
    ce.userData.ops = ops;
    
    fTbs = cell(size(sEvents));
    rTbs = cell(size(sEvents));
    for e = 1 : numel(sEvents)
        % Set the current grouping
        ops.eventVar = sEvents(e);
        en = ops.eventVar.varName;
        
        % Resample features and responses
        ops.rsVars = sFeats;
        ops.eventVar.readout = 0;
        fTbs{e} = NP.SE.ResampleEventFeatures(se, ops);
        
        ops.rsVars = sResp;
        ops.eventVar.readout = dt(e);
        rTbs{e} = NP.SE.ResampleEventFeatures(se, ops);
    end
    fTb = cat(1, fTbs{:});
    rTb = cat(1, rTbs{:});
    fTb.Properties.RowNames = evtNames;
    rTb.Properties.RowNames = evtNames;
    
    ce.SetTable("feat", fTb, 'timeSeries');
    ce.SetTable("resp", rTb, 'timeSeries');
    
    save(cePath, 'ce');
    fprintf("Saved the ce file at\n%s\n\n", cePath);
end
