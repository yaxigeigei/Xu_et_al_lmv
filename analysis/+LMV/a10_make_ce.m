%% Prepare features and responses for fitting TRF

ceDir = LMV.Data.GetAnalysisDir('trf', 'ce');
srcTb = LMV.Data.FindSource([]);

%% 

for i = 1 : height(srcTb)
    %% Load and standardize se
    
    % Check file
    cePath = fullfile(ceDir, strrep(srcTb.name{i}, 'se', 'ce'));
    if exist(cePath, 'file')
        fprintf("Not overwriting the existing ce file at\n%s\n\n", srcTb.path{i});
        continue
    end
    
    % Load recording
    se = NP.SE.LoadSession(srcTb.path{i});
    
    %% Add derived features
    
    % Standard enrichment
    ops = NP.Param.Enrich;
    ops.isSpkRate = true;
    ops.isMel = true;
    ops.isPitch = true; % this also adds 'trial' table, no need to run separately
    ops.isArtic = true;
    se = LMV.SE.Transform(se, 'enrich', ops); % this only keeps LMV trials
    
    % Set sentence onsets
    LMV.SE.AddSentenceOnsets(se);
    
    % Vectorize se
    se.SliceSession(0, 'absolute');
    
    % Add phoneme events
    phones = [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.consonants]';
    phGroups = { ...
        'high',         NP.Phone.high; ...
        'mid',          NP.Phone.mid; ...
        'low',          NP.Phone.low; ...
        'front',        NP.Phone.front; ...
        'back',         NP.Phone.back; ...
        'rounded',      NP.Phone.rounded; ...
        'plosives',     NP.Phone.plosives; ...
        'fricatives',   NP.Phone.fricatives; ...
        'nasals',       NP.Phone.nasals; ...
        'approximants', NP.Phone.approximants; ...
        'labial',       NP.Phone.labial; ...
        'velar',        NP.Phone.velar; ...
        'coronal',      NP.Phone.coronal; ...
        'glottal',      NP.Phone.glottal; ...
        'dental',       NP.Phone.dental; ...
        };
    phEvents = [repmat(phones, [1 2]); phGroups];
    NP.Phone.AddPhoneEvents(se, phEvents);
    
    % Add word level features
    LMV.Word.AddWordEvents(se, "lmv_syllabification_lookup.csv");
    LMV.Word.AddWordTimeseriesTable(se);
    
    % Add task phase events
    phases = struct;
    phases.atten = {'cue1On', 'stimOn'};
    phases.delay = {'stimOff', 'cue3On'};
    phases.init = {'cue3On', 'prodOn'};
    phases.iti = {'prodOff', 'cue1On'};
    NP.TaskBaseClass.AddEventObjects(se, phases, 'taskTime');
    
    %% Resampling
    
    % Initialize option struct from the enrichment option struct
    ops = NP.Param.Resample(ops);
    
    % Set the resmapling window and sampling rate
    probe = 1; % use first probe
    ksOps = se.userData.ksMeta(probe).ops;
    apMeta = se.userData.apMeta(probe);
    tRange = [ksOps.tstart ksOps.tend] / apMeta.imSampRate; % use the sorted period as the time window
    ops.rsWin = tRange;
    ops.rsBinSize = 0.01;
    
    % Specify the features to resample
    sCell = { ...
        'taskTime', {'cue1On', 'cue1Off', 'cue3On', 'cue3Off', 'stimOn', 'prodOn', 'sentOn', 'wordOn', 'syllOn'}, []; ...
        'taskTime', {'cue1', 'atten', 'stim', 'delay', 'cue3', 'init', 'prod', 'iti'}, 'span'; ...
        'mel',      {'mic', 'speaker1'}, []; ...
        'inten',    {'env', 'peakEnv', 'peakRate'}, []; ...
        'pitch',    {'rF0', 'drF0', 'brF0', 'voicing', 'phrase', 'accent'}, []; ...
        'artic',    [LMV.TRF.GetFeatureSet('artic') "d_"+LMV.TRF.GetFeatureSet('artic')], []; ...
        'word',     {'nSyll', 'blick', 'wordFreq'}, 'event'; ...
        'phone',    phEvents(:,1)', []; ...
        'trial',    cellstr(LMV.Param.stimIdList), 'span'; ...
        };
    ops.rsVars = struct('tableName', sCell(:,1), 'varNames', sCell(:,2), 'readout', sCell(:,3));
    
    
    % Resample features and responses
    fTb = NP.SE.ResampleFeatures(se, ops);
    rTb = NP.SE.ResampleResponses(se, ops);
    
    
    % Construct and save ce
    ce = NP.CodingExplorer;
    ce.userData = se.userData;
    ce.userData.ops = ops;
    ce.SetTable('feat', fTb, 'timeSeries');
    ce.SetTable('resp', rTb, 'timeSeries');
    save(cePath, 'ce');
    fprintf("Saved the ce file at\n%s\n\n", cePath);
    
end
