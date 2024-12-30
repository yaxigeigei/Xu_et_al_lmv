%% Resample features and spike counts for spike-triggered analyses

ceDir = LMV.Data.GetAnalysisDir("st_feat", "ce");

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
    sePath = fullfile(LMV.Data.GetAnalysisDir, "data", "se_m1", strrep(srcTb.name{i}, 'se', 'se_m1'));
    se = NP.SE.LoadSession(sePath);
    
    %% Add derived features
    
    % Remove spike rate table
    %   NP.SE.ResampleResponses below will resample from spike times with no smoothing
    se.RemoveTable('spikeRate');
    
    % Add sentence onsets
    LMV.SE.AddSentenceOnsets(se);
    
    % Add trial events
    LMV.SE.AddTrialEvents(se);
    
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
    
    % Initialize option struct
    ops = NP.Param.Resample;
    
    % Set the resmapling window and sampling rate
    rt = se.GetTable('taskTime').trialOn{1};
    ops.rsWin = [rt(1)-5 rt(end)+15];
    ops.rsBinSize = 0.02;
    
    % Specify the features to resample
    sCell = { ...
        'taskTime', {'cue1On', 'cue1Off', 'cue3On', 'cue3Off', 'stimOn', 'prodOn', 'sentOn', 'wordOn', 'syllOn'}, []; ...
        'taskTime', {'cue1', 'atten', 'stim', 'delay', 'cue3', 'init', 'prod', 'iti'}, 'span'; ...
        'mel',      {'mic', 'speaker1'}, []; ...
        'inten',    {'env', 'peakEnv', 'peakRate'}, []; ...
        'pitch',    {'rF0', 'drF0'}, []; ... % {'rF0', 'drF0', 'brF0', 'voicing', 'phrase', 'accent'}
        'artic',    {'ja', 'la', 'pro', 'ttcc', 'tbcc', 'tdcc', 'ttcl', 'tbcl', 'v_x', 'v_y'}, []; ...
%         'word',     {'nSyll', 'blick', 'wordFreq'}, 'event'; ...
        'phone',    phEvents(:,1)', 'span'; ...
        'trial',    cellstr(LMV.Param.stimIdList), 'span'; ...
        };
    ops.rsVars = struct('tableName', sCell(:,1), 'varNames', sCell(:,2), 'readout', sCell(:,3));
    
    % Resample features and responses
    fTb = NP.SE.ResampleFeatures(se, ops);
    rTb = NP.SE.ResampleResponses(se, ops);
    
    % Convert spike rates to spike counts
    rCell = rTb{:,2:end};
    rCell = cellfun(@(x) round(x*ops.rsBinSize), rCell, 'Uni', false);
    rTb{:,2:end} = rCell;
    
    
    % Construct and save ce
    ce = NP.CodingExplorer;
    ce.userData = se.userData;
    ce.userData.ops = ops;
    ce.SetTable('feat', fTb, 'timeSeries');
    ce.SetTable('resp', rTb, 'timeSeries');
    save(cePath, 'ce');
    fprintf("Saved the ce file at\n%s\n\n", cePath);
end

