%% Compute and cache unit PETHs without first 3 trials

anaDir = LMV.Data.GetAnalysisDir('pop_dynamics');
srcTb = LMV.Data.FindSource([]);

%% Compute M2 sentence PETHs and average features

% Get file paths
m2Dir = fullfile(LMV.Data.GetAnalysisDir, 'data', 'se_m2');
m2epPaths = fullfile(m2Dir, strrep(srcTb.name, '_se.mat', '_se_m2_ep.mat'));

ceVer = "ce_m2_ex3_sentence-avg";
ceDir = LMV.Data.GetAnalysisDir('data', ceVer);
cePaths = fullfile(ceDir, strrep(srcTb.name, '_se.mat', '_ce.mat'));

for i = 1 : height(srcTb)
    % Check computed
    if isfile(cePaths{i})
        fprintf("\nThe M2 ce has been computed and saved at '%s'\n", cePaths{i});
        continue
    end
    
    % Load epoch-resliced se
    se = NP.SE.LoadSession(m2epPaths{i});
    
    % Remove the first 3 trials
    se.RemoveEpochs(1:3);
    
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
    
    % Add task phase events
    phases = struct;
    phases.atten = {'cue1On', 'stimOn'};
    phases.delay = {'stimOff', 'cue3On'};
    phases.init = {'cue3On', 'prodOn'};
    NP.TaskBaseClass.AddEventObjects(se, phases, 'taskTime');
    
    % Split and group trials of the same sentences
    % (No need to use LMV.SE.SplitBySentence since it's already resliced and cleaned)
    senTb = NP.TaskBaseClass.SplitBySentence(se);
    
    for s = 1 : height(senTb)
        % Extract and resample features and responses
        fprintf("\n%s\n", senTb.stimText(s));
        seSen = senTb.se(s);
        tt = seSen.GetTable('taskTime');
        
        ops = NP.Param.Resample;
        ops.rsWin = [tt.trialOn(1) tt.matchOff(1)] + [-1 1]*0.5;
        ops.rsBinSize = 0.01;
        ops.rsArgs = {};
        
        sCell = { ...
            'taskTime', {'cue1On', 'cue1Off', 'cue3On', 'cue3Off', 'stimOn', 'prodOn'}, []; ...
            'taskTime', {'cue1', 'atten', 'stim', 'delay', 'cue3', 'init', 'prod'}, 'span'; ...
            'mel',      {'mic', 'speaker1'}, []; ...
            'inten',    {'env', 'peakEnv', 'peakRate'}, []; ...
            'pitch',    {'F0', 'rF0', 'drF0'}, []; ...
            'artic',    {'ja', 'la', 'pro', 'ttcc', 'tbcc', 'tdcc', 'ttcl', 'tbcl', 'v_x', 'v_y'}, []; ...
            'phone',    phEvents(:,1)', 'span'; ...
            'taskValue', {'stimNumId'}, []; ...
            };
        ops.rsVars = struct('tableName', sCell(:,1), 'varNames', sCell(:,2), 'readout', sCell(:,3));
        
        fTb = NP.SE.ResampleFeatures(seSen, ops);
        rTb = NP.SE.ResampleResponses(seSen, ops);
        
        % Compute trial average
        mfTb = NP.SE.MeanTimeseries(fTb);
        mrTb = NP.SE.MeanTimeseries(rTb);
        
        % Construct ce
        ce = NP.CodingExplorer();
        
        ce.SetTable('feat', mfTb(1,:), 'timeSeries');
        
        ce.SetTable('resp', mrTb(1,:), 'timeSeries');
        ce.SetTable('sd', mrTb(2,:), 'timeSeries');
        ce.SetTable('sem', mrTb(3,:), 'timeSeries');
        
        LMV.SE.AddTemplateTrial(ce, seSen);
        ce.SetColumn('taskValue', 'numTrial', senTb.numTrial(s));
        senTb.ce(s) = ce;
    end
    
    % Concatenate ce of all sentences
    ce = Merge(senTb.ce);
    ce.userData = se.userData;
    ce.userData.rsOps = ops;
    
    save(cePaths{i}, 'ce');
end

%% Combine M2 sentence PETHs with all units

% Load PETHs
cePaths = fullfile(ceDir, srcTb.recId+"_ce.mat");
ceArray = NP.CE.LoadSession(cePaths);

% Keep a standard set of sentences
stimIdList = LMV.Param.stimIdList14;

% Use NP46_B1 as the template ce
ce0 = ceArray(12).Duplicate({'resp', 'sd', 'sem', 'taskTime', 'taskValue'});
tv = ce0.GetTable("taskValue");
[~, I] = MMath.SortLike(tv.stimId, stimIdList);
ce0.SortEpochs(I);

ce14Array = ceArray; % just for preallocation

for i = 1 : numel(ceArray)
    ceSrc = ceArray(i);
    disp(NP.SE.GetID(ceSrc));
    
    % Create empty response tables
    ce = ce0.Duplicate;
    emptyTb = repmat(ce.GetTable("resp").time(1), ce.numEpochs, ceSrc.numResp+1);
    emptyTb(:,2:end) = cellfun(@(x) NaN(size(x)), emptyTb(:,2:end), 'Uni', false);
    emptyTb = cell2table(emptyTb, 'VariableNames', ["time" string(ceSrc.respNames)]);
    
    % Fill in sentence responses
    tbNames = ["resp", "sd", "sem"];
    tv = ceSrc.GetTable("taskValue");
    for j = 1 : numel(tbNames)
        tb0 = emptyTb;
        tb1 = ceSrc.GetTable(tbNames(j));
        for k = 1 : height(tb0)
            isSen = tv.stimId==stimIdList(k);
            if ~any(isSen)
                if j == 1
                    fprintf("Missing '%s'\n", stimIdList(k));
                end
                continue
            end
            tEdges = MMath.BinCenters2Edges(tb0.time{k});
            tb0(k,:) = ceSrc.ResampleTimeSeries(tbNames(j), tEdges, isSen, [], 'Method', 'linear', 'Extrapolation', 'nearest');
        end
        ce.SetTable(tbNames(j), tb0);
    end
    
    ce14Array(i) = ce;
end

ce = ce14Array.CatUnits;
ce.clusTb = cat(1, ceArray.clusTb);

save(fullfile(anaDir, ceVer+".mat"), 'ce');

%{

NP35_B2
Missing 'mbbr0_si2315'
NP43_B1
Missing 'mdlc2_si2244'
Missing 'fsjk1_si2285'
Missing 'mbbr0_si2315'
NP52_B1
Missing 'mbbr0_si2315'
NP52_B2
Missing 'fcaj0_si1479'

%}
