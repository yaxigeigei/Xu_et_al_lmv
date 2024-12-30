%% 

ceDir = LMV.Data.GetAnalysisDir('coding', 'phone_ce');
srcTb = LMV.Data.FindSource([]);

%% 

for i = 1 : height(srcTb)
    %% Load se
    
    % Check file
    cePath = fullfile(ceDir, strrep(srcTb.name{i}, 'se', 'ce'));
    if exist(cePath, 'file')
        fprintf("Not overwriting the existing ce file at\n%s\n\n", srcTb.path{i});
        continue
    end
    
    % Load recording
    se = NP.SE.LoadSession(srcTb.path{i});
    
    %% Add features
    
    % Standard enrichment
    ops = NP.Param.Enrich;
    ops.isMel = true;
    ops.isPitch = true; % this also adds 'trial' table, no need to run separately
    ops.isArtic = true;
    se = LMV.SE.Transform(se, 'enrich', ops);
    
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
    
    % Add biphone events
    
    
    %% Resampling
    
    % Set resampling events
    tt = se.GetTable('taskTime');
    targets = ["stim", "prod"];
    for t = 1 : numel(targets)
        tn = targets(t);
        tge = tt.(tn){1};
        tge = Cut(Cut(tge));
        tge = sort(tge);
        for n = 1 : numel(tge)
            s = tge(n).parentLabel;
            s(s>='0' & s<='9') = []; % remove the trailing digit (stress)
            tge(n).parentLabel = s;
        end
        tt.(tn+"Ph"){1} = tge;
    end
    se.SetTable('taskTime', tt);
    
    % Initialize option struct from the enrichment option struct
    ops = NP.Param.Resample(ops);
    
    % Specify the grouping event
    evtNames = {'stimPh', 'prodPh'};
    sEvents = struct('tableName', 'taskTime', 'varName', evtNames, 'readout', 0);
    ops.eventVar = sEvents;
    
    % Specify the features to resample
    cFeats = { ...
        'inten',    {'env', 'peakEnv', 'peakRate'}, []; ...
        'pitch',    {'rF0', 'drF0', 'brF0', 'voicing', 'phrase', 'accent'}, []; ...
        'artic',    {'ja', 'la', 'pro', 'ttcc', 'tbcc', 'tdcc', 'ttcl', 'tbcl', 'v_x', 'v_y'}, []; ...
        'phone',    phEvents(:,1)', []; ...
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

return
%% 

for i = 7 %1 : height(srcTb)
    % Load se
    se = NP.SE.LoadSession(srcTb.path{i});
    
    % Standard enrichment
    ops = NP.Param.Enrich;
    ops.isMel = true;
    ops.isPitch = true; % this also adds 'trial' table, no need to run separately
    ops.isArtic = true;
    se = LMV.SE.Transform(se, 'enrich', ops); % this only keeps LMV trials
    
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
    
    % Compute phoneme average
    ops = NP.Param.Resample;
    
    ops.groups = [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.consonants];
    
    repFeats = struct;
    repFeats.MannerPlace = {'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', 'nasals', 'approximants', 'labial', 'velar', 'coronal', 'dental'};
    repFeats.Artic = {'ja', 'la', 'pro', 'ttcc', 'tbcc', 'tdcc', 'ttcl', 'tbcl', 'v_x', 'v_y'};
    repFeats.Acous = [{'env', 'peakEnv', 'peakRate'}, {'voicing', 'phrase', 'accent', 'drF0'}];
    ops.repFeats = repFeats;
    
    mn = 'prod';
    tge = se.GetTable('taskTime').(mn){1};
    tge = Cut(Cut(tge));
    tge = sort(tge);
    for n = 1 : numel(tge)
        s = tge(n).parentLabel;
        s(s>='0' & s<='9') = []; % remove the trailing digit (stress)
        tge(n).parentLabel = s;
    end
    lbs = tge.GetParentLabel;
    
    fTbs = cell(numel(ops.groups), 1);
    rTbs = cell(numel(ops.groups), 1);
    
    for p = 1 : numel(ops.groups)
        phName = ops.groups{p};
        phEvts = tge(lbs == phName);
        if isempty(phEvts)
            continue
        end
        tWins = [phEvts.GetTfield('tmin') phEvts.GetTfield('tmax')];
        
        % Build feature table
        tbCat = table;
        tbCat.time = tWins(:,1);
        tbCat.label = phEvts.GetParentLabel;
        
        vn = repFeats.Artic;
        tb = se.SliceTimeSeries('artic', {tWins}, [], vn);
        arr = cellfun(@mean, tb{:,2:end});
        arr = double(arr);
        arr(isnan(arr)) = 0;
        tbCat = [tbCat, array2table(arr, 'VariableNames', vn)];
        
        vn = {'env', 'peakEnv', 'peakRate'};
        tb = se.SliceTimeSeries('inten', {tWins}, [], vn);
        arr = cellfun(@mean, tb{:,2:end});
        arr = double(arr);
        arr(isnan(arr)) = 0;
        tbCat = [tbCat, array2table(arr, 'VariableNames', vn)];
        
        vn = {'voicing', 'phrase', 'accent', 'drF0'};
        tb = se.SliceTimeSeries('pitch', {tWins}, [], vn);
        arr = cellfun(@mean, tb{:,2:end});
        arr = double(arr);
        arr(isnan(arr)) = 0;
        tbCat = [tbCat, array2table(arr, 'VariableNames', vn)];
        
        vn = repFeats.MannerPlace;
        tb = se.SliceEventTimes('phone', {tWins}, [], vn);
        arr = cellfun(@(x) ~isnan(x(1)), tb{:,:});
        arr = double(arr);
        tbCat = [tbCat, array2table(arr, 'VariableNames', vn)];
        
        fTbs{p} = tbCat;
        
        % Build response table
        tbCat = table;
        tbCat.time = tWins(:,1);
        tbCat.label = phEvts.GetParentLabel;
        
        tb = se.SliceEventTimes('spikeTime', {tWins - 0.15});
        dur = repmat(phEvts.GetDuration, [1 width(tb)]);
        arr = cellfun(@(x,d) numel(x(~isnan(x)))/d, tb{:,:}, num2cell(dur));
        tbCat = [tbCat, array2table(arr, 'VariableNames', se.GetTable('spikeTime').Properties.VariableNames)];
        
        rTbs{p} = tbCat;
    end
    
    fTb = cat(1, fTbs{:});
    rTb = cat(1, rTbs{:});
    
    % 
    ce = NP.CodingExplorer;
    ce.userData = se.userData;
    ce.userData.ops = ops;
    ce.SetTable('feat', fTb, 'timeSeries');
    ce.SetTable('resp', rTb, 'timeSeries');
end

%% 

% Load unit responsiveness table
sTest = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'ttest', "computed_ttest_clusTb.mat"));

%% Standalone RSA

mn = 'prod';
disp(NP.SE.GetID(ce));

% Select responsive units
isRec = sTest.clusTb.recId == NP.SE.GetID(ce);
pval = sTest.clusTb{isRec,mn};
uMask = any(pval < 0.05, 2);

% Configure analysis
ops = ce.userData.ops;
featFam = struct;
featFam.MannerPlace = {'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', 'nasals', 'approximants', 'labial', 'velar', 'coronal', 'dental'};
featFam.Artic = {'ja', 'la', 'pro', 'ttcc', 'tbcc', 'tdcc', 'ttcl', 'tbcl', 'v_x', 'v_y'};
featFam.Acous = [{'env', 'peakEnv', 'peakRate'}, {'voicing', 'phrase', 'accent', 'drF0'}];
ops.featFam = featFam;
ops.unitInd = find(uMask);
ops.maskName = mn;
ops.nBoot = 0; % not computing pval


% Prepare input
mfTb = groupsummary(fTb, 'label', 'mean');
mfTb.GroupCount = [];
mfTb.Properties.VariableNames = erase(mfTb.Properties.VariableNames, 'mean_');
mfTb = movevars(mfTb, 'label', 'After', 'time');

mrTb = groupsummary(rTb, 'label', 'mean');
mrTb.GroupCount = [];
mrTb.Properties.VariableNames = erase(mrTb.Properties.VariableNames, 'mean_');
mrTb.label = [];

ph2plot = [NP.Phone.approximants, NP.Phone.velar, NP.Phone.postalveolar, NP.Phone.alveolar, NP.Phone.dental, NP.Phone.labial];
ph2plot = unique(ph2plot, 'stable');
ph2plot = flip(ph2plot);
ph2plot = {'P', 'B', 'M', 'F', 'V' ,'TH' ,'DH', 'T', 'D', 'N', 'S', 'Z', 'CH', 'JH', 'SH', 'ZH', 'NG', 'K', 'G', 'Y', 'L', 'R', 'W'};
ph2plot = {'OY', 'OW', 'AO', 'AA', 'AW', 'AY', 'AE', 'EH', 'EY', 'IY', 'IH', 'AH', 'UW', 'ER', 'UH'};
isPh = ismember(mfTb.label, ph2plot);
mfTb = mfTb(isPh,:);
mrTb = mrTb(isPh,:);

R = mrTb{:,2:end};
R = R(:,uMask);
FF = cellfun(@(x) mfTb{:,x}, struct2cell(ops.featFam)', 'Uni', false);
XX = [{R}, FF];

sRez = ops;
sRez.ce = ce;
sRez.t = mfTb.time;
sRez.stimIdx = mfTb.label;
sRez.repNames = ["Neural", fieldnames(ops.featFam)'];
sRez.input = XX;


% % Compute RDM
% sRez = NP.RSA.ComputeSim(sRez);

%% Plot RDMs with phoneme grouping

sRez = NP.RSA.ComputeRDMs(sRez);

RDM = sRez.RDM;
repNames = sRez.repNames;

% Use phonemes for grouping
% ph2plot = [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.approximants, NP.Phone.velar, NP.Phone.postalveolar, NP.Phone.alveolar, NP.Phone.dental, NP.Phone.labial];
% ph2plot = [NP.Phone.approximants, NP.Phone.velar, NP.Phone.postalveolar, NP.Phone.alveolar, NP.Phone.dental, NP.Phone.labial];
% ph2plot = unique(ph2plot, 'stable');

cond = categorical(MLing.ARPA2IPA(mfTb.label), MLing.ARPA2IPA(ph2plot), 'Ordinal', true);

% Plot RDMs
f = MPlot.Figure(12); clf
tl = tiledlayout('flow');
tl.Padding = 'compact';
for i = 1 : numel(RDM)
    ax = nexttile;
    NP.RSA.PlotRDM(RDM{i}, cond);
    ax.Title.String = repNames(i);
%     colorbar;
end
% MPlot.Paperize(f, 'ColumnsWide', 1.4, 'Aspect', 0.667);
% exportgraphics(f, fullfile(anaDir, sprintf("rdms_%s_%s_group-ph.png", NP.SE.GetID(ce), mn)));

%% 

f = MPlot.Figure(2); clf
NP.RSA.PlotPSM(sRez);
h = gca;
h.Title = strrep(NP.SE.GetID(ce), '_', '-') + " " + mn;
% h.Interpreter = 'none'; % this option is only available in 2023b
MPlot.Paperize(f, 'ColumnsWide', .6, 'Aspect', .85);
% exportgraphics(f, fullfile(anaDir, sprintf("psm_%s_%s.png", NP.SE.GetID(ce), mn)));
