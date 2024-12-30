%% 

anaDir = LMV.Data.GetAnalysisDir('feat_qc');

%% Read data

mdlNames = {'F01_indep_Haskins_loss_90_filter_fix_bn_False_0_setting2', 'hprc_train_no_m1f2_h2tv_gru_tanh_nogan'};
aa = cell(size(mdlNames));
for i = 1 : numel(mdlNames)
    % Find AAI outputs
    articPattern = fullfile(anaDir, 'preproc_output', 'TIMIT_aai', 'artics', mdlNames{i}, '*_*.npy');
    articSearch = MBrowse.Dir2Table(articPattern);
    
    % Remove bad stim
    [~, stimIds] = cellfun(@fileparts, articSearch.name, 'Uni', false);
    badIds = {'mdwh0_si1925'};
    isBad = ismember(stimIds, badIds);
    articSearch(isBad,:) = [];
    stimIds(isBad) = [];
    
    % Read outputs
    aa{i} = cellfun(@readNPY, fullfile(articSearch.folder, articSearch.name), 'Uni', false);
end

% Read labels and waveforms
corpusDir = fullfile(NP.Data.GetProjectRoot, "code", "tasks", "LMV", "TIMIT");
tg = MLing.ReadTimitFeatures(corpusDir, stimIds);
tge = NP.TGEvent(tg);
[w, t] = MLing.ReadTimitWaveform(corpusDir, stimIds);

%% Construct se

se = MSessionExplorer();

% Make taskValue table
tv = table;
tv.stimIdx = (1:numel(stimIds))';
tv.stimId = stimIds;
se.SetTable('taskValue', tv, 'eventValues');

% Make taskTime table
tt = table;
tt.stimGT = tge;
tt.stimOn = tge.GetTfield('tmin');
tt.stimOff = tge.GetTfield('tmax');
se.SetTable('taskTime', tt, 'eventTimes');

% Make ni table
ni = table;
ni.time = t;
ni.speaker1 = w;
se.SetTable('ni', ni, 'timeSeries');

% Make artic table
tOn = zeros(size(tg));
tOff = cellfun(@(x) x(end), t);
isCombineEpochs = false;
for i = 1 : numel(aa)
    a = aa{i};
    switch size(a{1},2)
        case 18
            % For F01_indep_Haskins_loss_90_filter_fix_bn_False_0_setting2
            varNames = {'tt_x' 'tt_y' 'td_x' 'td_y' 'tb_x' 'tb_y' 'li_x' 'li_y' 'ul_x' 'ul_y' 'll_x' 'll_y' 'la' 'pro' 'ttcl' 'tbcl' 'v_x' 'v_y'};
        case 9
            % For hprc_train_no_m1f2_h2tv_gru_tanh_nogan
            varNames = {'la', 'lp', 'ja', 'ttcl', 'ttcd', 'tmcl', 'tmcd', 'trcl', 'trcd'};
        otherwise
            varNames = "akt" + (1:size(a{1},2));
    end
    artic = NP.Preproc.MakeTimeseriesTable(tOn, tOff, a, isCombineEpochs, varNames);
    se.SetTable("artic"+i, artic, 'timeSeries');
end

% Set dummy reference times
rt = [0; cumsum(tOff(1:end-1)+0.01)];
se.SetReferenceTime(rt);

sePath = fullfile(anaDir, 'timit_corpus_se.mat');
save(sePath, 'se');

return
%% Prepare inputs

% Vectorize se
seVec = se.Duplicate;
seVec.SliceSession(0, 'absolute');

% Resample artic
artic = seVec.GetTable('artic');
tEdges = artic.time{1}(1) : 0.01 : artic.time{1}(end);
artic = seVec.ResampleTimeSeries('artic', tEdges);
X = cell2mat(artic{:,2:end});

% Get phone events
tge = seVec.GetColumn('taskTime', 'stimGT');
phe = Cut(Cut(tge{1}));
[phList, phListId, phIds] = unique(phe.GetParentLabel);

% Create labels
Y = zeros(size(artic.time{1}));
for i = 1 : numel(phListId)
    isPh = phIds == phListId(i);
    m = phe(isPh).MaskTimestamps(artic.time{1});
    Y(m) = phListId(i);
end

%% Save inputs

s = struct;
s.artic = X;
s.phone_id = Y;
s.phone_list = cellstr(phList);
s.phone_list_id = phListId;

save(fullfile(anaDir, 'timit_corpus_feat.mat'), '-struct', 's');
