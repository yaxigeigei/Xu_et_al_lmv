%% This is the standard script for overview of LMV task

% Load se
[se, sePath] = NP.SE.LoadSession([], 'UserFunc', @(x) x.RemoveTable('LFP', 'niTime'));
if isempty(se)
    return
end

% Create overview folder
seDir = fileparts(string(sePath));
recId = NP.SE.GetID(se);
ovDir = fullfile(seDir, recId);
if ~exist(ovDir, 'dir')
    mkdir(ovDir);
end

% Remove tiny clusters (for auto clustering)
clusTb = NP.Unit.GetClusTb(se);
NP.Unit.RemoveUnits(se, clusTb.numSpikes < 100);

% Standard enrichment
ops = NP.Param.Enrich;
ops.isFiltSpeech = true;
ops.isMel = true;
ops.isSpkRate = true;
ops.isSpkSpan = true;
% ops.isPitch = true;
se = LMV.SE.Transform(se, "enrich", ops);

% Exclude LMV BR trials
tv = se.GetTable("taskValue");
isRm = startsWith(tv.stimId, "BR_") | ismissing(tv.stimId);
se.RemoveEpochs(isRm);

%% Save plots and audio clips for inspecting each trial

% Create figure folder
trialsDir = fullfile(ovDir, 'lmv_trials_prod');
if ~exist(trialsDir, 'dir')
    mkdir(trialsDir);
end

f = MPlot.Figure(124);
f.WindowState = 'maximized';
pause(0.5);

tt = se.GetTable('taskTime');

for tr = 1 : se.numEpochs
    % Time window around production
    tWin = [tt.prod{tr}(1).T.tmin tt.prod{tr}(end).T.tmax] + [-1 1]*0.1;
    if all(isnan(tWin))
        continue
    end
    
    % Plot trial
    LMV.Overview.OneTrial(se, tr, tWin);
    figName = sprintf("%s_trial%i_%.3f-%.3fs.png", recId, tr, tWin(1), tWin(2));
    exportgraphics(f, fullfile(trialsDir, figName));
    
    % Save audio clip
    NP.Audio.WriteTrial(se, tr, tWin, 'FolderPath', trialsDir);
end

%% Compute LMV responses

% Extract task phase responses
s = LMV.Resp.ExtractPhaseResponses(se);

% Set dropout spike rates to NaN
rTb = s.respTb;
mTb = s.dropoutTb;
rTb{:,:}(mTb{:,:}) = NaN;

% Compute responsiveness
rTestTb = NP.Resp.PhaseSignrank(rTb, s.phaseTb, s.stimIdTb(:,1:4));
rSigTb = LMV.Resp.GetSigTable(rTestTb);

% Compute sentence selectivity
phaseNames = rSigTb.Properties.VariableNames;
sTestCell = NP.Resp.StimSelectivityKW(s.respTb, s.phaseTb(:,phaseNames), s.stimIdTb(:,1:4), 'Mask', rSigTb{:,:});
pSelectSig = cellfun(@(x) x.p, sTestCell);
selectSig = sum(pSelectSig < permute([0.05 0.01 0.001], [1 3 2]), 3);
sSigTb = array2table(selectSig, 'VariableNames', phaseNames);

%% Plot the numbers of responsive and selective units

f = MPlot.Figure(159); clf
NP.Resp.PlotNumberOfSigUnits({rSigTb, sSigTb});
legend(["responsive", "selective"], 'Location', 'eastoutside');
title(sprintf("%s, %i total responsive units", recId, sum(any(rSigTb{:,:},2))), 'Interpreter', 'none');
MPlot.Axes(gca);
MPlot.Paperize(f, 'ColumnsWide', 1.2, 'ColumnsHigh', 0.5);
exportgraphics(f, fullfile(ovDir, recId+"_lmv_resp_summary.png"));

%% Prepare data for rasters

% Remove inactive units
isRm = ~any(rSigTb{:,:}, 2);
seResp = se.Duplicate;
NP.Unit.RemoveUnits(seResp, isRm);

% Time morphing
seMorph = LMV.SE.Transform(seResp, "morph");

% Compute sentence PETHs
[ce, senTb] = LMV.SE.ComputeSentencePETH(seMorph);

%% Plot rasters for all units

rasterDir = fullfile(ovDir, 'lmv_rasters');
if ~exist(rasterDir, 'dir')
    mkdir(rasterDir);
end

upp = 15;
uIdPages = NP.UnitPlot.MakeArrayID(ce.clusTb.clusId, upp);
for p = 1 : size(uIdPages,2)
    f = MPlot.Figure(123);
    LMV.Overview.Sentences(senTb, 'UnitId', uIdPages(:,p), 'Page', p, 'Folder', rasterDir);
end

%% Plot rasters of selected units across 4 sentences

[~, cid] = LMV.Param.GetSelectedClusId(NP.SE.GetID(se));
if isempty(cid)
    return
end

upp = 15;
uIdPages = NP.UnitPlot.MakeArrayID(cid, upp);
for p = 1 : size(uIdPages,2)
    f = MPlot.Figure(123);
    LMV.Overview.Sentences(senTb, 'UnitId', uIdPages(:,p), 'Page', p, 'Folder', ovDir);
end

%% Plot rasters of selected units across 12 sentences

[~, cid] = LMV.Param.GetSelectedClusId(NP.SE.GetID(se));
if isempty(cid)
    return
end

rasterDir = fullfile(ovDir, 'lmv_rasters_12');
if ~exist(rasterDir, 'dir')
    mkdir(rasterDir);
end

upp = 4;
uIdPages = NP.UnitPlot.MakeArrayID(cid, upp);
for p = 1 : size(uIdPages,2)
    f = MPlot.Figure(123); clf
    LMV.Overview.Sentences(senTb, 'Features', {'phone'}, 'StimIdList', LMV.Param.stimIdList12, 'UnitIds', uIdPages(:,p), 'Page', p, 'Folder', rasterDir);
end
