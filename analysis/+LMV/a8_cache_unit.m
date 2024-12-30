%% Cache basic data and results of each unit for plotting

srcTb = LMV.Data.FindSource([]);

morphType = "m1";

% We will use epoch-sliced M1 data for sentence alignment
seDir = fullfile(NP.Data.GetAnalysisRoot, "data", "se_"+morphType);
sePaths = fullfile(seDir, srcTb.recId + "_se_"+morphType+"_ep.mat");

% Cache location
cacheDir = fullfile(seDir, "unit_cache");

% Restore QC table
qcPaths = fullfile(LMV.Data.GetAnalysisDir, 'units', 'computed_qc', srcTb.recId+"_clusTb.mat");
qcData = arrayfun(@load, qcPaths);
clusTb = cat(1, qcData.clusTb);
% clusTb.clusId = clusTb.clusId + arrayfun(@(x) NP.Unit.GetBaseClusId(x), clusTb.recId);

%% Add unit info, rasters, and example trials to cache

NP.Unit.CreateUnitCache(sePaths, cacheDir);

%% Load unit properties

% Load task phase responsiveness
sResp = LMV.Resp.LoadPhaseResponseTest();
respTb = NP.Unit.AlignClusTb(clusTb, sResp.clusTb, true);

% Load sentence selectivity
sSelect = LMV.Resp.LoadSentenceSelectTest();
selectTb = NP.Unit.AlignClusTb(clusTb, sSelect.clusTb, true);

%% Load TRF models

sets = ["phone", "strf", "artic3"];
targets = LMV.TRF.targets;

trfTb = table;
for i = 1 : numel(sets)
    % Load model tables
    mdlNames = sets(i) + "_" + targets;
    mdlTb = LMV.TRF.LoadModels(mdlNames, srcTb.recId);
    
    % Append models
    mdlTb = cat(1, mdlTb{:});
    mdlTb = NP.Unit.AlignClusTb(clusTb, mdlTb, true);
    trfTb = [trfTb mdlTb];
end

%% 

s = struct;
% s.unitInfo = table2struct(clusTb);
% s.tResp = table2struct(respTb);
% s.kwSelect = table2struct(selectTb);
s.trf = table2struct(trfTb);

NP.Unit.AddData2UnitCache(cacheDir, clusTb.clusId, s);

return
%% Test plotting

[~, uId] = LMV.Param.GetSelectedClusId('NP41_B1');
NP.Unit.LoadUnitCache(uId);

%% 

MPlot.Figure(987); clf
LMV.Overview.SessionFromCache(uId(1:10), 'DataSource', morphType, 'StimIdList', LMV.Param.stimIdList12);

%%

MPlot.Figure(458); clf
LMV.Overview.SentencesFromCache(uId(1:15), 'DataSource', morphType);

%%

% MPlot.Figure(892); clf
% LMV.TRF.PlotWeightsComboFromCache(uId(1:15), 'm1', "phone_"+LMV.TRF.targets);
LMV.TRF.PlotWeightsComboFromCache(cid, 'm1', "phone_"+LMV.TRF.targets);

