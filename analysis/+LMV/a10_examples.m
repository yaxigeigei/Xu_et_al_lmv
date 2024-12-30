%% Sentence raster and TRF heatmaps of example units

mdlName = "combo3pm";
anaDir = LMV.Data.GetAnalysisDir('trf', mdlName);

%% Get example unit IDs

uid = LMV.Param.GetExampleClusId("selectivity");

% Sort units by depth
uData = NP.Unit.LoadUnitCache(uid);
d = cellfun(@(x) x.unitInfo.depth, uData);
[~, I] = sort(d, 'ascend');
uid = uid(I);

%% 

f = MPlot.Figure(9989); clf
LMV.Overview.SentencesFromCache(uid, 'Features', "phone");

%% 

f = MPlot.Figure(8989); clf
NP.TRF.PlotWeightsComboFromCache(uid, [], mdlName+"_"+["stim", "prod", "feedback"]);


