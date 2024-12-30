%% 

figDir = "C:\chang_lab\project_np\misc\20240425_simon_keynote";

%% Configure raster plots

% Load unit cache
uid = [410100249 410100245; 460100005 410100463; 440300575 410100450];
uData = NP.Unit.LoadUnitCache(uid(:), 'DataSource', 'm2');
uSE = LMV.SE.UnitCache2SE(uData, 'StimIdList', LMV.Param.stimIdList4);

% Configure plots
[nr, nc] = size(uid);
uInd = reshape(1:numel(uid), nr, nc);
panels1 = repmat({'raster'}, [nr nc]);
% panels1 = repmat({'rate_map'}, [nr nc]);
panels1args = num2cell(num2cell(uInd));

args = { ...
    ["cue1", "stim", "cue3", "prod"], ...
    'Colors', LMV.Param.GetTaskPhaseColors(["cue1", "stim", "cue3", "prod"]), ...
    'YRange', [-2 0] ...
    };
panels2 = repmat({@NP.TaskBaseClass.PlotEventWindows}, size(panels1));
panels2args = repmat({args}, size(panels1));

bbArgs = {"stimText", 'Label', false, 'Color', [0 0 0 .3]};
panels3 = repmat({'block_boundary'}, size(panels1));
panels3args = repmat({bbArgs}, size(panels1));

tWin = [-0.3 6];

%% Plot rasters

f = MPlot.Figure(9517); clf
f.Position(3:4) = [1000 700];
NP.TaskBaseClass.PlotSE(uSE, panels3, 'PanelArgs', panels3args, 'TimeWindow', tWin);
NP.TaskBaseClass.PlotSE(uSE, panels2, 'PanelArgs', panels2args, 'TimeWindow', tWin);
NP.TaskBaseClass.PlotSE(uSE, panels1, 'PanelArgs', panels1args, 'TimeWindow', tWin);
exportgraphics(f, fullfile(figDir, "example_unit_rasters.png"));
exportgraphics(f, fullfile(figDir, "example_unit_rasters.pdf"));

%% Plot classification accuracy of example units in a heatmap

% Load classification models
mdlDir = fullfile(NP.Data.GetAnalysisRoot, 'sent_resp', 'sen4', "computed_mdls");
clusTb = NP.Fit.LoadModels([], mdlDir);
clusTb = cat(1, clusTb{:});

%% 

uid = uid';
uid = uid(:);
[~, I] = MMath.SortLike(clusTb.clusId, uid);
subTb = clusTb(I,:);

%% 

phaseNames = ["stim", "delay", "init", "prod"];
R = NaN(numel(uid), numel(phaseNames));
for i = 1 : numel(uid)
    for j = 1 : numel(phaseNames)
        pn = phaseNames(j);
        k = subTb.clusId==uid(i);
        mdl = subTb.(pn){k};
        if isempty(mdl)
            continue
        end
        R(i,j) = 1 - mdl.loss;
    end
end
R(R<=0.25) = NaN;
R = round(R, 2);

f = MPlot.Figure(5803); clf
h = heatmap(phaseNames, "Unit "+(1:numel(uid)), R);
h.MissingDataColor = 'white';
h.MissingDataLabel = 'N.S.';
% h.ColorbarVisible = 'off';
h.ColorLimits(1) = 0.25;
h.FontSize = 12;
MPlot.Paperize(f, 'ColumnsWide', 0.7, 'ColumnsHigh', 0.5);
exportgraphics(f, fullfile(figDir, "example_r.png"));


return
%% Screen units for specific response type

% Load responsiveness
rTest = LMV.Resp.LoadPhaseResponseTest();

% Find units that are stim-only
sigTb = rTest.sigTb;
sigTb{:,:} = sigTb{:,:} > 0;
isStimOnly = sigTb.stim & sum(sigTb{:,:},2) == 1;

% 
isRegion = rTest.clusTb.region == "mPrCG";

isUnit = isRegion & isStimOnly;
uid = rTest.clusTb.clusId(isUnit);
rStim = max(rTest.clusTb.stimResp(isUnit,:), [], 2);

uTb = table;
uTb.clusId = uid;
uTb.rStim = rStim;

uTb = sortrows(uTb, 'rStim', 'descend');

%% 

f = MPlot.Figure(9517); clf
LMV.Overview.SessionFromCache(uTb.clusId((1:12)+12), 'DataSource', 'm2', 'StimIdList', LMV.Param.stimIdList4);

