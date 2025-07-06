%% Export unit PETHs as Pandas dataframes

ceVer = "ce_m2_ex3_sentence-avg";
outDir = LMV.Data.GetAnalysisDir("data", strrep(ceVer, "ce_", "pkl_"));

%% Load data

% Combined M2 sentence PETHs
load(fullfile(LMV.Data.GetAnalysisDir, "data", ceVer+".mat"), 'ce');

% Load responsiveness
rTest = LMV.Resp.LoadPhaseResponseTest();
clusTb = [rTest.clusTb(:,["clusId", "region", "recId"]) rTest.sigTb];
[~, I] = MMath.SortLike(clusTb.clusId, ce.clusTb.clusId);
clusTb = clusTb(I,:);

% Load linker info
lkTb = LMV.Linker.LoadClusTb("smooth_lm");
[~, I] = MMath.SortLike(clusTb.clusId, lkTb.clusId, false);
clusTb.hcGroup(I) = string(lkTb.hcGroup);
clusTb.hcGroup(ismissing(clusTb.hcGroup)) = "other";

%% Save dataframes

[respTb, semTb, ttTb, tvTb] = ce.GetTable("resp", "sem", "taskTime", "taskValue");

taskTb = [ ...
    tvTb(:, ["stimId", "stimText"]), ...
    ttTb(:, ["cue1On", "cue1Off", "stimMatchOn", "stimMatchOff", "cue3On", "cue3Off", "prodMatchOn", "prodMatchOff"]), ...
    ];

clusDf = py.pandas.DataFrame(clusTb);
respDf = py.pandas.DataFrame(respTb);
semDf = py.pandas.DataFrame(semTb);
taskDf = py.pandas.DataFrame(taskTb);

clusDf.to_pickle(fullfile(outDir, "clus.pkl"));
respDf.to_pickle(fullfile(outDir, "resp.pkl"));
semDf.to_pickle(fullfile(outDir, "sem.pkl"));
taskDf.to_pickle(fullfile(outDir, "task.pkl"));

return
%% Check if the input was extracted correctly

f = MPlot.Figure(84770); clf
tl = tiledlayout(7,2);
tl.Padding = "compact";

for i = 1 : size(Z,2)
    t = respTb.time{i};
    r = respTb.u410100245{i};
    
    ax = nexttile;
    plot(t, r);
    hold(ax, 'on');
    ax.YLim = [0 100];
    
    colNames = {'cue1', 'stim', 'prod', 'cue3'};
    cc = LMV.Param.GetTaskPhaseColors(colNames);
    NP.TaskBaseClass.PlotEventWindows(ax, ce, i, t([1 end])', colNames, 'YRange', ax.YLim, 'Colors', cc, 'Text', false);
    
    ax.Title.String = ce.GetTable("taskValue").stimText(i);
end

%% Baseline subtraction (not used)

[respTb, semTb] = ce.GetTable("resp", "sem");

respTb2 = respTb;
for i = 1 : height(respTb)
    t = respTb.time{i};
    isBase = t < 0 & t > -LMV.Param.respBaselineDur;
    for j = 2 : width(respTb)
        r = respTb.(j){i};
        respTb2.(j){i} = r - mean(r(isBase), 'omitmissing');
    end
end

