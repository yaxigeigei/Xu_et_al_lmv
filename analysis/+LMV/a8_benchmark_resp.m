%% Manually score responsiveness and use it to benchmark different responsiveness methods

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'benchmark');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end

srcTb = NP.Data.FindSource('lmv');

%% Save rasters

for i = 1 : height(srcTb)
    % Get selected units
    [~, uIds] = LMV.Param.GetSelectedClusId(srcTb.recId{i});
    if numel(uIds) < 5
        continue
    end
    
    nMax = 15;
    uInd = 1 : min(nMax, numel(uIds));
    uId4plot = uIds(uInd);
    srcTb.uId4plot{i} = uId4plot';
    
    % Plot and save rasters
    f = MPlot.Figure(5); clf
    LMV.Overview.SentencesFromCache(uId4plot, 'm1', 'UnitLabels', "depth");
    figPath = fullfile(anaDir, sprintf("%s_rasters.png", srcTb.recId{i}));
    exportgraphics(f, figPath);
end

%% Save a template spreadsheet for scoring

scoreTb = table;
scoreTb.clusId = cat(1, srcTb.uId4plot{:});
scoreTb.recId = repelem(srcTb.recId, cellfun(@numel, srcTb.uId4plot));
taskPhases = ["atten", "stim", "delay", "init", "prod"];
scoreTb{:,taskPhases} = NaN(height(scoreTb), numel(taskPhases));

writetable(scoreTb, fullfile(anaDir, "score_resp_template.csv"));

%% Load test results and scored spreadsheets

% Read scored spreadsheet
scoreTb = readtable(fullfile(anaDir, "score_resp_mx.csv"));
taskPhases = ["atten", "stim", "delay", "init", "prod"];
scoreTb(any(isnan(scoreTb{:,taskPhases}),2), :) = [];

% Load t responsiveness
sT = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'ttest', 'computed_ttest_clusTb.mat'));
tTb = NP.Unit.AlignClusTb(scoreTb, sT.clusTb, false);

% Load t responsiveness using refined windows
sT2 = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'ttest2', 'computed_ttest_clusTb.mat'));
t2Tb = NP.Unit.AlignClusTb(scoreTb, sT2.clusTb, false);

% Load Zeta responsiveness
sZeta = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'zeta', 'computed_zeta_clusTb.mat'));
zTb = NP.Unit.AlignClusTb(scoreTb, sZeta.clusTb, false);

%% Determine significance with different methods

func = @(tb) { ...
    NP.Resp.GetSigTable(tb, 'SubMask', 1), ...
    NP.Resp.GetSigTable(tb, 'Correction', true), ...
    NP.Resp.GetSigTable(tb, 'Correction', false), ... 
    };

tTbs = func(tTb);
t2Tbs = func(t2Tb);
zTbs = func(zTb);

aTbs = tTbs;
for i = 1 : numel(tTbs)
    aTbs{i}{:,:} = tTbs{i}{:,:} > 0 | zTbs{i}{:,:} > 0;
end

%% Benchmark different test methods

testTbs = [tTbs zTbs aTbs];
testNames = ["t-test", "ZETA-test", "t-OR-z"] + [" 1", " 5" ," 5 nc"]';
testNames = testNames(:);

f = MPlot.Figure(990); clf
tl = tiledlayout(2,1);
tl.Padding = 'compact';
NP.Resp.PlotBenchmark(scoreTb, testTbs, taskPhases);
legend(nexttile(1), testNames, 'Location', 'northeastoutside');
MPlot.Paperize(f, 'ColumnsHigh', 0.7, 'ColumnsWide', 1.2);
exportgraphics(f, fullfile(anaDir, "benchmark_resp_tests.png"));

%% Benchmark different versions of t-tests

testTbs = [tTbs t2Tbs];
testNames = ["t-test", "t-test rw"] + [" 1", " 5" ," 5 nc"]';
testNames = testNames(:);

f = MPlot.Figure(991); clf
tl = tiledlayout(2,1);
tl.Padding = 'compact';
NP.Resp.PlotBenchmark(scoreTb, testTbs, taskPhases(2:end));
legend(nexttile(1), testNames, 'Location', 'northeastoutside');
MPlot.Paperize(f, 'ColumnsHigh', 0.7, 'ColumnsWide', 1.2);
exportgraphics(f, fullfile(anaDir, "benchmark_t-tests.png"));

return
%% 

f = MPlot.Figure(990); clf
tl = tiledlayout(numel(taskPhases), 3);
tl.Padding = 'compact';

for i = 1 : numel(taskPhases)
    pn = taskPhases(i);
    y = scoreTb.(pn) > 0;
    yt = tTb.(pn) > 0;
    yz = zTb.(pn) > 0;
    
    ax = nexttile;
    h = confusionmat(y, yt, 'Normalization', 'total-normalized');
    h.Title = "t: " + pn;
    
    ax = nexttile;
    h = confusionchart(y, yz, 'Normalization', 'total-normalized');
    h.Title = "z: " + pn;
    
    ax = nexttile;
    h = confusionchart(y, yt|yz, 'Normalization', 'total-normalized');
    h.Title = "t|z: " + pn;
end

