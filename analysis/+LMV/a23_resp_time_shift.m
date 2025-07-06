%% Case studies of bridge responses during delay

mdlName = "smooth_lm";
mdlDir = LMV.Data.GetAnalysisDir("linker", mdlName);
anaDir = fullfile(mdlDir, "resp_time");

%% 

% Load linker info table
clusTb = LMV.Linker.LoadClusTb(mdlName);
mirTb = clusTb(clusTb.hcGroup=="mirror", :);
fbTb = clusTb(clusTb.hcGroup=="feedback", :);

% Load linked positions
load(fullfile(mdlDir, "linked_positions.mat"), "posTb");

% Load necessary recordings
recIds = unique([mirTb.recId; fbTb.recId]);
sePaths = fullfile(LMV.Data.GetAnalysisDir, "data", "se_m11", recIds+"_se.mat");
seArray = NP.SE.LoadSession(sePaths);

%% Compute and test response time shift

warning off
mirStimTb = LMV.Linker.TS.ComputeRespTimeShift(mirTb, posTb, seArray, "stim");
mirProdTb = LMV.Linker.TS.ComputeRespTimeShift(mirTb, posTb, seArray, "prod");
fbStimTb = LMV.Linker.TS.ComputeRespTimeShift(fbTb, posTb, seArray, "stim");
fbProdTb = LMV.Linker.TS.ComputeRespTimeShift(fbTb, posTb, seArray, "prod");
warning on

%% 

titles = ["mirror stim", "mirror prod", "feedback stim", "feedback prod"];
lkTbs = {mirStimTb, mirProdTb, fbStimTb, fbProdTb};
posTbs = cell(size(lkTbs));
for i = 1 : numel(lkTbs)
    cTb = lkTbs{i};
    pTbWidths = cellfun(@width, cTb.posTb);
    pTb = cat(1, cTb.posTb{pTbWidths==max(pTbWidths)});
    
    m = pTb.numResp > 3 & pTb.score > LMV.Linker.scoreTh & abs(pTb.diff) < 0.1;
    disp(sum(m))
    pTb = pTb(m,:);
    
    posTbs{i} = pTb;
end

%% Plot time shift effect

f = MPlot.Figure(2); clf
rowDist = [3 2 2];
nCols = numel(titles);
tl = tiledlayout(sum(rowDist), nCols);
tl.Padding = "compact";

for i = 1 : numel(titles)
    % Plot scatter
    ntArgs = MPlot.FindTileInd(rowDist, nCols, 1, i);
    ax = nexttile(ntArgs{:});
    LMV.Linker.TS.PlotTimeShiftEffectScatter(posTbs{i}, "Axes", ax, "XLim", [-100 100]);
    ax.Title.String = titles(i);
    
    % Plot time shift histogram
    ntArgs = MPlot.FindTileInd(rowDist, nCols, 2, i);
    ax = nexttile(ntArgs{:});
    LMV.Linker.TS.PlotTimeShiftHistogram(posTbs{i}, "Axes", ax, "Edges", -100:20:100);

    % Plot scatter histogram
    ntArgs = MPlot.FindTileInd(rowDist, nCols, 3, i);
    ax = nexttile(ntArgs{:});
    LMV.Linker.TS.PlotTimeShiftEffectHistogram(posTbs{i}, "Axes", ax);
end

MPlot.Paperize(f, 1.3, .7);
exportgraphics(f, fullfile(anaDir, "resp_time_shift_stats_100ms.png"));

%% Write stats to file

% Write stats to file
fid = fopen(fullfile(anaDir, "resp_time_shift_stats.txt"), 'w');

% Write header
fprintf(fid, "Response Time Shift Statistics\n");
fprintf(fid, "==============================\n\n");

% Write stats for each condition
for i = 1 : numel(titles)
    fprintf(fid, "%s\n", titles(i));
    fprintf(fid, "---------------------------\n");

    % Compute effect size stats
    r = posTbs{i}.diffSize;
    rQt = prctile(r, [25 50 75]);
    p = signrank(r);
    
    % Compute mean time shifts with r < -0.2
    t = posTbs{i}.diff * 1e3;
    tQt = prctile(t(r < -0.2), [25 50 75]);

    % Write stats
    fprintf(fid, "N = %d\n", numel(r));
    fprintf(fid, "Effect size: median = %.2f, IQR = %.2f to %.2f\n", rQt(2), rQt(1), rQt(3));
    fprintf(fid, "Wilcoxon signed rank test on effect size p = %.3f\n", p);
    fprintf(fid, "Time shift (r < -0.2): median = %.2f ms, IQR = %.2f to %.2f\n", tQt(2), tQt(1), tQt(3));
    fprintf(fid, "\n");
end

fclose(fid);

%% Plot time shift combo

for i = 1 : numel(titles)
    f = MPlot.Figure(10+i); clf
    axs = LMV.Linker.TS.PlotTimeShiftCombo(posTbs{i});
    axs(1).Title.String = titles(i);
    MPlot.Paperize(f, .4, .4);
    s = lower(strrep(titles(i), " ", "_"));
    exportgraphics(f, fullfile(anaDir, sprintf("time_shift_combo_%s.png", s)));
    exportgraphics(f, fullfile(anaDir, sprintf("time_shift_combo_%s.pdf", s)));
end

%% Plot example responses

% Sample and plot for each condition
for i = 1 : 2
    % Get positions with effect size < -0.2
    m = posTbs{i}.diffSize < -0.2;
    if ~any(m)
        continue
    end
    
    % Create figure for this condition
    f = MPlot.Figure(20+i); clf
    tl = tiledlayout(2, 8);
    tl.Padding = "compact";
    
    % Get indices and sort by effect size magnitude (most negative first)
    idx = find(m);
    [~, sortIdx] = sort(posTbs{i}.diffSize(idx), 'ascend');
    idx = idx(sortIdx);
    
    % Take top 10 positions
    if numel(idx) > 16
        idx = idx(1:16);
    end
    
    % Plot response ladder
    pTb = posTbs{i}(idx,:);
    LMV.Linker.TS.PlotRespLadder(pTb);
    
    % Save figure
    MPlot.Paperize(f, 2, .9);
    s = lower(strrep(titles(i), " ", "_"));
    % exportgraphics(f, fullfile(anaDir, sprintf("example_responses_%s.png", s)));
    % exportgraphics(f, fullfile(anaDir, sprintf("example_responses_%s.pdf", s)), ContentType="vector");
end

return
%% 

titles = ["Mirror listening", "Mirror speaking"];

f = MPlot.Figure(7); clf
rowDist = [3 2];
nCols = numel(titles);
tl = tiledlayout(sum(rowDist), nCols);
tl.Padding = "compact";

for i = 1 : numel(titles)
    % Plot scatter
    ntArgs = MPlot.FindTileInd(rowDist, nCols, 1, i);
    ax = nexttile(ntArgs{:});
    LMV.Linker.TS.PlotIntervalRelations(posTbs{i}, "Axes", ax);
    ax.Title.String = titles(i);
    
    % Plot histogram
    ntArgs = MPlot.FindTileInd(rowDist, nCols, 2, i);
    ax = nexttile(ntArgs{:});
    LMV.Linker.TS.PlotIntervalCorrHistogram(posTbs{i}, "Axes", ax);
end
MPlot.Paperize(f, 1, .6);
exportgraphics(f, fullfile(anaDir, "interval_dependency_100ms.png"));
