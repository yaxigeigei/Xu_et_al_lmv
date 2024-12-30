%% Case studies of bridge responses during delay

mdlName = "smooth_lm";
figDir = fullfile(LMV.Data.GetAnalysisDir, "linker", mdlName);

%% Specify units and sentences

uu = cell(0);

uu{end+1} = struct('clusId', 410100155, 'startText', "have you got");
uu{end+1} = struct('clusId', 410100463, 'startText', "have you got");
uu{end+1} = struct('clusId', 440300574, 'startText', "have you got");

uu{end+1} = struct('clusId', 410100441, 'startText', "we've got");
uu{end+1} = struct('clusId', 450100165, 'startText', "we've got");

uu{end+1} = struct('clusId', 440300575, 'startText', "junior");
uu{end+1} = struct('clusId', 450100165, 'startText', "junior");
uu{end+1} = struct('clusId', 450100323, 'startText', "junior");

caseTb = struct2table([uu{:}]);

%% 

clusTb = LMV.Linker.LoadClusTb(mdlName);
isUnit = ismember(clusTb.clusId, caseTb.clusId);
isBridge = clusTb.hcGroup=="bridge";
cTb = clusTb(isUnit & isBridge,:);

% isUnit = ismember(uu.clusId, cTb.clusId);
% uu = uu(isUnit,:);

%% Load se

recIdList = unique(cTb.recId);
srcTb = LMV.Data.FindSource(recIdList);
seArray = NP.SE.LoadSession(srcTb.path, 'UserFunc', @(x) x.RemoveTable('ni', 'LFP'));

arrayfun(@NP.SE.SetMorphTimes, seArray);
[~, senTbs] = arrayfun(@(x) LMV.SE.Transform(x, 'sentence'), seArray, 'Uni', false);

%% Extract sentence se for each case

for i = 1 : height(caseTb)
    % Find the recording
    recId = NP.SE.GetID(caseTb.clusId(i));
    isRec = recIdList==recId;
    senTb = senTbs{isRec};
    
    % Find the sentence
    isSen = startsWith(senTb.stimText, caseTb.startText(i));
    se = senTb.se(isSen).Duplicate;
    
    % Time alignment and sorting
    se.AlignTime('stimOff', 'taskTime');
    [~, I] = sort(se.GetTable('taskTime').cue3On, 'ascend');
    se.SortEpochs(I);
    
    caseTb.se(i) = se;
end

%% Collect samples

for i = 1 : height(caseTb)
    se = caseTb.se(i);
    uIdx = find(NP.Unit.GetClusTb(se).clusId == caseTb.clusId(i), 1);
    
    tt = se.GetTable('taskTime');
    tt.prodOn = cellfun(@(x) x(1), tt.prodOn);
    delayDur = tt.cue3On - tt.stimOff;
    RT = tt.prodOn - tt.cue3On;
    
    tWin = [tt.stimOff, tt.cue3On];
    st = se.SliceEventTimes('spikeTime', tWin, [], uIdx);
    rDelay = cellfun(@numel, st{:,:}) ./ delayDur;
    
    tWin = [tt.cue3On, tt.prodOn];
    st = se.SliceEventTimes('spikeTime', tWin, [], uIdx);
    rInit = cellfun(@numel, st{:,:}) ./ RT;
    
    caseTb.delayDur{i} = delayDur;
    caseTb.RT{i} = RT;
    caseTb.rDelay{i} = rDelay;
    caseTb.rInit{i} = rInit;
end

%% Plot

f = MPlot.Figure(4511); clf

nRows = height(caseTb);
colDist = [2 1 1 1];
nCols = numel(colDist);
tl = tiledlayout(nRows, sum(colDist));
tl.Padding = "compact";

tWin = [-.5 4];

for i = 1 : height(caseTb)
    % Raster
    se = caseTb.se(i);
    uIdx = find(NP.Unit.GetClusTb(se).clusId == caseTb.clusId(i), 1);
    
    ntArgs = MPlot.FindTileInd(nRows, colDist, i, 1);
    ax = nexttile(ntArgs{:});
    NP.UnitPlot.Raster(ax, se, [], tWin, uIdx);
    NP.TaskBaseClass.PlotEventBoundary(ax, se, [], tWin, ["stimOff", "cue3On", "prodOn"]);
    ax.XLabel.String = "Time from stim offset (s)";
    ax.XLim = [-.5 3];
    MPlot.Axes(ax);
    
    % Between delay spike rates and delay duration
    x = caseTb.delayDur{i};
    y = caseTb.rDelay{i};
    lm = fitlm(x, y);
    pval = lm.Coefficients.pValue(2);
    if pval > 0.05
        pval = NaN;
    end
    
    ntArgs = MPlot.FindTileInd(nRows, colDist, i, 2);
    ax = nexttile(ntArgs{:});
    plot(lm);
    legend off
    ax.YLim(1) = 0;
    ax.Title.String = sprintf("p = %.2f", pval);
    ax.XLabel.String = "Delay duration (s)";
    ax.YLabel.String = "Delay resp. (spk/s)";
    MPlot.Axes(ax);
    
    % Between delay spike rates and delay duration
    x = caseTb.RT{i};
    y = caseTb.rDelay{i};
    lm = fitlm(x, y);
    pval = lm.Coefficients.pValue(2);
    if pval > 0.05
        pval = NaN;
    end
    
    ntArgs = MPlot.FindTileInd(nRows, colDist, i, 3);
    ax = nexttile(ntArgs{:});
    plot(lm);
    legend off
    ax.YLim(1) = 0;
    ax.Title.String = sprintf("p = %.2f", pval);
    ax.XLabel.String = "RT (s)";
    ax.YLabel.String = [];
    MPlot.Axes(ax);
    
    % Between delay spike rates and init spike rate
    x = caseTb.rInit{i};
    y = caseTb.rDelay{i};
    lm = fitlm(x, y);
    pval = lm.Coefficients.pValue(2);
    if pval > 0.05
        pval = NaN;
    end
    
    ntArgs = MPlot.FindTileInd(nRows, colDist, i, 4);
    ax = nexttile(ntArgs{:});
    plot(lm);
    legend off
    % ax.XLim(1) = 0;
    ax.YLim(1) = 0;
    ax.Title.String = sprintf("p = %.2f", pval);
    ax.XLabel.String = "Init resp. (spk/s)";
    ax.YLabel.String = [];
    MPlot.Axes(ax);
    
    % return
end

MPlot.Paperize(f, 1.5, 2, 'FontSize', 5);
exportgraphics(f, fullfile(figDir, "bridge_delay_corr.png"));

return
%% Plot rasters

f = MPlot.Figure(4411); clf
tl = tiledlayout(3,3);
tl.Padding = "compact";

tWin = [-.5 4];

for i = 1 : height(caseTb)
    se = caseTb.se(i);
    uIdx = find(NP.Unit.GetClusTb(se).clusId == caseTb.clusId(i), 1);
    
    ax = nexttile;
    NP.UnitPlot.Raster(ax, se, [], tWin, uIdx);
    NP.TaskBaseClass.PlotEventBoundary(ax, se, [], tWin, ["stimOff", "cue3On", "prodOn"]);
    ax.XLabel.String = "Time from stim offset (s)";
    ax.XLim = [-.5 3];
    MPlot.Axes(ax);
end

%% Plot correlations

f = MPlot.Figure(4511); clf
tl = tiledlayout(3,9);
tl.Padding = "compact";

for i = 1 : height(caseTb)
    % Between delay spike rates and delay duration
    x = caseTb.delayDur{i};
    y = caseTb.rDelay{i};
    lm = fitlm(x, y);
    pval = lm.Coefficients.pValue(2);
    if pval > 0.05
        pval = NaN;
    end
    
    ax = nexttile;
    plot(lm);
    legend off
    ax.YLim(1) = 0;
    ax.Title.String = sprintf("u%i, p=%.2f", caseTb.clusId(i), pval);
    ax.XLabel.String = "Delay duration (s)";
    ax.YLabel.String = "Delay resp. (spk/s)";
    MPlot.Axes(ax);
    
    % Between delay spike rates and delay duration
    x = caseTb.RT{i};
    y = caseTb.rDelay{i};
    lm = fitlm(x, y);
    pval = lm.Coefficients.pValue(2);
    if pval > 0.05
        pval = NaN;
    end
    
    ax = nexttile;
    plot(lm);
    legend off
    ax.YLim(1) = 0;
    ax.Title.String = sprintf("u%i, p=%.2f", caseTb.clusId(i), pval);
    ax.XLabel.String = "RT (s)";
    ax.YLabel.String = "Delay resp. (spk/s)";
    MPlot.Axes(ax);
    
    % Between delay spike rates and init spike rate
    x = caseTb.rInit{i};
    y = caseTb.rDelay{i};
    lm = fitlm(x, y);
    pval = lm.Coefficients.pValue(2);
    if pval > 0.05
        pval = NaN;
    end
    
    ax = nexttile;
    plot(lm);
    legend off
    % ax.XLim(1) = 0;
    ax.YLim(1) = 0;
    ax.Title.String = sprintf("u%i, p=%.2f", caseTb.clusId(i), pval);
    ax.XLabel.String = "Init resp. (s)";
    ax.YLabel.String = "Delay resp. (spk/s)";
    MPlot.Axes(ax);
    
    % return
end
