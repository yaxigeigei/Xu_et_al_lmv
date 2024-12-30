%% 

anaDir = LMV.Data.GetAnalysisDir('units');
srcTb = LMV.Data.FindSource([]);

%% Compute unit quality stats

clusTbs = cell(height(srcTb), 1);

for i = 1 : height(srcTb)
    % Check cache
    cacheDir = fullfile(anaDir, 'computed_qc');
    cachePath = fullfile(cacheDir, srcTb.recId{i}+"_clusTb.mat");
    
    if exist(cachePath, 'file')
        % Load cached
        fprintf("\nLoad computed QC clusTb from cache:\n%s\n", cachePath);
        load(cachePath, 'clusTb');
    else
        % Compute from se
        sePath = srcTb.path{i};
        fprintf("\nCompute QC clusTb from se:\n%s\n", sePath);
        se = NP.SE.LoadSession(sePath, 'UserFunc', @(x) x.RemoveTable('LFP'));
        clusTb = NP.Unit.ComputeQualityTable(se, 'lmv');
        
        % Cache clusTb
        if ~exist(cacheDir, 'dir')
            mkdir(cacheDir);
        end
        save(cachePath, 'clusTb');
        fprintf("Saved QC clusTb to cache:\n%s\n", cachePath);
        
        % Create figure folder
        recId = NP.SE.GetID(se);
        recDir = fullfile(srcTb.folder{i}, recId);
        if ~exist(recDir, 'dir')
            mkdir(recDir);
        end
        
        % Convert to short cluster IDs
        clusTb.clusId = mod(clusTb.clusId, 1e5);
        se.userData.ksMeta.clusTb = clusTb;
        
        % Plot span
        spanFigPath = fullfile(recDir, recId+"_cluster_span");
        NP.UnitPlot.Span(se, [], spanFigPath);
        
        % Plot waveform
        wfFigBasePath = fullfile(recDir, recId+"_cluster_waveform");
        NP.UnitPlot.WaveformArray(clusTb, wfFigBasePath);
        
        % Plot ISI
        isiFigBasePath = fullfile(recDir, recId+"_cluster_isi");
        NP.UnitPlot.ISIArray(clusTb, isiFigBasePath);
    end
    
    clusTbs{i} = clusTb;
end

%% Session and area summaries

% Summarize unit quality by sessions
sqTb = cell(size(clusTbs));
for i = 1 : numel(clusTbs)
    clusTb = clusTbs{i};
    tbRow = clusTb(1, {'recId', 'subjectId', 'region'});
    tbRow.numUnits = height(clusTb);
    tbRow.numRPV = sum(clusTb.RPV_ND > NP.Param.maxRPV);
    tbRow.numContam = sum(clusTb.contamND > NP.Param.maxContam);
    tbRow.numMulti = sum(clusTb.contamND > NP.Param.maxContam | clusTb.RPV_ND > NP.Param.maxRPV);
    tbRow.numSingle = tbRow.numUnits - tbRow.numMulti;
    tbRow.numShort = sum(clusTb.fracSpan < NP.Param.minFracSpan);
    sqTb{i} = tbRow;
end
sqTb = vertcat(sqTb{:});

% Set order to the areas
sqTb.region = categorical(sqTb.region, {'mPrCG', 'vPrCG', 'IFG', 'STG', 'SMG', 'SFG'}, 'Ordinal', true);

% Summarize each area
[groupId, regTb] = findgroups(sqTb(:,'region'));
for i = 1 : width(sqTb)
    vn = sqTb.Properties.VariableNames{i};
    val = sqTb.(vn);
    if isnumeric(val)
        regTb.(vn) = splitapply(@sum, val, groupId);
    elseif ~strcmp(vn, 'region')
        regTb.(vn) = splitapply(@(x) {unique(x)}, val, groupId);
    end
end

% Add total #unit in each region back to session table
N = splitapply(@sum, sqTb.numUnits, groupId);
sqTb.numRegionUnits = N(groupId);

%% Summary stats in text

% Unpack variables
qcTbCat = vertcat(clusTbs{:});
rMean = qcTbCat.meanActiveRate;
RPV = qcTbCat.RPV_ND;
C = qcTbCat.contamND;
span = qcTbCat.fracSpan;
isRPV = RPV > NP.Param.maxRPV;
isContam = C > NP.Param.maxContam;
isSingle = ~(isRPV | isContam);
isShort = span < NP.Param.minFracSpan;

% Print overall percentages
fileID = fopen(fullfile(anaDir, 'unit_quality_stats.txt'), 'w');
fprintf(fileID, '%d sessions\n', height(sqTb));
fprintf(fileID, '%d units (%g ± %g per session)\n', ...
    sum(sqTb.numUnits), mean(sqTb.numUnits), std(sqTb.numUnits));
fprintf(fileID, '%d units cover less than %g%% of task time\n', ...
    sum(sqTb.numShort), NP.Param.minFracSpan*100);
fprintf(fileID, '%d single-units (%g ± %g per session)\n', ...
    sum(sqTb.numSingle), mean(sqTb.numSingle), std(sqTb.numSingle));
fprintf(fileID, '%d multi-units (%g ± %g per session)\n', ...
    sum(sqTb.numMulti), mean(sqTb.numMulti), std(sqTb.numMulti));
fprintf(fileID, '%.1f%% units failed refractory period violation threshold at %g%%\n', ...
    mean(isRPV)*100, NP.Param.maxRPV);
fprintf(fileID, '%.1f%% failed contamination threshold at %g%%\n', ...
    mean(isContam)*100, NP.Param.maxContam);
fprintf(fileID, '%.1f%% failed both\n\n', ...
    mean(isRPV | isContam)*100);

fprintf(fileID, formattedDisplayText(regTb, 'SuppressMarkup', true)+"\n\n");
fprintf(fileID, formattedDisplayText(sqTb, 'SuppressMarkup', true)+"\n\n");

fclose(fileID);

%% Single vs multi-unit metrics

f = MPlot.Figure(3111); clf

subplot(2,2,1);
RPVlim = 3;
PRVcut = MMath.Bound(RPV, [0 RPVlim]);
histogram(PRVcut, 0:.02:RPVlim, 'Normalization', 'cumcount', ...
    'EdgeColor', 'none', 'FaceColor', [0 0 0]); hold on
histogram(PRVcut(isSingle), 0:.02:RPVlim, 'Normalization', 'cumcount', ...
    'EdgeColor', 'none', 'FaceColor', [0 .7 0]);
ax = MPlot.Axes(gca);
ax.YLim = [0 numel(RPV)];
ax.XTick = 0:.5:RPVlim;
ax.YTick = linspace(0, numel(RPV), numel(0:.2:1));
ax.YTickLabel = 0:.2:1;
xlabel('Refractory period violation rate (%)');
ylabel('Fraction of units');

% subplot(2,2,2);
% histogram(rMean, 0:60, 'EdgeColor', 'none', 'FaceColor', [0 0 0]); hold on
% histogram(rMean(isSingle), 0:60, 'EdgeColor', 'none', 'FaceColor', [0 .7 0]);
% ax = MPlot.Axes(gca);
% xlabel('Mean spike rate (spk/s)');
% ylabel('# of units');

subplot(2,2,3);
MPlot.Blocks([0, RPVlim*1.01], [NP.Param.maxContam, 50*1.01], [0 0 0], 'FaceAlpha', .1); hold on
MPlot.Blocks([NP.Param.maxRPV, RPVlim*1.01], [0, 50*1.01], [0 0 0], 'FaceAlpha', .1);
plot(PRVcut(~isSingle), C(~isSingle), 'k.', 'MarkerSize', 3);
plot(PRVcut(isSingle), C(isSingle), '.', 'Color', [0 .7 0], 'MarkerSize', 5);
ax = MPlot.Axes(gca);
ax.XTick = 0:.5:RPVlim;
ax.YTick = 0:10:50;
axis tight
xlabel('Refractory period violation rate (%)');
ylabel('Contamination rate (%)');

subplot(2,2,4);
histogram(C, 0:.2:50, 'Normalization', 'cumcount', ...
    'EdgeColor', 'none', 'FaceColor', [0 0 0]); hold on
histogram(C(isSingle), 0:.2:50, 'Normalization', 'cumcount', ...
    'EdgeColor', 'none', 'FaceColor', [0 .7 0]);
ax = MPlot.Axes(gca);
ax.YLim = [0 numel(C)];
ax.XTick = 0:10:50;
ax.XTickLabel = ax.XTick;
ax.YTick = linspace(0, numel(C), numel(0:.2:1));
ax.YTickLabel = 0:.2:1;
xlabel('Contamination rate (%)');
ylabel('Fraction of units');


% % Number of units for each sessions
% regTbSorted = sortrows(regTb, 'numUnits', 'descend');
% sqTbSorted = sortrows(sqTb, {'numRegionUnits', 'numUnits'}, 'descend');
% 
% subplot(2,2,2);
% bar(sqTbSorted.numUnits, 'EdgeColor', 'none', 'FaceColor', 'k'); hold on
% bar(sqTbSorted.numSingle, 'EdgeColor', 'none', 'FaceColor', [0 .7 0]);
% [xGroups, yGroups] = MPlot.GroupRibbon(sqTbSorted.region, [400 500], ...
%     NP.Param.GetRegionColors(regTbSorted.region), ...
%     'Groups', regTbSorted.region);
% text(xGroups, 500*ones(size(xGroups)), string(regTbSorted.region), ...
%     'Horizontal', 'center', 'Vertical', 'bottom', 'Rotation', 0, 'FontSize', 12);
% ax = MPlot.Axes(gca);
% ax.TickLabelInterpreter = 'none';
% ax.XTick = 1 : height(sqTb);
% ax.XTickLabel = sqTbSorted.recId;
% ax.YScale = 'log';
% ax.YLim = [10 500];
% ax.YTick = [10 30 100 300];
% ylabel('# of units');

MPlot.Paperize(f, 'ColumnsWide', 1.2, 'AspectRatio', .7);
exportgraphics(f, fullfile(anaDir, "unit_quality.png"));
exportgraphics(f, fullfile(anaDir, "unit_quality.pdf"));

%% Unit yield plots

f = MPlot.Figure(3112); clf

% Show yield chronologically
sqTbSorted = sortrows(sqTb, 'recId', 'ascend');

ax = subplot(2,1,1);
bar(sqTbSorted.numUnits, 'EdgeColor', 'none', 'FaceColor', 'k'); hold on
bar(sqTbSorted.numSingle, 'EdgeColor', 'none', 'FaceColor', [0 .7 0]);

ax.XTick = 1 : height(sqTbSorted);
ax.XTickLabel = sqTbSorted.recId;
ax.TickLabelInterpreter = 'none';
ax.YLim = [10 500];
ax.YScale = 'log';
ax.YTick = [10 30 100 300];
ax.YLabel.String = '# of units';
MPlot.Axes(gca);

% Group recording by regions
sqTb.region = categorical(sqTb.region, LMV.Param.regions);
sqTbSorted = sortrows(sqTb, {'region', 'recId'}, {'ascend', 'ascend'});

ax = subplot(2,1,2);
bar(sqTbSorted.numUnits, 'EdgeColor', 'none', 'FaceColor', 'k'); hold on
bar(sqTbSorted.numSingle, 'EdgeColor', 'none', 'FaceColor', [0 .7 0]);
[xGroups, yGroups] = MPlot.GroupRibbon(sqTbSorted.region, [400 500], ...
    NP.Param.GetRegionColors(LMV.Param.regions), ...
    'Groups', LMV.Param.regions);
text(xGroups, 500*ones(size(xGroups)), string(LMV.Param.regions), ...
    'Horizontal', 'center', 'Vertical', 'bottom', 'Rotation', 0, 'FontSize', 12);

ax.XTick = 1 : height(sqTb);
ax.XTickLabel = sqTbSorted.recId;
ax.TickLabelInterpreter = 'none';
ax.YLim = [10 500];
ax.YScale = 'log';
ax.YTick = [10 30 100 300];
ax.YLabel.String = '# of units';
MPlot.Axes(gca);

MPlot.Paperize(f, 1, .75);
exportgraphics(f, fullfile(anaDir, "unit_yield.png"));
exportgraphics(f, fullfile(anaDir, "unit_yield.pdf"));

%% Plot example ISI histograms

% Randomly sample a subset of unit as examples
rng(61);
indEg = randsample(find(isSingle), 30);

f = MPlot.Figure(3116); clf
tl = tiledlayout(5, 6);
tl.Padding = "compact";
for i = 1 : numel(indEg)
    ax = nexttile;
    k = indEg(i);
    x = qcTbCat.isiEdges{k};
    x = x(1:end-1) + diff(x)/2;
    y = qcTbCat.isiCountND{k};
    
    h = bar(x, y, 'histc');
    h.FaceColor = [0 .7 0];
    h.EdgeColor = 'none';
    
    MPlot.Blocks([0 NP.Param.minISI], [0 max(y)], [1 0 0], 'FaceAlpha', .3);
    
    ax.XTick = [];
    ax.YTick = [];
    axis tight off
    xlim([0 0.02]);
    % ax.YLim(2) = ax.YLim(2)*1.5;
    text(ax.XLim(2), ax.YLim(2), sprintf("%.1f%%", qcTbCat.contamND(k)), "HorizontalAlignment", "right", "VerticalAlignment", "top");
end

MPlot.Paperize(f, 1, .6);
exportgraphics(f, fullfile(anaDir, "example_isi_histograms.png"));
exportgraphics(f, fullfile(anaDir, "example_isi_histograms.pdf"));

%% Plot example waveforms

f = MPlot.Figure(3118); clf
tl = tiledlayout(5, 6);
tl.Padding = "compact";

for i = 1 : numel(indEg)
    ax = nexttile;
    k = indEg(i);
    NP.UnitPlot.Waveform(ax, qcTbCat, k);
    
    ax.XTick = [];
    ax.YTick = [];
    ax.Title.String = [];
    axis tight off
    text(ax.XLim(1), ax.YLim(1), sprintf("%.1f", qcTbCat.SNR(k)), "HorizontalAlignment", "left", "VerticalAlignment", "top");
end

MPlot.Paperize(f, 1, .6);
exportgraphics(f, fullfile(anaDir, "example_waveforms.png"));
exportgraphics(f, fullfile(anaDir, "example_waveforms.pdf"), ContentType="vector");


