%% Characterize LMV task design

anaDir = LMV.Data.GetAnalysisDir('task');
srcTb = LMV.Data.FindSource([]);

%% Compute the distribution of sentence interval

stimTb = readtable(fullfile(NP.Data.GetProjectRoot, "code", "tasks", "lmv", "TIMIT_selected_221011.xlsx"));
stimId = string(stimTb.id);

% Get unique sentence IDs
[uniqueIds, ~, idx] = unique(stimId);
counts = accumarray(idx, 1);

% Calculate intervals for each unique sentence
senItvl = cell(length(uniqueIds), 1);
for i = 1:length(uniqueIds)
    % Find all occurrences of this sentence
    trialInd = find(stimId == uniqueIds(i));
    
    % Calculate intervals between consecutive occurrences, accounting for cycles
    itvl = diff(trialInd);
    itvl = [itvl; numel(stimId)-sum(itvl)];
    senItvl{i} = itvl - 1;
end

% Compute group stats
allItvl = cell2mat(senItvl);
mItvl = mean(allItvl);
sdItvl = std(allItvl);

sen4Itvl = cell2mat(senItvl(counts > 1));
mSen4Itvl = mean(sen4Itvl);
sdSen4Itvl = std(sen4Itvl);

% Create a box plot of intervals
f = MPlot.Figure; clf
for i = 1:length(uniqueIds)
    plot(i, senItvl{i}, 'ko');
    hold on
end
ax = gca;
xlim([0 length(uniqueIds)+1]);
ylim([0 numel(stimId)]);
xticks(1:length(uniqueIds));
xticklabels(uniqueIds);
xtickangle(45);
ax.TickLabelInterpreter = "none";
xlabel('Sentence ID');
ylabel('# of trials');
title('Intervals Between Repeated Sentences');
MPlot.Axes(ax);
MPlot.Paperize(f, 1, .7);
exportgraphics(f, fullfile(anaDir, 'sentence_repeat_intervals.png'));

% Write stats to file
fileId = fopen(fullfile(anaDir, 'task_design_stats.txt'), 'w');
fprintf(fileId, 'Number of unique sentences: %d\n', length(uniqueIds));
fprintf(fileId, 'Interval of adjacent repeats (mean ± SD):\n');
fprintf(fileId, '  All sentences: %.2f±%.2f trials\n', mItvl, sdItvl);
fprintf(fileId, '  8-repeat sentences: %.2f±%.2f trials\n', mSen4Itvl, sdSen4Itvl);
fclose(fileId);
