%% Quantify LMV performance across subjects

anaDir = LMV.Data.GetAnalysisDir('task_perf');
srcTb = LMV.Data.FindSource([]);

%%

perfTb = table;
perfTb.recId = srcTb.recId;

stimIdList = ["all" LMV.Param.stimIdList4' setdiff(LMV.Param.stimIdList14, LMV.Param.stimIdList4)'];

for i = 1 : height(srcTb)
    % Load se
    se = NP.SE.LoadSession(srcTb.path{i});
    se = LMV.SE.Transform(se);
    
    % % Split tasks and keep LMV only
    % seTaskTb = NP.SE.SplitConditions(se, 'conditionVars', {'taskName'});
    % disp(seTaskTb);
    % se = seTaskTb.se(seTaskTb.taskName=="lmv");
    % if isempty(se)
    %     continue
    % end
    
    % Compute average trial time
    rt = se.GetReferenceTime('taskTime');
    dur = rmoutliers(diff(rt)); % exclude pauses or inter-task time
    perfTb.avgTrialDur(i) = mean(dur);
    perfTb.taskDur(i) = sum(dur) + perfTb.avgTrialDur(i);
    
    % Match production to stim
    [tt, tv] = se.GetTable('taskTime', 'taskValue');
    [tt, tv] = LMV.SE.MatchProd2Stim(tt, tv, rt);
    [tt, tv] = LMV.SE.FindMorphTimes(tt, tv, rt);
    
    % Fraction of successful repetition (at least half of the sentence is correct, defined by normalized edit distance > 0.5)
    ra = tv.alignScore;
    rc = -Inf(size(ra));
    for j = 1 : numel(rc)
        if isempty(tv.alignInfo{j})
            continue
        else
            rc(j) = tv.alignInfo{j}.nscore;
        end
    end
    g = rc > 0.5;
    
    % Median reaction time (from cue3 onset to prod onset)
    RT = tv.tReact;
    
    % Relative speech rate (ratio of matched prod duration to stim duration)
    r = (tt.prodMatchOff - tt.prodMatchOn) ./ (tt.stimMatchOff - tt.stimMatchOn);
    
    % Other durations
    cue2stim = tt.stimOn - tt.cue1On;
    stimDur = tt.stim.GetDuration;
    delayDur = tt.cue3On - tt.stimOff;
    prodDur = cellfun(@(a,b) b(end)-a(1), tt.prodOn, tt.prodOff);
    
    % Compute metrics for each sentence
    for j = numel(stimIdList) : -1 : 1
        if stimIdList(j) == "all"
            m = true(size(tv.stimId));
        else
            m = tv.stimId == stimIdList(j);
        end
        perfTb.repAccuracy(i,j) = median(ra(m));
        perfTb.repConsistency(i,j) = median(rc(m));
        perfTb.includeRate(i,j) = mean(g(m));
        perfTb.relativeProdRate(i,j) = median(r(g & m));
        perfTb.cue2stim(i,j) = median(cue2stim(m));
        perfTb.stimDur(i,j) = mean(stimDur(m));
        perfTb.delayDur(i,j) = mean(delayDur(m));
        perfTb.RT(i,j) = median(RT(g & m), 'omitnan');
        perfTb.prodDur(i,j) = median(prodDur(g & m));
    end
end

cachePath = fullfile(anaDir, 'computed_lmv_perf.mat');
save(cachePath, 'srcTb', 'stimIdList', 'perfTb', 'tt', 'tv');

%% Report task time

fileID = fopen(fullfile(anaDir, 'task_time_stats.txt'), 'w');

fprintf(fileID, '%i sessions\n', height(perfTb));
fprintf(fileID, 'Median trial duration %.2g ± %.2gs (mean±SD)\n', mean(perfTb.avgTrialDur), std(perfTb.avgTrialDur));
fprintf(fileID, 'Median 72 trial duration %.2g ± %.2gmin (mean±SD)\n', mean(perfTb.avgTrialDur)*72/60, std(perfTb.avgTrialDur)*72/60);
fprintf(fileID, 'Task duration %.2g ± %.2gmin (mean±SD)\n', mean(perfTb.taskDur)/60, std(perfTb.taskDur)/60);
fprintf(fileID, 'Shortest task duration %.2gmin\n', min(perfTb.taskDur)/60);
fprintf(fileID, 'Longest task duration %.2gmin\n\n', max(perfTb.taskDur)/60);

fprintf(fileID, formattedDisplayText(perfTb(:,1:3), 'SuppressMarkup', true)+"\n\n");
fclose(fileID);

%% Plot performance

stimObj = arrayfun(@(x) tt.stim(find(tv.stimId==x,1)).Cut, stimIdList(2:end), 'Uni', false);
stimText = cellfun(@(x) x(1:2).GetAllParentLabel, stimObj);
stimText = ["all" stimText];

f = MPlot.Figure(8223); clf
f.WindowState = 'maximized';
tl = tiledlayout('flow');
var2plot = ["repAccuracy", "repConsistency", "includeRate", "relativeProdRate", "cue2stim", "stimDur", "delayDur", "RT", "prodDur"];
for i = 1 : numel(var2plot)
    ax = nexttile();
    x = 1 : numel(stimIdList);
    Y = perfTb.(var2plot(i));
    m = mean(Y, 'omitnan');
    plot(x', Y, 'Color', [0 0 0 .15]); hold on
%     plot(x, m, 'Color', 'k', 'LineWidth', 2);
    boxplot(Y);
    xticks(x);
    xticklabels(stimText);
    title(var2plot(i));
    switch var2plot(i)
        case "includeRate", ylim([0 1.1]); ylabel("Fraction");
        case "relativeProdRate", ylabel("Fold wrt stim");
        case {"cue2stim", "stimDur", "delayDur", "RT", "prodDur"}, ylabel("Second");
        case {"repAccuracy", "repConsistency"}, ylabel("Norm edit dist");
    end     
end

exportgraphics(f, fullfile(anaDir, 'lmv_perf.png'));
