%% Quantify kernel peaks

mdlName = LMV.Linker.currentModel;
anaDir = LMV.Data.GetAnalysisDir("linker", mdlName);
clusTb = LMV.Linker.LoadClusTb();

%% Compute and save peak stats

statsFile = fullfile(anaDir, "kernel_peak_stats.txt");
fileID = fopen(statsFile, "w");

% Peak stats within each group
pk = struct;
p0 = struct;
for g = ["mirror", "feedback"]
    % Find peaks
    m = clusTb.hcGroup==g;
    [~, pks] = cellfun(@(x) findpeaks(x.Beta, x.dt, NPeaks=1, SortStr="descend"), clusTb.linker(m));
    if g == "mirror"
        pks(15) = []; % true peak is out of bound
    elseif g == "feedback"
        pks(end) = 0; % detected the wrong peak
    end
    pk.(g) = pks;

    % Signrank test
    [~, pNorm] = lillietest(pks);        % Lilliefors test for normality
    [p0.(g), ~, stats] = signrank(pks);  % Signrank test for median difference

    % Effect size
    n = numel(pks);             % number of non-zero observations
    T = n*(n+1)/2;              % total of ranks 1…n  (needed below)
    Wplus = stats.signedrank;   % W+
    r_rb = (2*Wplus - T) / T;   % == (W+ - W-)/T
    
    % Print stats
    fprintf(fileID, "%s kernel peak time stats:\n", g);
    fprintf(fileID, "Mean ± SD: %.2f ± %.2f. Median (IQR): %.2f (%.2f-%.2f)\n", mean(pks), std(pks), median(pks), prctile(pks,25), prctile(pks,75));
    fprintf(fileID, "Lilliefors test p-val = %.2e\n", pNorm);
    fprintf(fileID, "Signrank test different from 0s with p-val = %.2e\n", p0.(g));
    fprintf(fileID, "Effect size r_rb = %.2f\n\n", r_rb);
end

% Peak difference between mirror and feedback
pDiff = ranksum(pk.mirror, pk.feedback);
fprintf(fileID, "Ranksum test mirror and feedback kernel peak times are different at p-val = %e\n", pDiff);

fclose(fileID);

%% Plot peak distributions

% Create box plot
f = MPlot.Figure(95); clf

boxdata = [pk.mirror; pk.feedback];
grouping = [ones(size(pk.mirror)); 2*ones(size(pk.feedback))];
boxplot(boxdata, grouping, 'Labels', {'Mirror', 'Feedback'}, 'Notch', 'on'); hold on;

plot([.5 2.5], [0 0], '-', 'Color', [0 0 0 .2]); hold on;

ylabel('Peak Time (s)');
title('Kernel Peak Times');

% Add individual points
violinArgs = {100, 'Span', 0.3, 'Color', 'k'};
MPlot.ViolinScatter(1, pk.mirror, violinArgs{:});
MPlot.ViolinScatter(2, pk.feedback, violinArgs{:});

% Format p-values for display
function str = format_pval(p)
    pTh = [0.05, 0.01, 0.001];
    iTh = find(p < pTh, 1, 'last');
    if isempty(iTh)
        str = sprintf('n.s.');
    else
        str = repelem('*', iTh); % sprintf('p=%.3f', p);
    end
end

% Add text annotations
text(1, 0.1, sprintf('%s', format_pval(p0.mirror)), 'HorizontalAlignment', 'center');
text(2, 0.1, sprintf('%s', format_pval(p0.feedback)), 'HorizontalAlignment', 'center');
text(1.5, mean([median(pk.feedback), median(pk.mirror)]), ...
    sprintf('%s', format_pval(pDiff)), 'HorizontalAlignment', 'center');

ax = gca;
ax.YLim = [-.4 .15];
ax.YTick = ax.YLim(1):0.1:ax.YLim(2);
ax.XTickLabel = ["tPP", "Aud."];

MPlot.Axes(gca);
MPlot.Paperize(f, .32, .27);
exportgraphics(f, fullfile(anaDir, "kernel_peak_dists.png"));
exportgraphics(f, fullfile(anaDir, "kernel_peak_dists.pdf"));