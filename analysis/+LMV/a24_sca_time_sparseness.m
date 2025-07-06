%% Analyze sparseness of units that contribute strongly to SCA components

anaDir = LMV.Data.GetAnalysisDir("pop_dynamics", "sca", "sparsity");

%% Load data

% SCA results
load(fullfile(fileparts(anaDir), "regTb.mat"), "regTb");
region = "mPrCG"; % Analyze motor cortex by default
nComp = 12; % Number of components to analyze
regTb = regTb(regTb.cond==region & regTb.nComp==nComp, :);

% Load responsiveness test results
rTest = LMV.Resp.LoadPhaseResponseTest();
isRegion = rTest.clusTb.region == region;
isResp = any(rTest.sigTb{:,:}, 2);
clusTb = rTest.clusTb(isRegion & isResp, :);

%% Sort SCA components by variance explained

[~, veInd] = sort(regTb.relVE{1}, 'descend');
% Reorder components if needed for specific brain regions
switch region
    case 'mPrCG'
        veInd = veInd([1:3 5 4 6:end]);
    otherwise
end
Z = regTb.Z{1}(:,veInd);
V = regTb.V{1}(veInd,:);
CI = regTb.ZnullCI{1}(:,veInd,:);
relVE = regTb.relVE{1}(veInd);

%% Load sparseness data

% Load pre-concatenated sparseness table
fprintf('Loading sparseness data...\n');
combinedPath = fullfile(LMV.Data.GetAnalysisDir, "units", "computed_sparsity_clusTb.mat");
sSP = load(combinedPath, 'clusTb');
spTb = sSP.clusTb;

% Filter for units in the selected region
spTb = spTb(ismember(spTb.clusId, clusTb.clusId), :);
fprintf('Loaded %d units with sparseness data\n', height(spTb));

% Calculate minimum kurtosis for each unit across time windows
timeWindows = ["stim", "delay", "prod", "trial"];
minKurtosis = zeros(height(spTb), 1);
for u = 1:height(spTb)
    kurtValues = zeros(numel(timeWindows), 1);
    for w = 1:numel(timeWindows)
        kurtValues(w) = spTb.sparsityTb{u}{timeWindows(w), 'kurtosis'};
    end
    minKurtosis(u) = min(kurtValues);
end

%% Identify units with strong vs weak contribution to different subspaces

loadingThreshHigh = 80; % Percentile threshold for high-contributing units
loadingThreshLow = 20;  % Percentile threshold for low-contributing units

% Define subspaces to analyze
subspaces = struct();
subspaces(1).name = 'Full Space';
subspaces(1).comps = 1:12;
subspaces(2).name = 'Components 1-4';
subspaces(2).comps = 1:4;
subspaces(3).name = 'Components 5-12';
subspaces(3).comps = 5:12;

% Create figure for all comparisons
f = MPlot.Figure(450); clf
tl = tiledlayout(1, numel(subspaces), 'Padding', 'compact', 'TileSpacing', 'compact');

% Process each subspace
results = struct();
for s = 1:numel(subspaces)
    % Calculate variance contribution for this subspace, scaled by relative variance explained
    varContrib = sum(V(subspaces(s).comps,:).^2 .* relVE(subspaces(s).comps)', 1)';
    
    % Identify high and low contributing units
    highThresh = prctile(varContrib, loadingThreshHigh);
    lowThresh = prctile(varContrib, loadingThreshLow);
    isHighContrib = varContrib > highThresh;
    isLowContrib = varContrib < lowThresh;
    
    % Get corresponding sparseness data
    highClusIds = clusTb.clusId(isHighContrib);
    lowClusIds = clusTb.clusId(isLowContrib);
    highUnitIdx = ismember(spTb.clusId, highClusIds);
    lowUnitIdx = ismember(spTb.clusId, lowClusIds);
    highKurtosis = minKurtosis(highUnitIdx);
    lowKurtosis = minKurtosis(lowUnitIdx);
    
    % Statistical test
    [~, pval] = ranksum(log10(highKurtosis), log10(lowKurtosis));
    
    % Create plot
    ax = nexttile;
    [highF, highX] = ecdf(log10(highKurtosis));
    [lowF, lowX] = ecdf(log10(lowKurtosis));
    plot(highX, highF, 'Color', 'r', 'LineWidth', 2); hold on
    plot(lowX, lowF, 'Color', 'b', 'LineWidth', 2);
    
    % Add labels and formatting
    title(subspaces(s).name);
    ylabel('Cumulative Probability');
    xlabel('log_{10}(Minimum Kurtosis)');
    
    % Add p-value annotation
    if pval < 0.05
        text(0.1, 0.9, sprintf('p=%.1f*', pval), 'Units', 'normalized', 'HorizontalAlignment', 'left');
    else
        text(0.1, 0.9, sprintf('p=%.1f', pval), 'Units', 'normalized', 'HorizontalAlignment', 'left');
    end
    
    MPlot.Axes(ax);
    
    % Store results
    results(s).name = subspaces(s).name;
    results(s).pval = pval;
    results(s).highMean = mean(log10(highKurtosis), 'omitnan');
    results(s).lowMean = mean(log10(lowKurtosis), 'omitnan');
end

% Add legend to the first subplot
legend({'High contribution', 'Low contribution'}, 'Location', 'eastoutside');

% Adjust figure
MPlot.Paperize(f, 1.5, .4);
exportgraphics(f, fullfile(anaDir, sprintf('%s_subspace_sparseness.png', region)));

% % Save statistical results
% statResults = table();
% statResults.Subspace = {results.name}';
% statResults.HighMean = [results.highMean]';
% statResults.LowMean = [results.lowMean]';
% statResults.PValue = [results.pval]';
% save(fullfile(anaDir, 'subspace_sparseness_stats.mat'), 'statResults');
