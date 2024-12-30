%% 

anaDir = LMV.Data.GetAnalysisDir('phase_resp', 'signrank', 'preference');
srcTb = LMV.Data.FindSource([]);

%% Load test data

resp = LMV.Resp.LoadPhaseResponse(srcTb.recId);
rTest = LMV.Resp.LoadPhaseResponseTest(srcTb.recId);

%% Compute phase tunings

phaseNames = ["atten", "stim", "delay", "init", "prod"];

nRec = height(srcTb);
nPhase = numel(phaseNames);
clusTbs = cell(nRec, 1);

for i = 1 : nRec
    fprintf("\nCompute phase tunings for %s\n", srcTb.recId{i});
    
    % Unpack variables
    respTb = resp(i).respTb;
    dropoutTb = resp(i).dropoutTb;
    phaseTb = resp(i).phaseTb;
    stimIdTb = resp(i).stimIdTb;
    clusTb = resp(i).clusTb;
    
    % Set dropout periods to NaN
    respTb{:,:}(dropoutTb{:,:}) = NaN;
    
    % Subtract baseline activity
    rBase = respTb{resp(i).phaseTb.baseline,:};
    respTb{:,:} = respTb{:,:} - repmat(rBase, [width(phaseTb) 1]);
    
    % Keep responses for the phases of interest
    phaseTb = phaseTb(:,phaseNames);
    isOut = ~any(phaseTb{:,:}, 2);
    respTb(isOut,:) = [];
    phaseTb(isOut,:) = [];
    
    % Compute tuning curves
    G = sum(phaseTb{:,:}.*(1:nPhase), 2);
    T = NaN(width(respTb), nPhase, 4);
    for j = 1 : width(respTb)
        T(j,:,1) = splitapply(@(x) MMath.MeanStats(x), respTb.(j), G);
        T(j,:,4) = splitapply(@(x) sum(~isnan(x)), respTb.(j), G)';
    end
    clusTb.tuningPhase = repmat(phaseNames, [height(clusTb) 1]);
    clusTb.tuningMean = T(:,:,1);
    clusTb.tuningN = T(:,:,4);
    
    % Normalize tuning curves
    clusTb.tuningScore = MMath.Normalize(clusTb.tuningMean', 'minmax')';
%     clusTb.tuningScore = exp(clusTb.tuningMean) ./ sum(exp(clusTb.tuningMean), 2, 'omitnan');
    
    % Add masks
    clusTb.isResp = rTest(i).sigTb{:,phaseNames}; % response not significant
    
    % Exclude values
    m = false(size(clusTb.tuningScore));
%     m = m | clusTb.tuningN < 20; % computed with too few samples
    m = m | rTest(i).sigTb{:,phaseNames} == 0; % response not significant
    m = m | abs(clusTb.tuningMean) < 1; % small effect size
    clusTb.tuningValid = ~m;
    
    clusTbs{i} = clusTb;
end

clusTbCat = cat(1, clusTbs{:});

%% Plot tuning violin-scatter plots

f = MPlot.Figure(6575); clf
LMV.Resp.PlotPhaseRespViolin(clusTbCat, phaseNames);
MPlot.Paperize(f, 'ColumnsWide', 2.5, 'ColumnsHigh', 0.3);
exportgraphics(f, fullfile(anaDir, "phase_tuning_scores.png"));

%% Plot distributions of the differential selectivity

regions = "mPrCG";
regions = ["mPrCG", "STG", "vPrCG", "IFG"];
yMax = [0.15 0.3 0.2 0.25];

% pairs = ["stim", "prod"; "stim", "delay"; "delay", "init"; "init", "prod"];
pairs = nchoosek(phaseNames(2:end), 2);
pairs = pairs([3 1 6 4 2 5], :);
pairs = flip(pairs, 2);

f = MPlot.Figure(2997); clf
tl = tiledlayout(4, size(pairs,1));
tl.Padding = 'compact';

for i = 1 : numel(regions)
    for j = 1 : size(pairs,1)
        ax = nexttile;
        LMV.Resp.PlotDiffPhaseResp(clusTbCat, regions(i), pairs(j,:));
        ax.YLim(2) = yMax(i);
        if i ~= 1
            ax.Title.String = [];
        end
        if j ~= 1
            ax.YLabel.String = [];
            ax.YTickLabel = [];
        end
        if i ~= numel(regions)
            ax.XLabel.String = [];
        end
    end
end

MPlot.Paperize(f, 'ColumnsWide', 2.5, 'ColumnsHigh', 0.8);
exportgraphics(f, fullfile(anaDir, sprintf("diff_response_histograms.png")));

return
%% Plot interactive scatter plots of response selectivity

regions = ["mPrCG", "vPrCG", "IFG", "STG"];

pairs = ["stim", "prod"; "stim", "delay"; "delay", "init"; "init", "prod"];
pairs = flip(pairs, 2);

f = MPlot.Figure(2987); clf
tl = tiledlayout("flow");
tl.Padding = 'compact';

for j = 1 : size(pairs,1)
    ax = nexttile;
    LMV.Resp.PlotPairedPhaseResp(clusTbCat, regions, pairs(j,:));
end

MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', 1);
exportgraphics(f, fullfile(anaDir, sprintf("paired_response_scatter.png")));

