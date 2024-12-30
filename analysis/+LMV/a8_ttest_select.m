%% 

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'selectivity');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end

%% Load data

% Load responsiveness
tStruct = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'ttest', 'computed_ttest_clusTb.mat'));
clusTb = tStruct.clusTb;

% Keep remove non-responsive units
clusTb(clusTb.tId1=="none",:) = [];

%% Use overall response magnitudes

% Find the biggest response across significant tests for each task phase
phaseNames = tStruct.phaseNames;
pTb = clusTb(:, phaseNames);
sTb = NP.Resp.GetSigTable(clusTb, 'Side', 'right', 'Minimum', false);
rTb = clusTb(:, phaseNames+"Resp");
for i = 1 : width(rTb)
    R = rTb.(i);
    R(~sTb.(i)) = 0;
%     rTb.(i) = R(:,1);
    [~, k] = max(abs(R), [], 2);
    k = sub2ind(size(R), (1:size(R,1))', k);
    rTb.(i) = R(k);
end

% Construct task phase indices
R = rTb{:,:};
I = cumsum(ones(size(R)), 2);
I(R==0) = 6;
% I(R<0) = I(R<0) + numel(phaseNames); % indexing inhibition types after activation
I = categorical(I, 1:6, ["atten", "stim", "delay", "init", "prod", "none"]);

% Find best task phase activation
for i = 1 : numel(phaseNames)
    % Find the indices of maximum response
    [~, k] = max(abs(R), [], 2);
    k = sub2ind(size(R), (1:size(R,1))', k);
    tR = R(k);
    tId = I(k);
    
    % Inherit previous tier in the absence of further tier
    if i > 1
        m = tId=="none";
        tId(m) = clusTb.("tId"+(i-1))(m);
        tR(m) = clusTb.("tR"+(i-1))(m);
    end
    
    % Add results to table
    clusTb.("tR"+i) = tR;
    clusTb.("tId"+i) = tId;
    
    % Exclude max resp for next round
    R(k) = 0;
    I(k) = "none";
end

%% Plot interactive scatter plots of response selectivity

regions = ["mPrCG", "vPrCG", "IFG", "STG"];

f = MPlot.Figure(2987); clf
tl = tiledlayout(2,2);
tl.Padding = 'compact';

phasePair = flip(["stim", "prod"]);
NP.Resp.PlotSelectivity(clusTb, regions, phasePair);

phasePair = flip(["stim", "delay"]);
NP.Resp.PlotSelectivity(clusTb, regions, phasePair);

phasePair = flip(["delay", "init"]);
NP.Resp.PlotSelectivity(clusTb, regions, phasePair);

phasePair = flip(["init", "prod"]);
NP.Resp.PlotSelectivity(clusTb, regions, phasePair);

MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', 1);
exportgraphics(f, fullfile(anaDir, sprintf("paired_spike_rate_scatter.png")));

%% Plot distributions of the differential selectivity

f = MPlot.Figure(2997); clf
tl = tiledlayout(2,2);
tl.Padding = 'compact';

phasePair = flip(["stim", "prod"]);
NP.Resp.PlotDiffSelectivity(clusTb, regions, phasePair);

phasePair = flip(["stim", "delay"]);
NP.Resp.PlotDiffSelectivity(clusTb, regions, phasePair);

phasePair = flip(["delay", "init"]);
NP.Resp.PlotDiffSelectivity(clusTb, regions, phasePair);

phasePair = flip(["init", "prod"]);
NP.Resp.PlotDiffSelectivity(clusTb, regions, phasePair);

MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', 1);
exportgraphics(f, fullfile(anaDir, sprintf("paired_diff_spike_rate_histograms.png")));

