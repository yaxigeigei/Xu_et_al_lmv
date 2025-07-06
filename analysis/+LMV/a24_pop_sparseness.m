%% Export unit PETHs as Pandas dataframes

anaDir = LMV.Data.GetAnalysisDir("pop_dynamics", "sca", "sparsity");

%% Load data

% Load sentence average population spike rates (same input to SCA)
load(fullfile(LMV.Data.GetAnalysisDir, "data", "ce_m2_ex3_sentence-avg.mat"), 'ce');

% Load responsiveness test results
rTest = LMV.Resp.LoadPhaseResponseTest();
[~, I] = MMath.SortLike(rTest.clusTb.clusId, ce.clusTb.clusId);
allClusTb = rTest.clusTb(I,:);
allSigTb = rTest.sigTb(I,:);

%% Data preparation

% Keep responsive units
isResp = any(allSigTb{:,:}, 2);
clusTb = allClusTb(isResp, :);

respTb = ce.GetTable("resp");
R = cell2mat(respTb{:,2:end});
R = R(:, isResp);

% Soft normalization to peak with 5 Hz padding
R = MMath.Normalize(R, 'maxsoft', 5);

%% Compute sparseness

% regions = LMV.Param.regions;
regions = "mPrCG";

isRegion = clusTb.region == regions;
regR = R(:, isRegion);



% Compute population sparseness using kurtosis - 3 for each time step
% Kurtosis is a measure of the "tailedness" of the distribution
% Subtracting 3 gives us excess kurtosis (relative to normal distribution)
% Higher values indicate more sparse activity
sparseness = kurtosis(regR, 1, 2) - 3;

%% Plot the sparseness time series

f = MPlot.Figure(4644); clf
tl = tiledlayout("flow");
tl.Padding = "compact";

LMV.Sparsity.PlotKurtosis(ce, sparseness);





