%% 

anaDir = LMV.Data.GetAnalysisDir("pca");

%% Load data

% Load grand average M2 PETHs
load(fullfile(LMV.Data.GetAnalysisDir, "data", "ce_m2_grand-avg.mat"), 'ce');

% Add task phase events
phases = struct;
phases.atten = {'cue1On', 'stimOn'};
phases.delay = {'stimOff', 'cue3On'};
phases.init = {'cue3On', 'prodOn'};
phases.iti = {'prodOff', 'cue1On'};
NP.TaskBaseClass.AddEventObjects(ce, phases, 'taskTime');

% Load responsiveness by recordings
recIds = unique(ce.clusTb.recId, 'stable');
sResp = LMV.Resp.LoadPhaseResponseTest(recIds);

% Remove non-responsive units
sigTb = cat(1, sResp.sigTb);
isResp = any(sigTb{:,:}, 2);
ce.RemoveUnits(~isResp);

%% Compute PCA

% Get responses
[t, R] = ce.GetArray('resp', 'Normalization', 'zscore');

% Weigh units inversely to the total number of units in each region
w = zeros(1, size(R,2));
regions = LMV.Param.regions;
for i = 1 : numel(regions)
    isReg = ce.clusTb.region == regions(i);
    w(isReg) = 1/sum(isReg);
end

% Compute PCA
[B, Z, V, ~, varExp] = pca(R, 'Economy', false, 'VariableWeights', w);
B = diag(sqrt(w)) * B;
Z = Z ./ mean(sqrt(w));

% % Compute PCA
% [B, Z, V, ~, varExp] = pca(R, 'Economy', false);

% Add projections to ce
trajTb = table;
trajTb.time = {t};
trajTb.all = {Z(:,1:12)};

for i = 1 : numel(regions)
    % Make region specific projections
    isReg = ce.clusTb.region == regions(i);
    regR = R;
    regR(:,~isReg) = 0;
    regR = regR ./ mean(isReg);
    regZ = regR * B;
    trajTb.(regions(i)) = {regZ(:,1:12)};
end

ce.SetTable('traj', trajTb, 'timeSeries');

%% Plot M2 trajectories projected from all units

trajNames = "all";

f = MPlot.Figure(7711); clf
LMV.Embed.PlotTrajectory(ce, 'Style', '1d', 'DimInd', 1:6, 'TrajectoryNames', trajNames);
MPlot.Paperize(f, 'ColumnsWide', 0.6, 'ColumnsHigh', 1.2);
exportgraphics(f, fullfile(anaDir, "traj_all_1d.png"));

f = MPlot.Figure(7712); clf
LMV.Embed.PlotTrajectory(ce, 'Style', '2d', 'DimInd', 1:3, 'TrajectoryNames', trajNames);
MPlot.Paperize(f, 'ColumnsWide', 1, 'ColumnsHigh', 0.8);
exportgraphics(f, fullfile(anaDir, "traj_all_2d.png"));

f = MPlot.Figure(7713); clf
LMV.Embed.PlotTrajectory(ce, 'Style', '3d', 'DimInd', 1:3, 'TrajectoryNames', trajNames);
MPlot.Paperize(f, 'ColumnsWide', 1, 'ColumnsHigh', 0.8);
exportgraphics(f, fullfile(anaDir, "traj_all_3d.png"));

%% 

trajNames = regions;
trajColors = LMV.Param.GetRegionColors(trajNames);
% trajColors(1,:) = NaN;

f = MPlot.Figure(7721); clf
LMV.Embed.PlotTrajectory(ce, 'Style', '1d', 'DimInd', 1:6, 'TrajectoryNames', trajNames, 'TrajectoryColors', trajColors);
MPlot.Paperize(f, 'ColumnsWide', 0.6, 'ColumnsHigh', 1.2);
exportgraphics(f, fullfile(anaDir, "traj_regions_1d.png"));

f = MPlot.Figure(7722); clf
LMV.Embed.PlotTrajectory(ce, 'Style', '2d', 'DimInd', 1:3, 'TrajectoryNames', trajNames, 'TrajectoryColors', trajColors);
MPlot.Paperize(f, 'ColumnsWide', 1, 'ColumnsHigh', 0.8);
exportgraphics(f, fullfile(anaDir, "traj_regions_2d.png"));

f = MPlot.Figure(7723); clf
LMV.Embed.PlotTrajectory(ce, 'Style', '3d', 'DimInd', 1:3, 'TrajectoryNames', trajNames, 'TrajectoryColors', trajColors);
MPlot.Paperize(f, 'ColumnsWide', 1, 'ColumnsHigh', 0.8);
exportgraphics(f, fullfile(anaDir, "traj_regions_3d.png"));

%% 



