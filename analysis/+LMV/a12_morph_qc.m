%% Quanlity control of time morphing by examining and quantifying aligned spectrograms

anaDir = LMV.Data.GetAnalysisDir("linker");

%% Load data

ceDir = fullfile(anaDir, 'computed_peth_sem-discounted');
seTbDir = fullfile(anaDir, 'computed_seTb');

ceSearch = MBrowse.Dir2Table(fullfile(ceDir, '*_ce.mat'));
seTbSearch = MBrowse.Dir2Table(fullfile(seTbDir, '*_seTb.mat'));
clear ceArray seTbs
for k = height(ceSearch) : -1 : 1
    load(fullfile(ceSearch.folder{k}, ceSearch.name{k}), 'ce');
    load(fullfile(seTbSearch.folder{k}, seTbSearch.name{k}), 'seTb');
    ceArray(k) = ce;
    seTbs{k} = seTb;
end

%% Inspect alignment with single-trial Mel-spectrograms

figDir = fullfile(anaDir, 'alignment');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

for k = 2 : numel(ceArray)
    seTb = seTbs{k};
    ce = ceArray(k);
    
    seTb = sortrows(seTb, {'numTrial', 'stimText'}, {'descend', 'ascend'});
    m = 1:8;
    
    f = MPlot.Figure(6190); clf
    % f.WindowState = 'maximized';
    
    LMV.Linker.PlotMelAlignment(seTb(m,:), ce);
    
    MPlot.Paperize(f, 2.5, 1.5); %return
    figPath = fullfile(figDir, sprintf("aligned4_%s", NP.SE.GetID(ce)));
    exportgraphics(f, figPath+".png");
    print(f, '-vector', figPath+".pdf", '-dpdf');
end

%% Quantify the overall alignment

% Compute mean across recordings
xrTb = struct2table(arrayfun(@(x) x.userData.xrMel, ceArray));
dt = xrTb.xrLag(1,:);
[rMean, ~, ~, rCI] = MMath.MeanStats(xrTb.xrCoef, 1, 'Alpha', 0.05);

% Flip time
dt = flip(-dt);
rMean = flip(rMean);
rCI = flip(rCI, 2);

% Find peak
[rPk, I] = max(rMean);
xPk = dt(I);

% Plot x-correlogram
f = MPlot.Figure(6191); clf
MPlot.ErrorShade(dt, rMean, rCI(1,:), rCI(2,:), 'IsRelative', false); hold on
plot(dt, rMean, 'k');
plot([xPk xPk], [-.2 rPk], 'k--');
xlim(dt([1 end]));
xlabel("Time shift in production (s)");
ylabel("Pearson's r");
title("Mean x-correlogram of Mel power");
MPlot.Paperize(f, .5, .4);
exportgraphics(f, fullfile(anaDir, "mean_mel_xcorr.png"));
exportgraphics(f, fullfile(anaDir, "mean_mel_xcorr.pdf"));
