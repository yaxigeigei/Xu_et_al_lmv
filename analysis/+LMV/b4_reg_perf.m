%% 

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'decode');

%% Load speech decoding models

mdlGroup = 'speech';
mdlDir = fullfile(anaDir, ['computed_' mdlGroup '_mdls']);
mdlNames = {'prod', 'stim'};
mdls = LMV.Decode.LoadRegModels(mdlDir, mdlNames);
mdlsHS = LMV.Decode.LoadRegModels(mdlDir, mdlNames+"_hs");

%% Plot hyperparam search

f = MPlot.Figure(885201); clf
f.WindowState = 'maximized';
LMV.Decode.PlotOptimalHyperparam(mdlsHS, 'LambdaMinMSE');
figPath = fullfile(anaDir, sprintf('%s_reg_hyper_search_lambda.png', mdlGroup));
% exportgraphics(f, figPath);

%% Plot goodness-of-fit

f = MPlot.Figure(885301); clf
f.WindowState = 'maximized';
LMV.Decode.PlotRegPerf(mdls);
figPath = fullfile(anaDir, sprintf('%s_reg_r.png', mdlGroup));
exportgraphics(f, figPath);
