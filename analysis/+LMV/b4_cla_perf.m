%% 

anaDir = LMV.Data.GetAnalysisDir('coding', 'phase');
srcTb = LMV.Data.FindSource('decode');


%% Load task phase decoding models

mdlDir = fullfile(anaDir, 'computed_phase_mdls');
mdlName = 'phase';
mdlTb = LMV.Decode.LoadClassModels(mdlDir, '*', mdlName);

%% 

f = MPlot.Figure(5500); clf
LMV.Decode.PlotClassAccuracy(mdlTb);
MPlot.Paperize(f, 'ColumnsWide', 0.5, 'AspectRatio', 4);
figPath = fullfile(anaDir, sprintf('%s_cla_accuracy.png', mdlName));
exportgraphics(f, figPath);

%% 

f = MPlot.Figure(5501); clf
LMV.Decode.PlotConfusionMat(mdlTb);
MPlot.Paperize(f, 'ColumnsWide', 0.8, 'AspectRatio', 2);
figPath = fullfile(anaDir, sprintf('%s_cla_confusion.png', mdlName));
exportgraphics(f, figPath);


%% Load 4-trial decoding models

mdlGroup = 'trial4-ovo';
mdlDir = fullfile(anaDir, ['computed_' mdlGroup '_mdls']);
mdlTbs = cellfun(@(x) LMV.Decode.LoadClassModels(mdlDir, x, '*'), srcTb.recId, 'Uni', false);

for i = 1 : numel(mdlTbs)
    mdlCats = categorical(mdlTbs{i}.mdlName, {'atten', 'stim', 'delay', 'init', 'prod'} + "4");
    [~, I] = sort(mdlCats);
    mdlTbs{i} = mdlTbs{i}(I,:);
end

%% 

f = MPlot.Figure(5600); clf
LMV.Decode.PlotGrossAccuracy(mdlTbs);
MPlot.Paperize(f, 'ColumnsWide', 0.5, 'AspectRatio', 4);
figPath = fullfile(anaDir, sprintf('%s_cla_accuracy.png', mdlGroup));
exportgraphics(f, figPath);

%% 

f = MPlot.Figure(5601);
t = tiledlayout(numel(mdlTbs), 5);
t.Padding = 'compact';
for i = 1 : numel(mdlTbs)
    LMV.Decode.PlotConfusionMat(mdlTbs{i}, ClassNames="S"+[3 4 8 13], Layout=t, SigOnly=true);
end
MPlot.Paperize(f, 'ColumnsWide', 1.3, 'ColumnsHigh', 2);
figPath = fullfile(anaDir, sprintf('%s_cla_confusion.png', mdlGroup));
exportgraphics(f, figPath);

%% 

f = MPlot.Figure(5602);
t = tiledlayout(numel(mdlTbs), 5);
t.Padding = 'compact';
for i = 1 : numel(mdlTbs)
    LMV.Decode.PlotConfusionMat(mdlTbs{i}, ClassNames="S"+[3 4 8 13], Layout=t, Diff=true);
end
MPlot.Paperize(f, 'ColumnsWide', 1.3, 'ColumnsHigh', 2);
figPath = fullfile(anaDir, sprintf('%s_cla_confusion_diff.png', mdlGroup));
% exportgraphics(f, figPath);


%% Load 14-trial decoding models

mdlGroup = 'trial14';
mdlDir = fullfile(anaDir, ['computed_' mdlGroup '_mdls']);
mdlTbs = cellfun(@(x) LMV.Decode.LoadClassModels(mdlDir, x, '*'), srcTb.recId, 'Uni', false);

for i = 1 : numel(mdlTbs)
    mdlCats = categorical(mdlTbs{i}.mdlName, {'atten', 'stim', 'delay', 'init', 'prod'} + "14");
    [~, I] = sort(mdlCats);
    mdlTbs{i} = mdlTbs{i}(I,:);
end

%% 

f = MPlot.Figure(5701); clf
LMV.Decode.PlotGrossAccuracy(mdlTbs);
MPlot.Paperize(f, 'ColumnsWide', 0.5, 'AspectRatio', 4);
figPath = fullfile(anaDir, sprintf('%s_cla_accuracy.png', mdlGroup));
exportgraphics(f, figPath);

