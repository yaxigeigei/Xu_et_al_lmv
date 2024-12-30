%% Plot phoneme tuning, sequence response heatmap, PETHs, feature timeseries for an example unit

% recId = 'NP46_B1';
% target = "stim";
% clusId = 460100250;

uidStr = "u"+clusId;
setName = "strf";
mdlName = setName+"_"+target;

srcTb = LMV.Data.FindSource(recId);
figDir = LMV.Data.GetAnalysisDir("peri_event", uidStr+"_"+target);

%% Load data

% Load se
sePath = srcTb.path;
se = NP.SE.LoadSession(sePath);
se = LMV.SE.Transform(se, "enrich");
if string(recId) == "NP54_B1"
    se.RemoveEpochs(72);  % exclude the incomplete last trial of NP54_B1
end
se.SliceSession(0, 'absolute');

clusTb = NP.Unit.GetClusTb(se);
uIdx = find(clusTb.clusId==clusId, 1);

% Load stRF
rfTb = LMV.TRF.LoadModels(mdlName, recId);
mdls = rfTb.(mdlName);
mdl = mdls{uIdx};
mdl = LMV.RF.UpdateModelR2(mdl);

% Load cached sequence data
cacheDir = LMV.Data.GetAnalysisDir("peri_event", "extracted_seq");
s = load(fullfile(cacheDir, recId+"_seqData.mat"), "se", "feats", "seqData");

%% Plot TRF weights as a heatmap

f = MPlot.Figure(985); clf
ax = nexttile;
LMV.TRF.PlotWeights2(mdl, 'ClusInfo', clusTb(uIdx,:), 'Parent', ax);
MPlot.Axes(ax);
MPlot.Paperize(f, 'ColumnsWide', 0.3, 'ColumnsHigh', .8);
exportgraphics(f, fullfile(figDir, uidStr+"_"+setName+"_stRF.png"));
print(f, fullfile(figDir, uidStr+"_"+setName+"_stRF.pdf"), '-dpdf');

%% Plot sliding time r2

f = MPlot.Figure(956); clf
LMV.RF.PlotR2Timecourse(mdl);
MPlot.Paperize(f, 'ColumnsWide', 0.35, 'ColumnsHigh', 0.3);
exportgraphics(f, fullfile(figDir, uidStr+"_"+setName+"_r2_timecourse.png"));
print(f, fullfile(figDir, uidStr+"_"+setName+"_r2_timecourse.pdf"), '-dpdf');

%% Plot RF weights at optimal time

f = MPlot.Figure(986); clf
ax = nexttile;

rfMdl = mdl.mdls{mdl.r2Idx};
rfMdl.feats = se.userData.melMeta.F/1000;
rfMdl.resps = mdl.resps;
% rfMdl.Beta = flip(rfMdl.Beta);
% rfMdl.Beta = flip(rfMdl.Beta);

LMV.RF.PlotWeights(rfMdl, 'Parent', ax, 'Orientation', 'vertical');

ax.YDir = "normal";
ax.YTick = linspace(1, numel(rfMdl.feats), 5);
ax.YTickLabel = round(interp1(1:numel(rfMdl.feats), rfMdl.feats, ax.YTick), 2);
ax.YTickLabelRotation = 90;
ax.YLabel.String = "Frequency (Hz)";
MPlot.Axes(ax);
MPlot.Paperize(f, 'ColumnsWide', .3, 'ColumnsHigh', .8);
exportgraphics(f, fullfile(figDir, uidStr+"_"+setName+"_opti-RF.png"));
print(f, fullfile(figDir, uidStr+"_"+setName+"_opti-RF.pdf"), '-dpdf');

%% Extract responses and features from se

% Find the top features
[~, I] = maxk(abs(rfMdl.Beta), 1);
featNames = mdl.feats(I);

% Define time window
switch target
    case 'stim'
        tWin = [-0.1 0.4];
        phaseName = 'stim';
    case 'prod'
        tWin = [-0.4 0.1];
        phaseName = 'prod';
    case 'feedback'
        tWin = [-0.1 0.4];
        phaseName = 'prod';
    otherwise
        error("'%s' is not a supported target name.", target);
end

C = s.seqData.(target){"phone"}.seq1;

for k = 1 : numel(C)
    seqTb = C{k}.seqTb;
    if isempty(seqTb)
        continue
    end
    seqTb(seqTb.nSample<3,:) = [];
    seqTb = LMV.PE.AddResp2SeqTb(seqTb, se, tWin, uIdx, mdls);
    seqTb = LMV.PE.AddFeat2SeqTb(seqTb, se, flip(-tWin), "mel", "speaker1");
    C{k}.seqTb = seqTb;
end

seqTb = cellfun(@(x) x.seqTb, C, 'Uni', false);
seqTb = cat(1, seqTb{:});

%% Plot peri-event raster and activation

f = MPlot.Figure(235); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";

LMV.Fig.PlotSpectralResp(mdl, seqTb, I, se.userData.melMeta.F(I));

MPlot.Paperize(f, 'ColumnsWide', .8, 'ColumnsHigh', .3);
exportgraphics(f, fullfile(figDir, uidStr+"_"+setName+"_resp.png"));
print('-vector', f, fullfile(figDir, uidStr+"_"+setName+"_resp.pdf"), '-dpdf');
