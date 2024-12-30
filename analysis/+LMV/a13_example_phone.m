%% Plot phoneme tuning, sequence response heatmap, PETHs, feature timeseries for an example unit

%{

recId = 'NP41_B1';
target = "prod";
clusId = 410100245;

recId = 'NP44_B3';
target = "stim";
clusId = 440300585;

%}

uidStr = "u"+clusId;
figDir = LMV.Data.GetAnalysisDir("peri_event", uidStr+"_"+target);

%% Load data

% Load cached sequence data
cacheDir = LMV.Data.GetAnalysisDir("peri_event", "extracted_seq");
s = load(fullfile(cacheDir, recId+"_seqData.mat"), "se", "feats", "seqData");
seeds = s.feats;
se = s.se;
clusTb = NP.Unit.GetClusTb(se);

% Load phone stRFs
trfTb = LMV.TRF.LoadModels("phone_"+target, recId);
[~, trfInd] = MMath.SortLike(trfTb.clusId, clusTb.clusId);
mdls = trfTb.("phone_"+target)(trfInd);
mdls = LMV.RF.UpdateModelR2(mdls);

%% 

% Extract responses and features from se
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

Cp = s.seqData.(target){"phone"}{:,:};

for k = 1 : numel(Cp)
    seqTb = Cp{k}.seqTb;
    if isempty(seqTb)
        continue
    end
    seqTb(seqTb.nSample<3,:) = [];
    seqTb = LMV.PE.AddResp2SeqTb(seqTb, se, tWin, [], mdls);
    Cp{k}.seqTb = seqTb;
end

% Mask out sequences with no further branching
nSeq = cellfun(@(x) height(x.seqTb), Cp);
for k = size(nSeq,2) : -1 : 2
    m = nSeq(:,k) <= nSeq(:,k-1);
    nSeq(m,k) = NaN;
end

% Get the optimal encoding time
uIdx = find(clusTb.clusId==clusId, 1);
mdl = mdls{uIdx};
tResp = Cp{1}.seqTb.tResp{1};
[~, iOpt] = min(abs(tResp - mdl.r2t));

% Find the spike rate at optimal encoding time
R = NaN(size(Cp));
for i = 1 : numel(Cp)
    seqTb = Cp{i}.seqTb;
    if isempty(seqTb)
        continue
    end
    R(i) = max(cellfun(@(y) y(iOpt,uIdx), seqTb.hh));
end
R(isnan(nSeq)) = NaN;

% Find the max response for each seed and sort seeds in descending order
maxR = max(R, [] ,2);
[maxR, iR] = sort(maxR, 'descend');

%% Plot sliding time r2

f = MPlot.Figure(956); clf
LMV.RF.PlotR2Timecourse(mdl);
MPlot.Paperize(f, 'ColumnsWide', 0.35, 'ColumnsHigh', 0.3);
exportgraphics(f, fullfile(figDir, uidStr+"_phone_r2_timecourse.png"));
print(f, fullfile(figDir, uidStr+"_phone_r2_timecourse.pdf"), '-dpdf');

%% Plot RF weights at optimal time

f = MPlot.Figure(986); clf
ax = nexttile;
rfMdl = mdl.mdls{mdl.r2Idx};
rfMdl.feats = mdl.feats;
LMV.RF.PlotWeights(rfMdl, 'ClusInfo', clusTb(uIdx,:), 'Parent', ax, 'Orientation', 'vertical');
ax.YTickLabel = MLing.ARPA2IPA(ax.YTickLabel);
MPlot.Axes(ax);
MPlot.Paperize(f, 'ColumnsWide', .2, 'ColumnsHigh', .8);
exportgraphics(f, fullfile(figDir, uidStr+"_phone_opti-RF.png"));
print(f, fullfile(figDir, uidStr+"_phone_opti-RF.pdf"), '-dpdf');

%% Plot heatmap of maximum sequence responses

f = MPlot.Figure(987); clf
h = heatmap("n+"+(1:3), MLing.ARPA2IPA(seeds(iR)), round(R(iR,:)));
h.Title = uidStr;
h.Colormap = colormap('copper');
MPlot.Paperize(f, 'ColumnsWide', 0.6, 'ColumnsHigh', 1);
exportgraphics(f, fullfile(figDir, uidStr+"_seq_tuning.png"));

%% Plot top sequence responses

for i = 1 : min(10, numel(iR))
    k = iR(i);
    s2plot = Cp(k,:);
    nBranch = cellfun(@(x) height(x.seqTb), s2plot);
    
    figName = sprintf("%s_%i_%s.png", uidStr, i, seeds(k));
    figPath = fullfile(figDir, figName);
    % if exist(figPath, 'file')
    %     continue
    % end
    
    f = MPlot.Figure(234); clf
    LMV.Fig.PlotPeriEventSeqArray(s2plot, 'UnitIdx', uIdx);
    MPlot.Paperize(f, 'ColumnsWide', .8, 'ColumnsHigh', 0.1+0.07*max(nBranch)); return
    exportgraphics(f, figPath);
    print(f, strrep(figPath, ".png", ".pdf"), '-dpdf');
end

%% Plot top phoneme rasters

switch clusId
    case 410100245
        [~, ind] = MMath.SortLike(seeds, ["P", "JH"]);
    case 440300585
        [~, ind] = MMath.SortLike(seeds, ["AE", "S", "Z", "HH"]);
    otherwise
        return
end

s2plot = Cp{1};
s2plot.seed = "";
seqTbs = cellfun(@(x) x.seqTb, Cp(ind,1), 'Uni', false);
s2plot.seqTb = cat(1, seqTbs{:});

f = MPlot.Figure(235); clf
LMV.Fig.PlotSeqResp(nexttile, s2plot, 'UnitIdx', uIdx, 'Sort', false);
MPlot.Paperize(f, 0.25, 0.55);
figName = sprintf("%s_%s.png", uidStr, strjoin(seeds(ind), "_"));
figPath = fullfile(figDir, figName);
% return
exportgraphics(f, figPath);
print(f, strrep(figPath, ".png", ".pdf"), '-dpdf');

%% Plot specific sequences

switch clusId
    case 410100245
        k = find(seeds=="P", 1);
        h = [0.4 0.6];
    case 440300585
        k = find(seeds=="S", 1);
        h = [0.4 0.6 1.2];
    otherwise
        return
end
s2plot = Cp(k,:);

f = MPlot.Figure(234); clf
LMV.Fig.PlotPeriEventSeqArray(s2plot, 'UnitIdx', uIdx);

for i = 1 : numel(h)
    MPlot.Paperize(f, 0.8, h(i));
    figName = sprintf("%s_%s_scale%i.png", uidStr, seeds(k), i);
    figPath = fullfile(figDir, figName);
    return
    exportgraphics(f, figPath);
    print(f, strrep(figPath, ".png", ".pdf"), '-dpdf');
end