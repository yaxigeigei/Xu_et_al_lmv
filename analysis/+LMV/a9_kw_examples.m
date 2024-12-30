%% Single-unit classification of sentences at each task phase

testList = ["sen14", "sen4"];
verName = join(testList, '-');
anaDir = LMV.Data.GetAnalysisDir("sent_resp", verName);

%% Load test results

% Load selectivity test results
sTest = LMV.Resp.LoadSentenceSelectTest([], testList);
clusTb = sTest.clusTb;

%% Find most strongly selective units in each region

phases = ["stim", "delay", "init", "prod"];
cs = cellfun(@(x) x.tbl{2,5}, clusTb{:,phases});
p = cellfun(@(x) x.p, clusTb{:,phases});
cs(p > 0.05) = NaN;
[clusTb.maxCS, I] = max(cs, [], 2, "omitmissing");
clusTb.maxPhase = phases(I)';
clusTb.maxPhase(isnan(clusTb.maxCS)) = "none";

regions = LMV.Param.regions;
nR = numel(regions);
nE = 5;
cids = zeros(nR, nE);
pns = strings(size(cids));
for i = 1 : nR
    ind = find(clusTb.region==regions(i) & ~isnan(clusTb.maxCS));
    cs = clusTb.maxCS(ind);
    [~, I] = maxk(cs, size(cids,2));
    cids(i,:) = clusTb.clusId(ind(I));
    pns(i,:) = clusTb.maxPhase(ind(I));
end

%% Plot 14-sentence rasters

exampleDir = LMV.Data.GetAnalysisDir("sent_resp", verName);
stimIdList = [LMV.Param.stimIdList4; setdiff(LMV.Param.stimIdList, LMV.Param.stimIdList4)];

f = MPlot.Figure(234); clf
tl = tiledlayout(nE, nR);
tl.Padding = "compact";
for i = 1 : nE
    for j = 1 : nR
        ax = nexttile;
        LMV.Fig.SessionRaster(ax, cids(j,i), 'DataSource', 'm2', 'StimIdList', stimIdList, 'PlotFun', 'rate_map');
        title(ax, sprintf('%s: u%i; %s', regions(j), cids(j,i), pns(j,i)));
        % return
    end
end

MPlot.Paperize(f, 2, 1.4);
figName = sprintf("top_session_rasters_s12");
exportgraphics(f, fullfile(anaDir, figName+".png"));
exportgraphics(f, fullfile(anaDir, figName+".pdf"), ContentType="vector");

return
%% 

anaDir = LMV.Data.GetAnalysisDir('coding', LMV.UnitDecode.senMdl);
srcTb = LMV.Data.FindSource([]);

%% Load data

load(fullfile(anaDir, "perf_clusTb.mat"), "clusTb", "phaseNames");

%% Choose examples to plot

uu = LMV.Fig.GetExampleInfo("Fig1G");
uid = uu.clusId;

%% Plot 14-sentence rasters

f = MPlot.Figure(234); clf
stimIdList = [LMV.Param.stimIdList4; setdiff(LMV.Param.stimIdList, LMV.Param.stimIdList4)];
LMV.Fig.SessionFromCache(uid, 'DataSource', 'm2', 'StimIdList', stimIdList, 'PlotFun', 'rate_map');
MPlot.Paperize(f, 2.5, .7);
exportgraphics(f, fullfile(anaDir, "example_session_rasters_s12.png"));

%% Plot 4-sentence rasters

uid1 = uid(end);

f = MPlot.Figure(345); clf
LMV.Fig.SessionFromCache(uid1, 'DataSource', 'm2', 'StimIdList', LMV.Param.stimIdList4, 'PlotFun', 'rate_map');
MPlot.Paperize(f, .8, 0.3);
exportgraphics(f, fullfile(anaDir, "u"+uid1+"_session_rasters_s4.png"));
exportgraphics(f, fullfile(anaDir, "u"+uid1+"_session_rasters_s4.pdf"), ContentType="vector");

%% Plot classification accuracy of example units in a heatmap

R = NaN(numel(uid), numel(phaseNames));
for i = 1 : numel(uid)
    for j = 1 : numel(phaseNames)
        pn = phaseNames(j);
        k = clusTb.clusId==uid(i);
        R(i,j) = clusTb.(pn+"R")(k);
    end
end
R(R==1 | R<0.25) = NaN;
R = round(R, 2);

f = MPlot.Figure(5802); clf
h = heatmap(phaseNames, "u"+uid, R);
% h.Colormap = gray;
h.ColorbarVisible = 'off';
h.FontSize = 12;
MPlot.Paperize(f, 'ColumnsWide', 0.6, 'ColumnsHigh', 0.4);
exportgraphics(f, fullfile(anaDir, "example_r.png"));

%% Find units from each region with top performance

phases = ["stim", "delay", "init", "prod"];
regions = LMV.Param.regions;
nr = numel(regions);
uuu = cell(1, nr);
for i = 1 : nr
    R = clusTb{:,phases+"R"};
    R(clusTb.region~=regions(i),:) = NaN;
    [rMax, I] = maxk(R(:), 20);
    [a, b] = ind2sub(size(R), I);
    tb = table;
    tb.clusId = clusTb.clusId(a);
    tb.meanActiveRate = clusTb.meanActiveRate(a);
    tb.phase = phases(b)';
    tb.r = rMax;
    uuu{i} = tb;
end

%% Plot 4-sentence rasters

for i = 1 : nr
    f = MPlot.Figure(440); clf
    LMV.Fig.SessionFromCache(uuu{i}.clusId, 'DataSource', 'm2', 'StimIdList', LMV.Param.stimIdList4, 'NRows', 4, 'PlotFun', 'rate_map');
    MPlot.Paperize(f, 2.5, 1.2);
    exportgraphics(f, fullfile(anaDir, sprintf("example_best_perf_%s.png", regions(i))));
end

return
%% Plot sentence overview

f = MPlot.Figure(123); clf
LMV.Overview.SentencesFromCache(uid, 'StimIds', LMV.Param.stimIdList4, 'NColumns', 2, 'Features', "phone", 'UnitLabels', "depth", 'MinTaskHight', 12);
MPlot.Paperize(f, 'ColumnsWide', 2.2, 'ColumnsHigh', 1.4, 'FontSize', 5);
exportgraphics(f, fullfile(anaDir, "example_sentence_rasters_s4.png"));

%% Plot confusion matrices of example units

k = find(clusTb.clusId==uid(4)); pn = 'prod'; % u245
% k = find(clusTb.clusId==cid(11)); pn = 'init'; % u441
% k = find(clusTb.clusId==cid(21)); pn = 'init'; % u463

for i = 1 : numel(phaseNames)
    fprintf("u%i %s: %.2f\n", clusTb.clusId(k), phaseNames(i), 1-clusTb.(phaseNames(i)+"R")(k));
end

f = MPlot.Figure(2); clf
for i = 1 : numel(uid)
    k = find(clusTb.clusId==uid(4));
    if isempty(clusTb.(pn){k})
        continue
    end
    ax = nexttile;
    confusionchart(clusTb.(pn){k}.confMat);
    title("u"+clusTb.clusId(k)+" "+pn);
end

%% 

% [~, cid] = LMV.Param.GetSelectedClusId('NP41_B1');
[~, uid] = LMV.Param.GetSelectedClusId('NP44_B3');

figure(456); clf
LMV.Overview.SentencesFromCache(uid);
