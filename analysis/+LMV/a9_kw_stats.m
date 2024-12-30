%% Plot the fraction of sentence selective units and the strength of selectivity

testList = ["sen14", "sen4"];
verName = join(testList, '-');
anaDir = LMV.Data.GetAnalysisDir("sent_resp", verName);

%% Prepare significance table

% Load selectivity test results
sTest = LMV.Resp.LoadSentenceSelectTest([], testList);
clusTb = sTest.clusTb;
phases = string(sTest.sigTb.Properties.VariableNames);
for i = 1 : numel(phases)
    pn = phases(i);
    clusTb.(pn) = sTest.sigTb.(pn);
    clusTb.(pn+"CS") = cellfun(@(x) x.tbl{2,5}, sTest.clusTb.(pn));
end

% Exclude non-responsive units
rTest = LMV.Resp.LoadPhaseResponseTest();
[~, I] = MMath.SortLike(rTest.clusTb.clusId, sTest.clusTb.clusId);
isNR = ~any(rTest.sigTb{I,:}, 2);
clusTb(isNR,:) = [];

%% Compute fractions by recordings

regions = LMV.Param.regions;
nReg = numel(regions);
nPha = numel(phases);
ss = cell(nPha, nReg);
ff = cell(nPha, nReg);
tt = cell(nPha, nReg);
for i = 1 : nPha
    pn = phases(i);
    for j = 1 : nReg
        isUnit = clusTb.region==regions(j);
        tb = clusTb(isUnit,:);
        
        isSig = tb.(pn) > 0;
        ss{i,j} = tb.(pn+"CS")(isSig);
        tt{i,j} = tb(isSig,:);
        
        G = findgroups(tb.recId);
        ff{i,j} = splitapply(@mean, isSig, G);
    end
end

%% Bar plots of fraction selective

f = MPlot.Figure(9001); clf
tl = tiledlayout("flow");
tl.Padding = "compact";
ax = nexttile;
phaseLabels = ["Atten", "Listen", "Delay", "Init", "Speak"];
LMV.Resp.PlotFractions(ff, 'GroupNames', phaseLabels, 'ChildrenNames', regions, 'Color', LMV.Param.GetRegionColors(regions));
MPlot.Paperize(f, 0.6, 0.4);
% exportgraphics(f, fullfile(anaDir, "frac_selective.png"));
% exportgraphics(f, fullfile(anaDir, "frac_selective.pdf"));

%% Violin plots of Chi-squared

f = MPlot.Figure(9002); clf
tl = tiledlayout("flow");
tl.Padding = "compact";
ax = nexttile;
phaseLabels = ["Atten", "Listen", "Delay", "Init", "Speak"];
LMV.Resp.PlotViolinSets(ss, tt, 'GroupNames', phaseLabels, 'ChildrenNames', regions, 'Color', LMV.Param.GetRegionColors(regions));
% ax.YLim = [0 .5];
ax.YLabel.String = "Chi-squared";
MPlot.Paperize(f);
MPlot.Paperize(f, 0.8, 0.4);
% exportgraphics(f, fullfile(anaDir, "chi-sq_violin.png"));
% exportgraphics(f, fullfile(anaDir, "chi-sq_violin.pdf"));

return
%% Prepare the data to plot Venn diagrams of sentence selectivity in Python

vennDir = LMV.Data.GetAnalysisDir("sent_resp", verName, "venn");

%% Prepare significance table

% Load selectivity test results
sTest = LMV.Resp.LoadSentenceSelectTest([], testList);
clusTb = sTest.clusTb;
phaseNames = sTest.sigTb.Properties.VariableNames;
for i = 1 : numel(phaseNames)
    clusTb.(phaseNames{i}) = sTest.sigTb.(phaseNames{i});
end

% Exclude non-responsive units
rTest = LMV.Resp.LoadPhaseResponseTest();
[~, I] = MMath.SortLike(rTest.clusTb.clusId, sTest.clusTb.clusId);
isNR = ~any(rTest.sigTb{I,:}, 2);
clusTb(isNR,:) = [];

%% Venn diagrams of response types

% Export clusTb for python
colMask = ~varfun(@(x) iscell(x) | isstruct(x) | size(x,2)>1, clusTb, 'OutputFormat', 'uniform');
parquetwrite(fullfile(vennDir, "clusTb.parquet"), clusTb(:,colMask), ...
    'VariableCompression', 'uncompressed', 'VariableEncoding', 'plain');

% Run python scirpt
% \babble\lmv\venn\venn_sentence_selectivity.py
