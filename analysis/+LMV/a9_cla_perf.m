%% Single-unit single-trial classification of sentences in each task phase

anaDir = LMV.Data.GetAnalysisDir('coding', LMV.UnitDecode.senMdl);
srcTb = LMV.Data.FindSource([]);

%% Load data

% Classifiers
mdlDir = fullfile(anaDir, "computed_mdls");
clusTb = NP.Fit.LoadModels([], mdlDir);
clusTb = cat(1, clusTb{:});

% Task phase responsiveness and sentence selectivity
rTest = LMV.Resp.LoadPhaseResponseTest();
sTest = LMV.Resp.LoadSentenceSelectTest();
[~, I] = MMath.SortLike(rTest.clusTb.clusId, clusTb.clusId);
rSigTb = rTest.sigTb(I,:);
sSigTb = sTest.sigTb(I,:);

%% Extract accuracies and cache clusTb

phaseNames = ["atten", "stim", "delay", "init", "prod"];
for i = 1 : numel(phaseNames)
    pn = phaseNames(i);
    clusTb.(pn+"R")(:) = NaN;
    for j = 1 : height(clusTb)
        mdl = clusTb.(pn){j};
        if isempty(mdl)
            continue
        end
        clusTb.(pn+"R")(j) = 1 - mdl.loss;
    end
end

cachePath = fullfile(anaDir, "perf_clusTb.mat");
if ~isfile(cachePath)
    save(cachePath, "clusTb", "phaseNames");
else
    fprintf("Not overwriting existing cache:\n%s\n", cachePath);
end

%% Plot classification accuracy CDFs

phaseNames = ["atten", "stim", "delay", "init", "prod"];
regions = ["mPrCG", "vPrCG", "IFG", "STG"];

f = MPlot.Figure(7470); clf
tl = tiledlayout(1, numel(phaseNames));
tl.Padding = 'compact';
for i = 1 : numel(phaseNames)
    pn = phaseNames(i);
    ax = nexttile;
    for j = 1 : numel(regions)
        % Get classification accuracy
        r = clusTb.(pn+"R");
        
        % Set accuracy of non-selective units to NaN
        r(sSigTb.(pn)==0) = NaN;
        
        % Select units
        isReg = clusTb.region==regions(j);
        isResp = rSigTb.(pn) > 0; % task responsiveness
        m = isReg & isResp;
        r = r(m);
        
        % Exclude below chance level units
        rTh = 0.25;
        r(r <= rTh) = NaN;
        
        % Drop trace if too few samples
        if sum(~isnan(r)) < 3
            r(:) = NaN;
        end
        
        % Compute CDF
        rEdges = rTh:0.005:1;
        rCenters = MMath.BinEdges2Centers(rEdges);
        N = histcounts(r, rEdges, 'Normalization', 'cdf');
        
        % Plot CDF
        stairs(rCenters, N, 'Color', NP.Param.GetRegionColors(regions(j)));
        hold(ax, 'on');
    end
    if pn == "atten"
        ldg = legend(ax, regions);
        ldg.Location = 'southeast';
        ldg.EdgeColor = 'none';
    end
%     ax.XLim = [rTh 0.9];
%     ax.YLim = [0 0.3];
    ax.XLabel.String = "Accuracy";
    ax.YLabel.String = "Frac. of responsive";
    ax.Title.String = pn;
    MPlot.Axes(ax);
end
MPlot.Paperize(f, 1.6, 0.3);
exportgraphics(f, fullfile(anaDir, "sen4_cdf_phases.png"));
exportgraphics(f, fullfile(anaDir, "sen4_cdf_phases.pdf"));

%% Compute mean stats on fractions

phaseNames = ["stim", "delay", "init", "prod"];

F = zeros(numel(phaseNames), numel(regions));
rCell = cell(size(F));
for i = 1 : numel(phaseNames)
    pn = phaseNames(i);
    for j = 1 : numel(regions)
        % Get classification accuracy
        r = clusTb.(pn+"R");
        
        % Set accuracy of non-selective units to NaN
        r(sSigTb.(pn)==0) = NaN;
        
        % Select units
        isReg = clusTb.region==regions(j);
        isResp = rSigTb.(pn) > 0; % task responsiveness
        m = isReg & isResp;
        r = r(m);
        
        F(i,j) = mean(r>0.25);
        rCell{i,j} = r(r>0.25);
    end
end

fTb = array2table(F, "RowNames", phaseNames, "VariableNames", regions);
rTb = cell2table(rCell, "RowNames", phaseNames, "VariableNames", regions);

%% 

fileID = fopen(fullfile(anaDir, "perf_stats.txt"), "w");
fprintf(fileID, formattedDisplayText(varfun(@mean, fTb), 'SuppressMarkup', true)+"\n\n");
fprintf(fileID, formattedDisplayText(varfun(@std, fTb), 'SuppressMarkup', true)+"\n\n");
fclose(fileID);

%% 

phaseNames = ["stim", "delay", "init", "prod"];

% Get fractions
fCell = cell(numel(phaseNames), numel(regions));
rCell = cell(size(fCell));
for i = 1 : numel(phaseNames)
    pn = phaseNames(i);
    for j = 1 : numel(regions)
        % Get classification accuracy
        r = clusTb.(pn+"R");
        
        % Set accuracy of non-selective units to NaN
        r(sSigTb.(pn)==0) = NaN;
        
        % Select units
        isReg = clusTb.region==regions(j);
        isResp = rSigTb.(pn) > 0; % task responsiveness
        m = isReg & isResp;
        r = r(m);
        
        % Compute fraction for each recording
        G = findgroups(clusTb.recId(m));
        fRec = splitapply(@(x) mean(x>0.25), r, G);
        rRec = splitapply(@(x) median(x(x>0.25)), r, G);
        
        fCell{i,j} = fRec;
        rCell{i,j} = rRec;
    end
end

% % Compute mean stats
% [fMean, fSD] = cellfun(@MMath.MeanStats, fCell);

% Test difference in fractions between pairwise regions
P = NaN(numel(regions), numel(regions), numel(phaseNames));
for i = 1 : numel(phaseNames)
    for j = 1 : numel(regions)
        for k = j+1 : numel(regions)
            p = ranksum(fCell{i,j}, fCell{i,k});
            P(j,k,i) = p;
        end
    end
end

f = MPlot.Figure(8470); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";
ax = nexttile;
xx = meshgrid(1:4)' + linspace(-.25, .25, 4);
regionColors = LMV.Param.GetRegionColors(regions);

for i = 1 : numel(regions)
    yy = cat(2, fCell{:,i})';
    plot(xx(:,i), yy, 'o', 'Color', regionColors(i,:)); hold on
end

for i = 1 : numel(phaseNames)
    y = 0.6;
    for j = 1 : numel(regions)
        for k = j+1 : numel(regions)
            p = P(j,k,i);
            if p < 0.05
                c = [0 0 0];
            else
                c = [0 0 0]+.7;
            end
            plot([xx(i,j) xx(i,k)], [y y], 'Color', c);
            text(xx(i,end), y, sprintf("%.1e", p), 'FontSize', 8, 'Color', c);
            y = y + 0.05;
        end
    end
end

% hh = errorbar(xx, fMean, fSD);
% for k = 1 : numel(hh)
%     hh(k).Color = regionColors(k,:);
%     hh(k).LineStyle = "none";
% end

ax.XLim = [0.5 4.5];
ax.XTick = 1 : numel(phaseNames);
% ax.YLim(2) = [];
ax.XTickLabel = ["Listening", "Delay", "Initiation", "Speaking"];
MPlot.Axes(ax);
MPlot.Paperize(f, 0.7, 0.4);

%% Plot classification accuracy violin scatter

% Prepare data
nReg = numel(regions);
nPh = numel(phaseNames);
rr = cell(nPh, nReg);
tt = cell(nPh, nReg);
for i = 1 : nPh
    pn = phaseNames(i);
    for j = 1 : nReg
        % Overall accuracy (only responsive units have decoding models)
        r = NaN(size(clusTb.(pn)));
        isMdl = ~cellfun(@isempty, clusTb.(pn));
        r(isMdl) = cellfun(@(x) 1-x.loss, clusTb.(pn)(isMdl));
        clusTb.(pn+"R") = r;
        
        % Set accuracy of non-selective units to NaN
        r(sSigTb.(pn)==0) = NaN;
        
        % Select units
        isReg = clusTb.region==regions(j);
        isResp = rSigTb.(pn) > 0; % task responsiveness
        
        % Exclude below chance level units
        rTh = 0.25;
        isChance = r <= 0.25; % below chance level
        
        m = isReg & isResp & ~isChance & ~isnan(r);
        r = r(m);
        ctb = clusTb(m,:);
        
        rr{i,j} = r;
        tt{i,j} = ctb;
    end
end
hasUnit = cellfun(@(x) any(~isnan(x)), rr);

% Make plot
f = MPlot.Figure(7480); clf
ax = nexttile;
hh = cell(size(rr));
xx = cell(size(rr));
phSp = nReg + 2;
for i = 1 : numel(phaseNames)
    pn = phaseNames(i);
    for j = 1 : nReg
        if ~hasUnit(i,j)
            continue
        end
        x = (i-1)*phSp + j;
        rEdges = rTh:0.05:1;
        hh{i,j} = MPlot.ViolinScatter(x, rr{i,j}, rEdges, 'Color', NP.Param.GetRegionColors(regions(j)));
        hold(ax, 'on');
        xx{i,j} = hh{i,j}.XData';
        rr{i,j} = hh{i,j}.YData'; % get the sorted
    end
end

% Set up brushing callback
X = cat(1, xx{:});
R = cat(1, rr{:});
hCB = plot(X, R, 'w.');
hCB.UserData.forBrush = true;
ax.Children = circshift(ax.Children, -1); % move callback plot to the bottom
ax.UserData.clusTb = cat(1, tt{:});
ax.UserData.taskPhase = "full";
brushObj = brush(gcf);
brushObj.ActionPostCallback = @NP.UnitPlot.RasterOnBrush;

% Formatting
ldg = legend([hh{end,:}], regions);
ldg.Location = 'eastoutside';
ldg.EdgeColor = 'none';

ax.XLim = [0 numel(phaseNames)*phSp];
ax.XTick = ((1:numel(phaseNames))-1)*phSp + 2.5;
ax.XTickLabel = phaseNames;

ax.YLim = [rTh 1];
ax.YLabel.String = "Accuracy";
MPlot.Axes(ax);

MPlot.Paperize(f, 'ColumnsWide', 1, 'ColumnsHigh', 0.3);
exportgraphics(f, fullfile(anaDir, "sen4_violin_phases.png"));
