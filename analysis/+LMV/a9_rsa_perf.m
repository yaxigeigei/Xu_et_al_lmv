%% Single-unit single-trial classification of sentences in each task phase

anaDir = LMV.Data.GetAnalysisDir('coding', LMV.RSA.senMdl);
srcTb = LMV.Data.FindSource([]);

%% Load data

% RSA models
mdlDir = fullfile(anaDir, "computed_mdls");
clusTb = NP.Fit.LoadModels([], mdlDir);
clusTb = cat(1, clusTb{:});

% Task phase responsiveness and sentence selectivity
rTest = LMV.Resp.LoadPhaseResponseTest();
sTest = LMV.Resp.LoadSentenceSelectTest();
[~, I] = MMath.SortLike(rTest.clusTb.clusId, clusTb.clusId);
rSigTb = rTest.sigTb(I,:);
sSigTb = sTest.sigTb(I,:);

%% Extract similarity and cache clusTb

phases = ["atten", "stim", "delay", "init", "prod"];
for i = 1 : numel(phases)
    pn = phases(i);
    clusTb.(pn+"R")(:) = NaN;
    for j = 1 : height(clusTb)
        mdl = clusTb.(pn){j};
        if isempty(mdl)
            continue
        end
        clusTb.(pn+"R")(j) = mdl.r;
        clusTb.(pn+"P")(j) = mdl.pval;
    end
end

cachePath = fullfile(anaDir, "perf_clusTb.mat");
if ~isfile(cachePath)
    save(cachePath, "clusTb", "phases");
else
    fprintf("Not overwriting existing cache:\n%s\n", cachePath);
end

%% Group scores by phase and region

phases = ["atten", "stim", "delay", "init", "prod"];
regions = ["mPrCG", "vPrCG", "IFG", "STG"];
nReg = numel(regions);
nPh = numel(phases);

% Split clusTb of responsive units by regions
cCell = cell(size(regions));
for i = 1 : nReg
    isReg = clusTb.region==regions(i);  % in this region
    isResp = rSigTb.(pn) > 0;           % task responsiveness
    m = isReg & isResp;
    cCell{i} = clusTb(m,:);
end

% Thresholding based on atten
rCell = cell(nPh, nReg);
rTh = cellfun(@(x) prctile(x.attenR, 95), cCell);
rTh(:) = 0;
for i = 1 : nReg
    for j = 1 : nPh
        pn = phases(j);
        r = cCell{i}.(pn+"R");
        p = cCell{i}.(pn+"P");
        m = r > rTh(i) & p < 0.05;
        rCell{j,i} = r(m);
        cCell{i}.(pn+"S") = m;
    end
end

% % Compute recording stats
% fRecCell = cell(size(rCell));
% rRecCell = cell(size(rCell));
% for i = 1 : nPh
%     pn = phases(i);
%     for j = 1 : nReg
%         isSig = cCell{j}.(pn+"S");
%         recIds = cCell{j}.recId;
%         G = findgroups(recIds);
%         fRec = splitapply(@(x) mean(x), isSig, G);
% 
%         r = cCell{j}.(pn+"R");
%         G = findgroups(recIds(isSig));
%         rRec = splitapply(@(x) mean(x), r(isSig), G);
% 
%         fRecCell{i,j} = fRec;
%         rRecCell{i,j} = rRec;
%     end
% end

%% Plot similarity violin scatter

% Prepare data
hasUnit = cellfun(@(x) any(~isnan(x)), rCell);

% Make plot
f = MPlot.Figure(7480); clf
ax = nexttile;

xc = meshgrid(1:nPh, 1:nReg)' + linspace(-.25, .25, nReg);
hh = cell(size(rCell));
xx = cell(size(rCell));
rr = cell(size(rCell));

for i = 1 : nPh
    ys = 0.3;
    for j = 1 : nReg
        if ~hasUnit(i,j)
            continue
        end
        
        % Plot all
        rEdges = 0:0.03:0.5;
        % rEdges = [];
        hh{i,j} = MPlot.ViolinScatter(xc(i,j), rCell{i,j}, rEdges, 'Color', NP.Param.GetRegionColors(regions(j)), 'Span', 0.12);
        hold(ax, 'on');
        x = hh{i,j}.XData';
        r = hh{i,j}.YData'; % get the sorted
        
        % % Plot subthreshold
        % [~, I] = MMath.SortLike(rCell{i,j}, r);
        % p = pCell{i,j}(I);
        % m = r < rTh(j) | p > 0.05;
        % plot(x(m), r(m), '.', 'Color', [0 0 0]+0.7);
        
        xx{i,j} = x;
        rr{i,j} = r; % save the sorted
        
        for k = j+1 : nReg
            ys = ys + 0.03;
            
            % Test significance
            if all(isnan(rCell{i,j})) || all(isnan(rCell{i,k}))
                continue
            end
            p = ranksum(rCell{i,j}, rCell{i,k});
            
            % Plot marker
            if p < 0.05
                c = [0 0 0];
            else
                c = [0 0 0]+.7;
            end
            plot([xc(i,j) xc(i,k)], [ys ys], 'Color', c);
            text(xc(i,end), ys, sprintf("%.1e", p), 'FontSize', 8, 'Color', c);
        end
    end
end

% Set up brushing callback
X = cat(1, xx{:});
R = cat(1, rr{:});
hCB = plot(X, R, 'w.');
hCB.UserData.forBrush = true;
ax.Children = circshift(ax.Children, -1); % move callback plot to the bottom
ax.UserData.clusTb = cat(1, cCell{:});
ax.UserData.taskPhase = "full";
brushObj = brush(gcf);
brushObj.ActionPostCallback = @NP.UnitPlot.RasterOnBrush;

% Formatting
ldg = legend([hh{end,:}], regions);
ldg.Location = 'eastoutside';
ldg.EdgeColor = 'none';

ax.XLim = [0.5 numel(phases)+.5];
ax.XTick = 1:nPh;
ax.XTickLabel = phases;

ax.YLim = [0 .5];
ax.YLabel.String = "Similarity";
MPlot.Axes(ax);

MPlot.Paperize(f, 'ColumnsWide', 1, 'ColumnsHigh', 0.4);
exportgraphics(f, fullfile(anaDir, "sen4_violin_phases.png"));

%% Plot fractions of content encoding neurons

f = MPlot.Figure(8470); clf
tl = tiledlayout("horizontal");
tl.Padding = "compact";
ax = nexttile;

xx = meshgrid(1:nPh, 1:nReg)' + linspace(-.25, .25, nReg);
regionColors = LMV.Param.GetRegionColors(regions);

for i = 1 : numel(phases)
    ys = 0.7;
    for j = 1 : numel(regions)
        % Plot points
        
        plot(xx(i,j), fRecCell{i,j}, 'o', 'Color', regionColors(j,:)); hold on
        
        % Plot significance
        for k = j+1 : numel(regions)
            % Test significance
            p = ranksum(fRecCell{i,j}, fRecCell{i,k});
            
            % Plot marker
            if p < 0.05
                c = [0 0 0];
            else
                c = [0 0 0]+.7;
            end
            plot([xx(i,j) xx(i,k)], [ys ys], 'Color', c);
            text(xx(i,end), ys, sprintf("%.1e", p), 'FontSize', 8, 'Color', c);
            ys = ys + 0.05;
        end
    end
end

ax.XLim = [0.5 5.5];
ax.XTick = 1 : numel(phases);
ax.XTickLabel = ["Attention", "Listening", "Delay", "Initiation", "Speaking"];
ax.YLabel.String = "Frac.";
ax.Title.String = "Fractions of sentence encoding neurons";
MPlot.Axes(ax);
MPlot.Paperize(f, 1, 0.4);
exportgraphics(f, fullfile(anaDir, "frac_by_recordings.png"));

return
%% Plot similarity CDFs

f = MPlot.Figure(7470); clf
tl = tiledlayout(1, numel(phases));
tl.Padding = 'compact';
for i = 1 : numel(phases)
    pn = phases(i);
    ax = nexttile;
    for j = 1 : numel(regions)
        % Get similarity
        r = rCell{i,j};
        
        % % Drop trace if too few samples
        % if sum(~isnan(r)) < 3
        %     r(:) = NaN;
        % end
        
        % Compute CDF
        rEdges = 0:0.002:0.3;
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
    ax.XLabel.String = "Similarity";
    ax.YLabel.String = "Frac. of responsive";
    ax.Title.String = pn;
    MPlot.Axes(ax);
end
MPlot.Paperize(f, 1.6, 0.3);
exportgraphics(f, fullfile(anaDir, "cdf_phases.png"));
exportgraphics(f, fullfile(anaDir, "cdf_phases.pdf"));

%% Compute mean stats on fractions

phases = ["stim", "delay", "init", "prod"];

F = zeros(numel(phases), numel(regions));
rCell = cell(size(F));
for i = 1 : numel(phases)
    pn = phases(i);
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
        
        F(i,j) = mean(r>0);
        rCell{i,j} = r(r>0);
    end
end

fTb = cell2table(F, "RowNames", phases, "VariableNames", regions);
rTb = cell2table(rCell, "RowNames", phases, "VariableNames", regions);

%% 

fileID = fopen(fullfile(anaDir, "perf_stats.txt"), "w");
fprintf(fileID, formattedDisplayText(varfun(@mean, fTb), 'SuppressMarkup', true)+"\n\n");
fprintf(fileID, formattedDisplayText(varfun(@std, fTb), 'SuppressMarkup', true)+"\n\n");
fclose(fileID);

%% 




