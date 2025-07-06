%% Analyze chunking patterns in speech production using morphed speed data

anaDir = LMV.Data.GetAnalysisDir("chunking");

srcTb = LMV.Data.FindSource([]);

%% 

recIds = string(srcTb.recId);
stimIds = string(LMV.Param.stimIdList4);

fprintf('Analyzing chunking patterns for %d recordings and %d sentences...\n', ...
    numel(recIds), numel(stimIds));

% Initialize results storage
allResults = cell(numel(recIds), numel(stimIds));

% Process each recording
for i = 1 : numel(recIds)
    recId = recIds(i);
    fprintf('Processing recording %s (%d/%d)...\n', recId, i, numel(recIds));

    % Check if results already exist
    cachePath = fullfile(anaDir, sprintf('%s_chunking.mat', recId));
    if exist(cachePath, 'file')
        fprintf('Skipping recording %s (results already exist)...\n', recId);
        continue
    end
    
    % Load session data
    se = NP.SE.LoadSession(srcTb.path(i));
    ops = NP.Param.Enrich;
    [se, senTb] = LMV.SE.Transform(se, ["enrich", "morph", "sentence"], ops);

    % Keep only the long sentences
    senTb = senTb(ismember(senTb.stimId, stimIds), :);

    % Analyze each sentence
    for j = 1 : height(senTb)
        stimId = senTb.stimId(j);
        seSen = senTb.se(j);
        
        % Slice and resample speed data
        [tt, tv] = seSen.GetTable('taskTime', 'taskValue');
        isTemp = tv.trialNum == tv.tempTrialNum(1);
        tWin = [tt.prodMatchOn(isTemp), tt.prodMatchOff(isTemp)];
        tEdges = tWin(1) : 0.01 : tWin(2);
        morphTb = seSen.ResampleTimeSeries("morph", tEdges);

        % Get speed data from morph table for all epochs
        t = morphTb.time{1};
        S = cat(2, morphTb.speed{:});

        % Get the template sentence as an example
        if ~iscell(tt.prod)
            tt.prod = num2cell(tt.prod);
        end
        sen = tt.prod{isTemp};

        % Analyze word boundaries
        gapData = LMV.Chunk.AnalyzeWordBoundaries(seSen, morphTb, sen);
        summary = LMV.Chunk.CalculateSummaryStats(gapData);

        % Store results
        result = struct;
        result.recId = recId;
        result.stimId = stimId;
        result.sen = sen;
        result.time = t;
        result.speed = S;
        result.tWin = tWin;
        result.gapData = gapData;
        result.summary = summary;
        allResults{i, j} = result;
    end

    % Save individual results
    recResults = allResults(i, :);
    save(cachePath, 'recResults');

    f = MPlot.Figure(2291); clf
    LMV.Chunk.PlotSpeedMaps(recResults);
    MPlot.Paperize(f, 1.5, 1.7);
    exportgraphics(f, fullfile(anaDir, sprintf('%s_speed_maps.png', recId)));

    f = MPlot.Figure(2292); clf
    f.WindowState = "maximized";
    LMV.Chunk.PlotGapChanges(recResults);
    exportgraphics(f, fullfile(anaDir, sprintf('%s_duration_changes.png', recId)));
end

%% 

% Save combined results
combinedResults = struct();
combinedResults.recIds = recIds;
combinedResults.stimIds = stimIds;
combinedResults.results = allResults;
combinedResults.timestamp = datetime('now');

save(fullfile(anaDir, 'combined_chunking_results.mat'), 'combinedResults');
fprintf('Combined results saved to %s\n', fullfile(anaDir, 'combined_chunking_results.mat'));

% Generate summary plots
LMV.Chunk.PlotSummaryResults(allResults, recIds, stimIds);











return
%% 

recId = "NP46_B1";
anaDir = LMV.Data.GetAnalysisDir("chunking", recId);

%% 

% Load original se
srcTb = LMV.Data.FindSource(recId);
se = NP.SE.LoadSession(srcTb.path(1));

% Morph se
ops = NP.Param.Enrich;
se = LMV.SE.Transform(se, ["enrich", "morph"], ops);

%% 

% List the IDs of the long sentences
stimIdLong = LMV.Param.stimIdList4(2);

% Remove all epochs except the long sentence
isSen = se.GetTable('taskValue').stimId == stimIdLong;
seSen = se.Duplicate;
seSen.RemoveEpochs(~isSen);

% Slice and resample speed data
[tt, tv] = seSen.GetTable('taskTime', 'taskValue');
isTemp = tv.trialNum == tv.tempTrialNum(1);
tWin = [tt.prodMatchOn(isTemp), tt.prodMatchOff(isTemp)];
tEdges = tWin(1) : 0.01 : tWin(2);
morphTb = seSen.ResampleTimeSeries("morph", tEdges);

% Get speed data from morph table for all epochs
t = morphTb.time{1};
S = cat(2, morphTb.speed{:});

%% Plot speed map

f = MPlot.Figure(2291); clf
tl = tiledlayout("vertical");
tl.Padding = 'compact';
tl.TileSpacing = "none";

sen = tt.prod{isTemp};
ax = nexttile;
sen.PlotTiers();
axis(ax, 'off');
ax.XLim = tWin;
MPlot.Axes(ax);

ax = nexttile;
imagesc(t-t(1), 1:size(S,2), S');
colormap('jet');
c = colorbar;
ylabel(c, 'Speed (X)');
xlabel('Time (s)');
MPlot.Axes(gca);

MPlot.Paperize(f, 1.5, .6);
exportgraphics(f, fullfile(anaDir, recId+"_speed_map.png"));

%% Evaluate changes in gap duration

% Get TGEvent objects for the sentence
sen = tt.prod{isTemp};
words = Cut(sen);
phones = arrayfun(@(x) Cut(Trim(x)), words, 'Uni', false);
gapWins = cellfun(@(x,y) [x(end), y(1)], phones(1:end-1), phones(2:end), 'Uni', false);
gapWins = double(cat(1, gapWins{:}));

% For each word boundary window
nGaps = size(gapWins, 1);
nRep = size(S, 2);
gapData = cell(nGaps, 1);
for i = 1 : nGaps
    % Slice original timestamps
    gapTb = seSen.SliceTimeSeries(morphTb, gapWins(i,:));

    % Get the speed data for this window
    x = (1 : height(gapTb))';
    y = cellfun(@(x) x(end) - x(1), gapTb.time0);
    
    % Fit sigmoid using fit function
    [ft, gof] = fit(x, y, 'poly1');

    % Store the data
    s = struct;
    s.win = gapWins(i,:);
    s.words = words(i:i+1);
    s.repInd = x;
    s.dur = y;
    s.durHat = ft(x);
    s.r2 = gof.adjrsquare;
    gapData{i} = s;
end

%% Plot changes in speed at word boundaries

f = MPlot.Figure(2292); clf
tl = tiledlayout("flow");
tl.Padding = 'compact';

nGaps = numel(gapData);
for i = 1:nGaps
    ax = nexttile;
    data = gapData{i};
    
    plot(data.repInd, data.dur, 'k.', 'MarkerSize', 8);
    hold on;
    plot(data.repInd, data.durHat, 'r-', 'LineWidth', 1.5);
    
    ax.XTick = data.repInd;
    ax.XLim = [0 data.repInd(end)+1];
    ax.YLim(1) = 0;
    ax.YLim(2) = max([ax.YLim(2) 0.15]);
    title(sprintf('%s %s (R² = %.3f)', data.words(1).GetParentLabel, data.words(2).GetParentLabel, data.r2));
    xlabel('Repeats');
    ylabel('Duration (s)');
    MPlot.Axes(ax);
end

MPlot.Paperize(f, 1.5, 1);
exportgraphics(f, fullfile(anaDir, recId+"_gap_changes.png"));



