%% Compute and cache unit lifetime sparseness measured by kurtosis of firing rate distribution

cacheDir = LMV.Data.GetAnalysisDir("units", "computed_sparsity");
srcTb = LMV.Data.FindSource([]);

%% Load ce

cePaths = fullfile(LMV.Data.GetAnalysisDir, "data", "ce_m1_sentence-avg", srcTb.recId+"_ce_m1.mat");
ceArray = NP.CE.LoadSession(cePaths);

%% Calculate sparsity for each recording and cache to files

for i = 1 : numel(ceArray)
    ce = ceArray(i);
    recId = string(NP.SE.GetID(ce));
    fprintf('Processing %s (%d/%d)\n', recId, i, numel(ceArray));
    
    % Define output path
    outPath = fullfile(cacheDir, recId + "_clusTb.mat");
    
    % Skip if already computed
    if exist(outPath, 'file')
        fprintf('  Already computed, skipping\n');
        continue
    end
    
    % Get unit table
    clusTb = NP.Unit.GetClusTb(ce);
    nUnits = height(clusTb);

    % Create task time windows
    tt = ce.GetTable('taskTime');
    tWins = { ...
        [tt.stimOn, tt.stimOff], ...
        [tt.stimOff, tt.cue3On], ...
        [tt.prodMatchOn, tt.prodMatchOff], ...
        [tt.stimOn, tt.prodMatchOff] ...
    };
    tWins = cellfun(@(x) mat2cell(x, ones(size(x,1),1), 2), tWins, 'Uni', false); % convert to a cell array of 2-element windows
    tWinNames = ["stim", "delay", "prod", "trial"];

    % Slice spike rates
    srTbs = cellfun(@(x) ce.SliceTimeSeries('spikeRate', x), tWins, 'Uni', false);

    % Initialize 
    spTbs = cell(nUnits, 1);

    % Loop through units to compute sparsity
    for u = 1 : nUnits
        ss = cell(numel(tWins), 1);
        for p = 1 : numel(tWins)
            % Concatenate spike rates across all slices
            r = cell2mat(srTbs{p}.(u+1));
            r = r(~isnan(r)); % shouldn't be any but just in case
            
            % Define bin edges and centers
            binEdges = linspace(0, max(r), 51)';
            binCenters = binEdges(1:end-1) + diff(binEdges)/2;
            
            % Compute histogram of spike rates
            frHist = histcounts(r, binEdges, "Normalization", "probability");
            
            % Calculate kurtosis of spike rate distribution
            kurtValue = kurtosis(r);
            
            ss{p}.winName = tWinNames(p);
            ss{p}.binCenters = binCenters;
            ss{p}.binHeights = frHist;
            ss{p}.kurtosis = kurtValue;
        end
        spTb = struct2table([ss{:}], 'RowNames', tWinNames);
        spTbs{u} = spTb;
    end
    
    % Add results to clusTb
    clusTb.sparsityTb = spTbs;
    
    % Save results
    save(outPath, 'clusTb');
    fprintf('  Saved to %s\n', outPath);
end

%% Optional: Load and combine results from all recordings

% Load all cached tables
clusTbList = cell(height(srcTb), 1);
for i = 1 : height(srcTb)
    recId = srcTb.recId(i);
    tbPath = fullfile(cacheDir, recId + "_clusTb.mat");
    if exist(tbPath, 'file')
        data = load(tbPath);
        clusTbList{i} = data.clusTb;
    end
end

% Combine all tables
clusTb = vertcat(clusTbList{:});

% Optional: Save combined table
combinedPath = fullfile(fileparts(cacheDir), "computed_sparsity_clusTb.mat");
save(combinedPath, 'clusTb');
fprintf('Combined results saved to %s\n', combinedPath);
