%% Determine unit responsiveness using t-test

anaDir = LMV.Data.GetAnalysisDir('phase_resp', 'ttest');
srcTb = LMV.Data.FindSource([]);

%% 

% Specify the task phases of interest
phaseNames = ["atten", "stim", "delay", "init", "prod"];

% Make folder to save cache files
cacheDir = fullfile(anaDir, "computed_tests");
if ~exist(cacheDir, 'dir')
    mkdir(cacheDir);
end

% Compute responsiveness for each recording
for i = 1 : height(srcTb)
    % Check cache status
    recId = srcTb.recId{i};
    cachePath = fullfile(cacheDir, recId+"_clusTb.mat");
    if exist(cachePath, 'file')
        fprintf("\nTest result of %s has been computed and cached\n", recId);
        continue
    end
    
    % Load response data
    load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'extracted_resp', recId+"_phase-resp.mat"), 'clusTb', 'respTb', 'dropoutTb', 'phaseTb', 'stimIdTb');
    respTb{:,:}(dropoutTb{:,:}) = NaN; % set dropout responses to NaN
    
    % Compute responsiveness
    fprintf("\nCompute task phase response for %s\n", recId);
    pTb = NP.Resp.PhaseTTest(respTb, [phaseTb stimIdTb], phaseNames, LMV.Param.stimIdList4);
    
    % Save result
    clusTb = [clusTb, pTb];
    save(cachePath, 'clusTb', 'phaseNames');
end

%% Make group labels (for activation only)

% Load and concatenate test tables
nRec = height(srcTb);
clusTbs = cell(nRec, 1);
for i = 1 : nRec
    load(fullfile(cacheDir, srcTb.recId{i}+"_clusTb.mat"));
    clusTbs{i} = clusTb;
end
clusTb = cat(1, clusTbs{:});

% Find the biggest response across significant tests for each task phase
pTb = clusTb(:, phaseNames);
sTb = LMV.Resp.GetSigTable(clusTb, 'Side', 'right', 'Minimum', false);
rTb = clusTb(:, phaseNames+"Resp");
for i = 1 : width(rTb)
    R = rTb.(i);
    R(~sTb.(i)) = 0;
    [~, k] = max(abs(R), [], 2);
    k = sub2ind(size(R), (1:size(R,1))', k);
    rTb.(i) = R(k);
end

% Construct task phase indices
R = rTb{:,:};
I = cumsum(ones(size(R)), 2);
I(R==0) = 6;
% I(R<0) = I(R<0) + numel(phaseNames); % indexing inhibition types after activation
I = categorical(I, 1:6, ["atten", "stim", "delay", "init", "prod", "none"]);

% Find best task phase activation
for i = 1 : numel(phaseNames)
    % Find the indices of maximum response
    [~, k] = max(abs(R), [], 2);
    k = sub2ind(size(R), (1:size(R,1))', k);
    tR = R(k);
    tId = I(k);
    
    % Inherit previous tier in the absence of further tier
    if i > 1
        m = tId=="none";
        tId(m) = clusTb.("tId"+(i-1))(m);
        tR(m) = clusTb.("tR"+(i-1))(m);
    end
    
    % Add results to table
    clusTb.("tR"+i) = tR;
    clusTb.("tId"+i) = tId;
    
    % Exclude max resp for next round
    R(k) = 0;
    I(k) = "none";
end

%% Save concatenated table

save(fullfile(anaDir, 'computed_test_clusTb.mat'), 'clusTb', 'phaseNames');
