%% RSA between sentence identity and single-unit responses in each task phase

srcTb = LMV.Data.FindSource([]);
mdlName = LMV.RSA.senMdl;

%% Load data

% Load unit responsiveness ressults
rTest = LMV.Resp.LoadPhaseResponseTest();

%% 

for i = 1 : height(srcTb)
    % Check computed
    cacheDir = LMV.Data.GetAnalysisDir('coding', mdlName, "computed_mdls");
    recId = srcTb.recId{i};
    cachePath = fullfile(cacheDir, sprintf("%s_clusTb.mat", recId));
    if exist(cachePath, 'file')
        fprintf("\nSkip computed models\n%s\n", cachePath);
        continue
    end
    
    % Load response data
    fprintf("\n%s\n", recId);
    dataPath = fullfile(LMV.Data.GetAnalysisDir, 'phase_resp', 'extracted_resp', recId+"_phase-resp.mat");
    load(dataPath, 'clusTb', 'respTb', 'dropoutTb', 'phaseTb', 'stimIdTb');
    
    % Mask out responses when units are absent
    respTb{:,:}(dropoutTb{:,:}) = NaN;
    
    % Compute for each task phase
    phaseNames = ["atten", "stim", "delay", "init", "prod"];
    stimIds = string(stimIdTb.Properties.VariableNames);
    for p = 1 : numel(phaseNames)
        pn = phaseNames(p);
        
        % Select responsive units
        isRec = rTest.clusTb.recId == recId;
        uMask = rTest.sigTb.(pn)(isRec) > 0;
        
        % % Select all units
        % uMask = true(width(rTb), 1);
        
        % Classification
        fprintf("\nClassify %i sentences during %s by each of %i units.\n", numel(stimIds), pn, sum(uMask));
        clusTb.(pn)(uMask) = LMV.RSA.Sentences(respTb(:,uMask), [phaseTb stimIdTb], stimIds, pn);
    end
    
    % Save computed
    save(cachePath, 'clusTb');
end
