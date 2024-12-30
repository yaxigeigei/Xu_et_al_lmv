%% Single-unit classification of sentences at each task phase

anaDir = LMV.Data.GetAnalysisDir('coding', LMV.UnitDecode.sen4Mdl);
srcTb = LMV.Data.FindSource([]);

%% 

% Load unit sentence selectivity results
sTest = LMV.Resp.LoadSentenceSelectTest();

% Load unit responsiveness ressults
rTest = LMV.Resp.LoadPhaseResponseTest();

%%

for i = 1 : height(srcTb)
    % Load response data
    recId = srcTb.recId{i};
    fprintf("\n%s\n", recId);
    dataPath = fullfile(LMV.Data.GetAnalysisDir, 'phase_resp', 'extracted_resp', recId+"_phase-resp.mat");
    load(dataPath, 'clusTb', 'respTb', 'dropoutTb', 'phaseTb', 'stimIdTb');
    
    % Mask out responses when units are absent
    respTb{:,:}(dropoutTb{:,:}) = NaN;
    
    % 
    nSen = 4;
    phaseNames = ["atten", "stim", "delay", "init", "prod"];
    for k = 1 : numel(nSen)
        % Find the most repeated sentences
        ns = nSen(k);
        stimIds = string(stimIdTb.Properties.VariableNames(1:ns));
        
        % Check computed
        cacheDir = fullfile(anaDir, "computed_mdls");
        if ~exist(cacheDir, 'dir')
            mkdir(cacheDir);
        end
        cachePath = fullfile(cacheDir, sprintf("%s_clusTb.mat", recId));
        if exist(cachePath, 'file')
            fprintf("\nSkip computed models\n%s\n", cachePath);
            continue
        end
        
        for p = 1 : numel(phaseNames)
            pn = phaseNames(p);
            
            % Select responsive units
            isRec = rTest.clusTb.recId == recId;
            uMask = rTest.sigTb.(pn)(isRec) > 0;
            
%             % Select all units
%             uMask = true(width(rTb), 1);
            
            % Classification
            fprintf("\nClassify %i sentences during %s by each of %i units.\n", nSen(k), pn, sum(uMask));
            [clusTb.(pn)(uMask), clusTb.(pn+"Null")(uMask)] = LMV.UnitDecode.Sentences(respTb(:,uMask), [phaseTb stimIdTb], stimIds, pn);
        end
        
        % Save computed
        save(cachePath, 'clusTb');
    end
end
