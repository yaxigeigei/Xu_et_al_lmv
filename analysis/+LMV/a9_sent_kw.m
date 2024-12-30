%% Single-unit classification of sentences at each task phase

% Analysis folder
anaDir = LMV.Data.GetAnalysisDir('sent_resp');
srcTb = LMV.Data.FindSource([]);

%%

nSen = [4 14];
phaseNames = ["atten", "stim", "delay", "init", "prod"];

for i = 1 : height(srcTb)
    % Load response data
    recId = srcTb.recId{i};
    fprintf("\n%s\n", recId);
    dataPath = fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'extracted_resp', recId+"_phase-resp.mat");
    load(dataPath, 'clusTb', 'respTb', 'dropoutTb', 'phaseTb', 'stimIdTb', 'stimIdInfoTb');
    
    % Mask out responses when units are absent
    respTb{:,:}(dropoutTb{:,:}) = NaN;
    
    % Find responsive unit-phase
    % rTest = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'signrank', 'computed_tests', recId+"_clusTb.mat"));
    % rTest.sigTb = LMV.Resp.GetSigTable(rTest.clusTb);
    rTest = LMV.Resp.LoadPhaseResponseTest(recId);
    
    for k = 1 : numel(nSen)
        % Check computed
        ns = nSen(k);
        cacheDir = fullfile(anaDir, "sen"+ns, "computed_tests");
        if ~exist(cacheDir, 'dir')
            mkdir(cacheDir);
        end
        cachePath = fullfile(cacheDir, sprintf("%s_clusTb.mat", recId));
        if exist(cachePath, 'file')
            fprintf("\nSkip computed models\n%s\n", cachePath);
            continue
        end
        
        % Statistical tests
        fprintf("\nTest across %i sentences during %i task phases for %i units.\n", nSen(k), numel(phaseNames), width(respTb));
        clusTb{:,phaseNames} = NP.Resp.StimSelectivityKW(respTb, phaseTb(:,phaseNames), stimIdTb(:,1:ns), 'Mask', rTest.sigTb{:,:});
        
        % Save computed
        save(cachePath, 'clusTb');
    end
end

%% Save concatenated table

for k = 1 : numel(nSen)
    ns = nSen(k);
    
    % Load computed tests
    clusTbs = cell(height(srcTb), 1);
    for i = 1 : height(srcTb)
        recId = srcTb.recId{i};
        fprintf("Load cache for %s\n", recId);
        cachePath = fullfile(anaDir, "sen"+ns, "computed_tests", sprintf("%s_clusTb.mat", recId));
        load(cachePath, 'clusTb');
        clusTbs{i} = clusTb;
    end
    
    % Save concatenated table
    clusTb = cat(1, clusTbs{:});
    savePath = fullfile(anaDir, "sen"+ns, "computed_test_clusTb.mat");
    save(savePath, 'clusTb');
end
