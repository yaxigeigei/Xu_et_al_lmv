%% Extract single-trial spike rate for each unit in each task phase

anaDir = LMV.Data.GetAnalysisDir('phase_resp', 'extracted_resp');
srcTb = LMV.Data.FindSource([]);

%% 

for i = 1 : height(srcTb)
    % Check cache
    cachePath = fullfile(anaDir, srcTb.recId{i}+"_phase-resp.mat");
    if exist(cachePath, 'file')
        fprintf("\nSkip extracted responses for %s.\n", srcTb.recId{i});
        continue
    end
    
    % Load se
    se = NP.SE.LoadSession(srcTb.path{i});
    se = LMV.SE.Transform(se);
    
    % Extract responses
    s = LMV.Resp.ExtractPhaseResponses(se);
    
    % Save cache
    clusTb = NP.Unit.GetClusTb(se);
    clusTb = NP.Unit.AddRecMeta(se, clusTb);
    s.clusTb = clusTb;
    save(cachePath, '-struct', 's');
end
