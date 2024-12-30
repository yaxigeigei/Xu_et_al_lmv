%% Compute sentence PETHs from M2 unit cache

mdlName = LMV.Linker.currentModel;
anaDir = LMV.Data.GetAnalysisDir("linker", mdlName);

%% 

cachePath = fullfile(anaDir, "linker_clusTb_peth.mat");

if ~isfile(cachePath)
    % Find linker units
    clusTb = LMV.Linker.LoadClusTb(mdlName);
    
    % Load cached M2 data
    clusTb.uCache = NP.Unit.LoadUnitCache(clusTb.clusId, 'DataSource', 'm2');
    
    % Compute sentence responses
    clusTb = LMV.Linker.ComputeSentenceResponseFromCache(clusTb);
    
    save(cachePath, "clusTb");
end
