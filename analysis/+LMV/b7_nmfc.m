%% NMF clustering on grand average PETHs

% Job command
% submit_job -q pia-batch.q -c 12 -m 48 -o /userdata/dxu/project_np/code/babble/np_a2_nmfc.txt -n np_a2_nmfc -x /data_store2/MATLAB/R2022a/bin/matlab /userdata/dxu/project_np/code/babble/np_a2_nmfc.m

if ~ispc && ~ismac
    addpath(genpath("/userdata/dxu/project_np/code/third_party_matlab"));
    addpath("/userdata/dxu/project_np/code/babble");
end

anaDir = LMV.Data.GetAnalysisDir('phase_resp', 'nmf');

%% Load data

% Load grand average (M2) PETHs
load(fullfile(LMV.Data.GetAnalysisDir, 'data', 'ce_m2_session-avg.mat'), 'ce');

% Remove non-responsive units from ce
sTest = LMV.Resp.LoadPhaseResponseTest();
[~, I] = MMath.SortLike(sTest.clusTb.clusId, ce.clusTb.clusId);
isResp = any(sTest.sigTb{I,:}, 2);
ce.RemoveUnits(~isResp);

%% NMF clustering with bootstrap evaluation

regions = ["mPrCG", "vPrCG", "IFG", "STG"];
kList = 1 : 15; % the number of clusters to test

for r = 1 : numel(regions)
    regionCachePath = fullfile(anaDir, "computed_nmfc_"+regions(r)+".mat");
    if exist(regionCachePath, 'file')
        fprintf("\nCache for '%s' already exists.\n", regions(r))
        continue
    end
    
    % Clustering
    s = LMV.NMF.BatchClustering(ce, kList, regions(r));
    
    % Cache results
    save(regionCachePath, '-struct', 's');
end

return
%% Save data for python

% s = struct;
% [~, s.peth] = ce.GetArray('resp', 'Normalization', 'minmaxsoft');
% dataPath = fullfile(anaDir, 'computed_peth_nmf_all_py.mat');
% save(dataPath, '-struct', 's');
