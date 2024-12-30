%% 

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'phasic');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end
srcTb = NP.Data.FindSource('lmv');

%% Compute sentence (M1) PETHs and average features

% Get file paths
m1Dir = fullfile(NP.Data.GetAnalysisRoot, 'data', 'se_m1');
ceDir = fullfile(anaDir, 'ce_m1_sentence-avg');
if ~exist(ceDir, 'dir')
    mkdir(ceDir);
end
m1epPaths = fullfile(m1Dir, strrep(srcTb.name, '_se.mat', '_se_m1_ep.mat'));
cePaths = fullfile(ceDir, strrep(srcTb.name, '_se.mat', '_ce_m1.mat'));

for i = 1 : height(srcTb)
    % Check computed
    if exist(cePaths{i}, 'file')
        fprintf("\nThe M1 ce has been computed and saved at '%s'\n", cePaths{i});
        continue
    end
    
    % Load epoch-resliced M1 se
    se = NP.SE.LoadSession(m1epPaths{i});
    
    % Remove spikeRate table to make NP.SE.ResampleResponses recompute spike rates by binning spike times
    se.RemoveTable('spikeRate');
    
    % Split and group trials of the same sentences
    % (No need to use LMV.SE.SplitBySentence since it's already resliced and cleaned)
    senTb = NP.TaskBaseClass.SplitBySentence(se);
    
    % Compute PETHs
    [ce, senTb] = LMV.SE.ComputeSentencePETH(senTb);
    
    save(cePaths{i}, 'ce');
end

return
%% Compute sentence (M1) PETHs and average features

% Get file paths
m1Dir = fullfile(NP.Data.GetAnalysisRoot, 'data', 'se_m1');
ceDir = fullfile(anaDir, 'ce_m1_sentence-avg');
if ~exist(ceDir, 'dir')
    mkdir(ceDir);
end
m1epPaths = fullfile(m1Dir, strrep(srcTb.name, '_se.mat', '_se_m1_ep.mat'));
cePaths = fullfile(ceDir, strrep(srcTb.name, '_se.mat', '_ce_m1.mat'));

for i = 1 : height(srcTb)
    % Check computed
    if exist(cePaths{i}, 'file')
        fprintf("\nThe M1 ce has been computed and saved at '%s'\n", cePaths{i});
        continue
    end
    
    % Load epoch-resliced M1 se
    se = NP.SE.LoadSession(m1epPaths{i});
    
    % Remove spikeRate table to make NP.SE.ResampleResponses recompute spike rates by binning spike times
    se.RemoveTable('spikeRate');
    
    % Split and group trials of the same sentences
    % (No need to use LMV.SE.SplitBySentence since it's already resliced and cleaned)
    senTb = NP.TaskBaseClass.SplitBySentence(se);
    
    for s = 1 : height(senTb)
        % Extract and resample features and responses
        fprintf("\n%s\n", senTb.stimText(s));
        seSen = senTb.se(s);
        tt = seSen.GetTable('taskTime');
        
        ops = NP.Param.Resample;
        ops.rsWin = [tt.trialOn(1) tt.matchOff(1)] + [-1 1]*0.5;
        ops.rsBinSize = 0.025;
        ops.rsArgs = {};
        rTb = NP.SE.ResampleResponses(seSen, ops);
        
        % Compute trial average
        mrTb = NP.SE.MeanTimeseries(rTb);
        mrTb2 = NP.SE.MeanTimeseries(rTb, 'StatFunc', @MMath.MedianStats); % used for identifying phasic responses
        
        % Construct ce
        ce = NP.CodingExplorer();
        ce.SetTable('resp', mrTb(1,:), 'timeSeries');
        ce.SetTable('respMed', mrTb2(1,:), 'timeSeries');
        LMV.SE.AddTemplateTrial(ce, seSen);
        ce.SetColumn('taskValue', 'numTrial', senTb.numTrial(s));
        senTb.ce(s) = ce;
    end
    
    % Concatenate ce of all sentences
    ce = Merge(senTb.ce);
    ce.userData = se.userData;
    ce.userData.rsOps = ops;
    
    save(cePaths{i}, 'ce');
end
