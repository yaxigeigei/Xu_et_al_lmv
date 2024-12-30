%% 

% Job command
% submit_job -q pia-batch.q -c 12 -m 48 -o /userdata/dxu/project_np/code/babble/np_a4_cla_phase.txt -n np_a4_cla_phase -x /data_store2/MATLAB/R2022a/bin/matlab /userdata/dxu/project_np/code/babble/np_a4_cla_phase.m

if ~ispc && ~ismac
    addpath(genpath("/userdata/dxu/pkg/matlab"));
    addpath("/userdata/dxu/project_np/code/babble");
end

%% Load ce

ceSearch = MBrowse.Dir2Table(fullfile(LMV.Data.GetAnalysisDir, 'coding', 'stRF', 'ce', '*_ce.mat'));
cePaths = fullfile(ceSearch.folder, ceSearch.name);
ceArray = NP.CE.LoadSession(cePaths);

%% Fit trial classification models

anaDir = LMV.Data.GetAnalysisDir('coding', 'phase');
mdlDir = fullfile(anaDir, 'computed_phase_mdls');
if ~exist(mdlDir, 'dir')
    mkdir(mdlDir);
end

mdlNames = {'phase'};

for k = 1 : numel(mdlNames)
    mn = mdlNames{k};
    for i = 1 : numel(ceArray)
        % Check computed
        mdlFile = strrep(ceSearch.name{i}, '_ce.mat', ['_mdlTb_' mn '.mat']);
        mdlPath = fullfile(mdlDir, mdlFile);
        if exist(mdlPath, 'file')
            fprintf("\nSkip computed models\n%s\n", mdlPath);
            continue
        end
        
        % Fit models
        mdlTb = LMV.Decode.ClassifyPhases(ceArray(i), mn);
        
        % Save result
        fprintf("\nSave computed models\n%s\n", mdlPath);
        save(mdlPath, 'mdlTb');
    end
end

return
%% Fit phase classification models

mdlNames = {'phase'};

for k = 1 %: numel(mdlNames)
    mn = mdlNames{k};
    for i = 1 %: numel(ceArray)
        mdlTb = LMV.Decode.ClassifyPhases(ceArray(i), mn);
    end
end

%% 

mdl = mdlTb.phase{1};

f = MPlot.Figure(1); clf
cm = confusionchart(f, mdl.Y, mdl.kfoldPredict);
cm.Normalization = 'row-normalized';
title(sprintf("Accuracy %.2f", 1-mdl.kfoldLoss));


s = mdlTb.phaseNull{1};

for i = 1 : 5
    f = MPlot.Figure(1+i); clf
    cm = confusionchart(f, s.confMat(:,:,i));
    cm.Normalization = 'row-normalized';
    title(sprintf("Accuracy %.2f", 1-s.kfoldLoss(i)));
end





