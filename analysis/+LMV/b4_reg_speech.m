%% 

% Job command
% submit_job -q pia-batch.q -c 12 -m 48 -o /userdata/dxu/project_np/code/babble/np_a4_reg_speech.txt -n np_a4_reg_speech -x /data_store2/MATLAB/R2022a/bin/matlab /userdata/dxu/project_np/code/babble/np_a4_reg_speech.m

if ~ispc && ~ismac
    addpath(genpath("/userdata/dxu/pkg/matlab"));
    addpath("/userdata/dxu/project_np/code/babble");
end

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'decode');
mdlDir = fullfile(anaDir, 'computed_speech_mdls');
if ~exist(mdlDir, 'dir')
    mkdir(mdlDir);
end

%% Load ce

ceSearch = MBrowse.Dir2Table(fullfile(LMV.Data.GetAnalysisDir, 'coding', 'stRF', 'ce', '*_ce.mat'));
cePaths = fullfile(ceSearch.folder, ceSearch.name);
ceArray = NP.CE.LoadSession(cePaths);

%% Search speech model hyperparameter space

mdlNames = {'stim', 'prod'};

mdlTbs = cell(numel(ceArray), numel(mdlNames));

for k = 1 : numel(mdlNames)
    mn = mdlNames{k};
    for i = 1 : numel(ceArray)
        % Check computed
        mdlFile = strrep(ceSearch.name{i}, '_ce.mat', ['_mdlTb_' mn 'HS.mat']);
        mdlPath = fullfile(mdlDir, mdlFile);
        if exist(mdlPath, 'file')
            fprintf("\nLoad computed hyperparam search from\n%s\n", mdlPath);
            load(mdlPath, 'mdlTb');
            mdlTbs{i,k} = mdlTb;
            continue
        end
        
        % Fit models
        tic
        mdlTb = LMV.Decode.FitSpeechDecoding(ceArray(i), mn);
        toc
        
        % Save result
        mdlTbs{i,k} = mdlTb;
        save(mdlPath, 'mdlTb');
    end
end

%% Fit models

% Find mean lambda
mdlTbCat = cat(1, mdlTbs{:});
L = NP.Fit.GetOptimalLambda(mdlTbCat{:,2:end}, 'ConvexOnly', true);
mL = 10.^median(log10(L), 'all', 'omitnan');
fprintf("\nMedian lambda of speech decoding models is %.2f (10^%.2f).\n", mL, log10(mL));

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
        tic
        mdlTb = LMV.Decode.FitSpeechDecoding(ceArray(i), mn, 'Lambda', mL);
        toc
        
        % Save result
        save(mdlPath, 'mdlTb');
    end
end
