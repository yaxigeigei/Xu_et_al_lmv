%% 

% Job command
% submit_job -q pia-batch.q -c 12 -m 48 -o /userdata/dxu/project_np/code/babble/f4_fit_speech.txt -n f4_fit_speech -x /data_store2/MATLAB/R2022a/bin/matlab /userdata/dxu/project_np/code/babble/f4_fit_speech.m

if ~ispc && ~ismac
    addpath(genpath("/userdata/dxu/project_np/code/third_party_matlab"));
    addpath("/userdata/dxu/project_np/code/babble");
end

%% Load data

% Load ce
ceSearch = MBrowse.Dir2Table(fullfile(NP.Data.GetAnalysisRoot, 'trf', 'ce', '*_ce.mat'));
cePaths = fullfile(ceSearch.folder, ceSearch.name);
ceArray = NP.CE.LoadSession(cePaths);

% Load unit responsiveness table
sTest = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'ttest', "computed_ttest_clusTb.mat"));
sigTb = NP.Resp.GetSigTable(sTest.clusTb, 'VariableNames', ["stim" "prod"]);

%% Fit models for hyperparam search

sets = ["combo3pm", "combo3akt"];
phases = ["stim", "prod", "feedback"];

clusTbs = cell(numel(ceArray), numel(sets), numel(phases));

for s = 1 : numel(sets)
    % Make a subfolder for the feature set
    sn = sets(s);
    mdlDir = fullfile(NP.Data.GetAnalysisRoot, 'trf', sn, 'mdls');
    if ~exist(mdlDir, 'dir')
        mkdir(mdlDir);
    end
    
    for p = 1 : numel(phases)
        pn = phases(p);
        mn = sn+"_"+pn;
        for d = 1 : numel(ceArray)
            % Check computed
            ce = ceArray(d);
            mdlFile = sprintf("%s_clusTb_%s_HS.mat", NP.SE.GetID(ce), mn);
            mdlPath = fullfile(mdlDir, mdlFile);
            if exist(mdlPath, 'file')
                fprintf("\nLoad computed hyperparam search from\n%s\n", mdlPath);
                load(mdlPath, 'clusTb');
                clusTbs{d,s,p} = clusTb;
                continue
            end
            
            % Select all units
%             uInd = [];
            
            % Select speech responsive units
            isRec = sTest.clusTb.recId == NP.SE.GetID(ce);
            if pn == "stim"
                uInd = find(sigTb{isRec,"stim"} > 0);
            else
                uInd = find(sigTb{isRec,"prod"} > 0);
            end
            
            % Fit model
            tic
            clusTb = NP.TRF.FitSpeechEncoding(ceArray(d), uInd, 'FeatSet', sn, 'Target', pn);
            toc
            
            % Save result
            clusTbs{d,s,p} = clusTb;
            save(mdlPath, 'clusTb');
        end
    end
end

%% Fit models with optimal hyperparam and bootstrap null

for s = 1 : numel(sets)
    sn = sets(s);
    mdlDir = fullfile(NP.Data.GetAnalysisRoot, 'trf', sn, 'mdls');
    
    for p = 1 : numel(phases)
        pn = phases(p);
        mn = sn+"_"+pn;
        
        % Find mean lambda
        clusTbCat = cat(1, clusTbs{:,s,p});
        L = NP.Fit.GetOptimalLambda(clusTbCat.(mn+"_HS"), 'ConvexOnly', true);
        mL = 10.^median(log10(L), 'omitnan');
        fprintf("\nMedian lambda of %s models is %.2f (10^%.2f).\n", mn, mL, log10(mL));
        
        for d = 1 : numel(ceArray)
            % Check computed
            ce = ceArray(d);
            mdlFile = sprintf("%s_clusTb_%s.mat", NP.SE.GetID(ce), mn);
            mdlPath = fullfile(mdlDir, mdlFile);
            if exist(mdlPath, 'file')
                fprintf("\nSkip computed models from\n%s\n", mdlPath);
                continue
            end
            
            % Select units
            uInd = find(~cellfun(@isempty, clusTbs{d,s,p}.(mn+"_HS"))); % those used in hyperparam search
            % uInd = []; % all units
            
            % Fit model
            tic
            clusTb = NP.TRF.FitSpeechEncoding(ce, uInd, 'FeatSet', sn, 'Target', pn, 'Lambda', mL);
            toc
            
            % Save result
            save(mdlPath, 'clusTb');
        end
    end
end
