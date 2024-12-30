%% Fit time-independent TRF encoding models

% Job command
% submit_job -q pia-batch.q -c 12 -m 48 -o /userdata/dxu/project_np/code/babble/f4_fit_ti.txt -n f4_fit_ti -x /data_store2/MATLAB/R2022a/bin/matlab /userdata/dxu/project_np/code/babble/f4_fit_ti.m

if ~ispc && ~ismac
    addpath(genpath("/userdata/dxu/project_np/code/third_party_matlab"));
    addpath("/userdata/dxu/project_np/code/babble");
end

%% Load data

ceSearch = MBrowse.Dir2Table(fullfile(LMV.Data.GetAnalysisDir, 'coding', 'stRF', 'ce', '*_ce.mat'));
cePaths = fullfile(ceSearch.folder, ceSearch.name);
ceArray = NP.CE.LoadSession(cePaths);

%% Fit models

sets = "artic3"; %["phone", "strf", "artic3"];
phases = ["prod", "stim", "feedback"];

for s = 1 : numel(sets)
    % Make a subfolder for the feature set
    sn = sets(s);
    mdlDir = LMV.Data.GetAnalysisDir('coding', 'stRF', sn, 'mdls');
    
    for p = 1 : numel(phases)
        pn = phases(p);
        mn = sn+"_"+pn;
        for d = 1 : numel(ceArray)
            % Check computed
            ce = ceArray(d);
            mdlFile = sprintf("%s_clusTb_%s.mat", NP.SE.GetID(ce), mn);
            mdlPath = fullfile(mdlDir, mdlFile);
            if exist(mdlPath, 'file')
                fprintf("\nSkip computed models from\n%s\n", mdlPath);
                continue
            end
            
            % Select all units
%             uInd = [];
            
            % Select speech responsive units
            sTest = LMV.Resp.LoadPhaseResponseTest(NP.SE.GetID(ce));
            if pn == "stim"
                uInd = find(sTest.sigTb.stim > 0);
            else
                uInd = find(sTest.sigTb.prod > 0);
            end
            
            % Use fixed lambda
            switch sn
                case 'phone'
                    L = 0;
                case 'artic'
                    L = 1;
                case 'strf'
                    L = 50;
                otherwise
                    L = 0;
            end
            
            % Fit models
            tic
            clusTb = LMV.RF.FitSlidingTimeRF(ceArray(d), uInd, 'FeatSet', sn, 'Target', pn, 'Lambda', L);
            toc
            
            % Save result
            save(mdlPath, 'clusTb');
        end
    end
end
