%% Extract and cache single-trial neuronal responses

anaDir = LMV.Data.GetAnalysisDir('data');
srcTb = LMV.Data.FindSource([]);

%% Compute M2 sentence PETHs and average features

% Get file paths
m2Dir = fullfile(anaDir, 'se_m2');
m2epPaths = fullfile(m2Dir, strrep(srcTb.name, '_se.mat', '_se_m2_ep.mat'));

tbDir = LMV.Data.GetAnalysisDir('data', 'resp_m2_ex3');
tbPaths = fullfile(tbDir, strrep(srcTb.name, '_se.mat', '.mat'));

for i = 1 : height(srcTb)
    % Check computed
    if isfile(tbPaths{i})
        fprintf("\nThe responses has been computed and saved:\n'%s'\n", tbPaths{i});
        continue
    end
    
    % Load epoch-resliced se
    se = NP.SE.LoadSession(m2epPaths{i});
    
    % De-trend spike rates
    NP.Unit.DetrendSpikeRates(se);
    
    % Remove the first 3 trials
    se.RemoveEpochs(1:3);
    
    % Split and group trials of the same sentences
    % (No need to use LMV.SE.SplitBySentence since it's already resliced and cleaned)
    senTb = NP.TaskBaseClass.SplitBySentence(se);
    
    % Initialize table with the 14 sentences
    sen14Tb = table;
    sen14Tb.stimId = LMV.Param.stimIdList14;
    sen14Tb.stimText = LMV.Param.stimTextList14;
    sen14Tb.recId(:) = senTb.recId(1);
    sen14Tb.subjectId(:) = senTb.subjectId(1);
    
    for s = 1 : height(sen14Tb)
        % Check if the sentence is present
        isSen = find(sen14Tb.stimId(s)==senTb.stimId);
        if ~any(isSen)
            fprintf("\n%s - not found\n", sen14Tb.stimText(s));
            continue
        else
            fprintf("\n%s\n", sen14Tb.stimText(s));
        end
        
        % Resample responses
        seSen = senTb.se(isSen);
        tt = seSen.GetTable('taskTime');
        
        ops = NP.Param.Resample;
        ops.rsWin = [tt.trialOn(1) tt.matchOff(1)] + [-1 1]*0.5;
        ops.rsBinSize = 0.01;
        ops.rsArgs = {};
        
        rTb = NP.SE.ResampleResponses(seSen, ops);
        
        % Reshape responses to array
        R = permute(rTb{:,2:end}, [3 2 4 1]);
        R = cell2mat(R);
        
        sen14Tb.numTrial(s) = senTb.numTrial(isSen);
        sen14Tb.resp{s} = R; % time-by-unit-by-sentence-by-trial
    end
    
    % Save to cache
    clusTb = NP.Unit.GetClusTb(se);
    save(tbPaths{i}, "clusTb", "sen14Tb");
end

%% Create pseudopopulation for SCA

% Load for each recording
tbPaths = fullfile(anaDir, 'resp_m2_ex3', srcTb.recId+".mat");
sResp = arrayfun(@load, tbPaths);

% Concatenate data across sentences and recordings
nTargetRep = 8;
rSenCat = cell(size(sResp));
for i = 1 : numel(sResp)
    % Evenly resample trials to the target number of repeats
    sen14Tb = sResp(i).sen14Tb;
    for j = 1 : height(sen14Tb)
        n = sen14Tb.numTrial(j);
        if ~n
            sen14Tb.resp{j} = [];
            continue
        end
        ind = repmat(1:n, [1 nTargetRep]);
        ind = ind(1:nTargetRep);
        R = sen14Tb.resp{j}(:,:,:,ind);
        sen14Tb.resp{j} = R;
    end
    
    % Set NaN array to missing sentences
    isEpt = cellfun(@isempty, sen14Tb.resp);
    sen14Tb.resp(isEpt) = {NaN(size(R))};
    
    % Concatenate sentence dimension
    rSenCat{i} = cat(3, sen14Tb.resp{:});
end

% Cache combined data
sCat = struct;
sCat.clusTb = cat(1, sResp.clusTb);
sCat.resp = cat(2, rSenCat{:});
sCat.dimNames = ["time", "unit", "sentence", "repeat"];
save(fullfile(anaDir, "resp_m2_ex3_rep8_cat.mat"), '-struct', "sCat");
