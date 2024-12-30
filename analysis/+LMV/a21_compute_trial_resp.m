%% Extract and cache single-trial neuronal responses

anaDir = LMV.Data.GetAnalysisDir('data');
srcTb = LMV.Data.FindSource([]);

%% Compute M2 sentence PETHs and average features

% Get file paths
m2Dir = fullfile(anaDir, 'se_m2');
m2epPaths = fullfile(m2Dir, strrep(srcTb.name, '_se.mat', '_se_m2_ep.mat'));

tbDir = LMV.Data.GetAnalysisDir('data', 'resp_m2_detrend_ex3');
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

%% Prepare pseudopopulation for decoding with leave-one-out CV

tbPaths = fullfile(anaDir, 'resp_m2_detrend_ex3', srcTb.recId+".mat");
cacheDir = LMV.Data.GetAnalysisDir('data', 'resp_m2_detrend_ex3_sen4_loo');

% Prepare data
nTrain = 7;
nBoot = 1000;
for i = 1 : numel(tbPaths)
    % Load for each recording
    load(tbPaths(i), "sen14Tb", "clusTb");
    
    % Keep four most repeated sentences
    [~, I] = maxk(sen14Tb.numTrial, 4);
    sen4Tb = sen14Tb(I,:);
    
    for j = 1 : height(sen4Tb)
        % Randomly generate test trial indices
        n = sen4Tb.numTrial(j);
        indTest = randi(n, [1 nBoot]);
        sen4Tb.indTest{j} = indTest;
        
        % Use the rest of trials for training
        indTrain = arrayfun(@(x) setdiff(1:n, x), indTest, 'Uni', false);
        indTrain = cat(1, indTrain{:})';
        
        % Resample training trials to target number
        I = repmat(1:size(indTrain,1), [1 nTrain]);
        I = I(1:nTrain);
        sen4Tb.indTrain{j} = indTrain(I,:);
    end
    
    % Convert repeat indices to trial indices
    indBase = cumsum([0; sen4Tb.numTrial(1:end-1)]);
    indTest = cellfun(@(x,y) x+y, sen4Tb.indTest, num2cell(indBase), 'Uni', false);
    indTrain = cellfun(@(x,y) x+y, sen4Tb.indTrain, num2cell(indBase), 'Uni', false);
    
    % Concatenate trials across sentences
    indTest = cat(1, indTest{:});
    indTrain = cat(1, indTrain{:});
    R = permute(cat(4, sen4Tb.resp{:}), [4 2 1 3]); % trial-by-unit-by-time
    
    % Only keep unique indexing
    [~, I] = unique([indTest; indTrain]', 'rows', 'stable');
    indTest = indTest(:,I);
    indTrain = indTrain(:,I);
    
    % Cache data
    s = struct;
    s.clusTb = clusTb;
    s.sen4Tb = sen4Tb(:,1:5);
    s.resp = R;
    s.respDimNames = ["trial", "unit", "time"];
    s.indTrain = indTrain;
    s.indTest = indTest;
    save(fullfile(cacheDir, srcTb.recId(i)+".mat"), '-struct', "s");
end
