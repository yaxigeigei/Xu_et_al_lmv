%% Compute PETHs

srcTb = LMV.Data.FindSource([]);

%% 

seTbDir = LMV.Data.GetAnalysisDir("linker", 'computed_seTb');

ceDir1 = LMV.Data.GetAnalysisDir("linker", "computed_ce");

verName = "sem-discounted";
ceDir2 = LMV.Data.GetAnalysisDir("linker", "computed_ce_"+verName);

for k = 1 : height(srcTb)
    % Check computed
    cePath = fullfile(ceDir1, strrep(srcTb.name{k}, '_se.mat', '_ce.mat'));
    cePath2 = fullfile(ceDir2, strrep(srcTb.name{k}, '_se.mat', '_ce.mat'));
    if isfile(cePath) && isfile(cePath2)
        continue
    end
    
    % Load seTb
    seTbPath = fullfile(seTbDir, strrep(srcTb.name{k}, '_se.mat', '_seTb.mat'));
    load(seTbPath, "seTb");
    
    % Compute PETHs
    if isfile(cePath)
        fprintf("\nSkip computed ce\n%s\n", cePath);
    else
        ce = LMV.SE.ComputeSentencePETH(seTb);
        ce.clusTb.mi = []; % remove variable sized column to enable later concatnation
        save(cePath, 'ce', '-v7.3');
    end
    
    % Compute error-discounted PETHs
    if isfile(cePath2)
        fprintf("\nSkip computed ce\n%s\n", cePath2);
    else
        ce = LMV.SE.ComputeSentencePETH(seTb, 'Modification', verName);
        ce.clusTb.mi = []; % remove variable sized column to enable later concatnation
        save(cePath2, 'ce', '-v7.3');
    end
end
