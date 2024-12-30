%% Read template script

% Command
% /data_store2/MATLAB/R2022a/bin/matlab -nodesktop -nosplash -r "run('/userdata/dxu/project_np/analysis_lmv/coding/pop_sen4_ecoc/a21_resp_cla_jobs.m'); exit"

if ~ispc && ~ismac
    addpath(genpath("/userdata/dxu/project_np/code/third_party_matlab"));
    addpath("/userdata/dxu/project_np/code/babble");
end

thisFile = [mfilename("fullpath") '.m'];
tempFile = erase(thisFile, '_jobs');
tempTxt = fileread(tempFile);

anaDir = LMV.Data.GetAnalysisDir("coding", "pop_sen4_ecoc");
[~, baseFileName] = fileparts(tempFile);

%% Create scripts for different model types

regions = "all"; % ["all", "mPrCG"];
groups = ["feedback", "mirror", "bridge"];

scriptFiles = cell(numel(regions), numel(groups));

for i = 1 : numel(regions)
    for j = 1 : numel(groups)
        % Modify set
        [a, b] = regexp(tempTxt, "(?<=regionName = "")[^\r]+(?="";\r)", "once");
        fprintf("\nReplace <%s> ", tempTxt(a:b));
        
        scriptTxt = [tempTxt(1:a-1) char(regions(i)) tempTxt(b+1:end)];
        [a, b] = regexp(tempTxt, "(?<=regionName = "")[^\r]+(?="";\r)", "once");
        fprintf(" with <%s>.\n", scriptTxt(a:b));
        
        % Modify target
        [a, b] = regexp(scriptTxt, "(?<=groupName = "")[^\r]+(?="";\r)", "once");
        fprintf("Replace <%s> ", scriptTxt(a:b));
        
        scriptTxt = [scriptTxt(1:a-1) char(groups(j)) scriptTxt(b+1:end)];
        [a, b] = regexp(scriptTxt, "(?<=groupName = "")[^\r]+(?="";\r)", "once");
        fprintf(" with <%s>.\n", scriptTxt(a:b));
        
        % Write file
        scriptFile = sprintf("%s_%s_%s.m", baseFileName, regions(i), groups(j));
        scriptFiles{i,j} = fullfile(anaDir, scriptFile);
        fileID = fopen(scriptFiles{i,j}, 'w');
        fprintf(fileID, '%s', scriptTxt);
        fclose(fileID);
    end
end

scriptFiles = permute(scriptFiles, [2 1]);

%% Send job commands to console

qList = {'pia-batch.q', 'skull-batch.q', 'spirit-batch'};
qList = qList([1 2 3 1 3 1 3]);
qList = repmat(qList, 1, 100);

for k = 1 : numel(scriptFiles)
    jobCmd = NP.Data.GetJobCmd(scriptFiles{k}, true, 'Queue', qList{k}, 'RAM', 8*8);
    fprintf("\n%s\n", jobCmd);
end

return
%% Create scripts for different model types

regions = "mPrCG";
groups = ["bridge", "mirror"];
shuffles = ["false", "true"];

scriptFiles = cell(numel(regions), numel(groups), numel(shuffles));

for i = 1 : numel(regions)
    for j = 1 : numel(groups)
        for k = 1 : numel(shuffles)
            % Modify set
            [a, b] = regexp(tempTxt, "(?<=regionName = "")[^\r]+(?="";\r)", "once");
            fprintf("\nReplace <%s> ", tempTxt(a:b));
            
            scriptTxt = [tempTxt(1:a-1) char(regions(i)) tempTxt(b+1:end)];
            [a, b] = regexp(tempTxt, "(?<=regionName = "")[^\r]+(?="";\r)", "once");
            fprintf(" with <%s>.\n", scriptTxt(a:b));
            
            % Modify target
            [a, b] = regexp(scriptTxt, "(?<=groupName = "")[^\r]+(?="";\r)", "once");
            fprintf("Replace <%s> ", scriptTxt(a:b));
            
            scriptTxt = [scriptTxt(1:a-1) char(groups(j)) scriptTxt(b+1:end)];
            [a, b] = regexp(scriptTxt, "(?<=groupName = "")[^\r]+(?="";\r)", "once");
            fprintf(" with <%s>.\n", scriptTxt(a:b));
            
            % Modify shuffling
            [a, b] = regexp(scriptTxt, "(?<='Shuffle', )[^,]+(?=,)", "once");
            fprintf("Replace <%s> ", scriptTxt(a:b));
            
            scriptTxt = [scriptTxt(1:a-1) char(shuffles(k)) scriptTxt(b+1:end)];
            [a, b] = regexp(scriptTxt, "(?<='Shuffle', )[^,]+(?=,)", "once");
            fprintf(" with <%s>.\n", scriptTxt(a:b));
            
            if shuffles(k) == "false"
                scriptFile = sprintf("%s_%s", regions(i), groups(j));
            else
                scriptFile = sprintf("%s_%s_null", regions(i), groups(j));
            end
            
            % Write file
            scriptFiles{i,j,k} = strrep(thisFile, 'jobs', scriptFile);
            % fileID = fopen(scriptFiles{i,j}, 'w');
            % fprintf(fileID, '%s', scriptTxt);
            % fclose(fileID);
        end
    end
end

scriptFiles = permute(scriptFiles, [3 2 1]);

