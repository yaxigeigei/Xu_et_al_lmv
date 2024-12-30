%% Read template script

% Command
% /data_store2/MATLAB/R2022a/bin/matlab -nodesktop -nosplash -r "run('/userdata/dxu/project_np/code/babble/f4_fit_trf_jobs.m'); exit"

if ~ispc && ~ismac
    addpath(genpath("/userdata/dxu/project_np/code/third_party_matlab"));
    addpath("/userdata/dxu/project_np/code/babble");
end

thisFile = [mfilename("fullpath") '.m'];
tempFile = erase(thisFile, '_jobs');
tempTxt = fileread(tempFile);

%% Create scripts for different model types

sets = ["combo4akt", "combo4pm", "strf"];
targets = ["stim", "prod", "feedback"];
scriptFiles = cell(numel(sets), numel(targets));

for i = 1 : numel(sets)
    for j = 1 : numel(targets)
        % Modify set
        [a, b] = regexp(tempTxt, "(?<=sets = \["")[^\r]+(?=""\];\r)", "once");
        fprintf("\nReplace <%s> ", tempTxt(a:b));
        
        scriptTxt = [tempTxt(1:a-1) char(sets(i)) tempTxt(b+1:end)];
        [a, b] = regexp(scriptTxt, "(?<=sets = \["")[^\r]+(?=""\];\r)", "once");
        fprintf(" with <%s>.\n", scriptTxt(a:b));
        
        % Modify target
        [a, b] = regexp(scriptTxt, "(?<=phases = \["")[^\r]+(?=""\];\r)", "once");
        fprintf("Replace <%s> ", scriptTxt(a:b));
        
        scriptTxt = [scriptTxt(1:a-1) char(targets(j)) scriptTxt(b+1:end)];
        [a, b] = regexp(scriptTxt, "(?<=phases = \["")[^\r]+(?=""\];\r)", "once");
        fprintf(" with <%s>.\n", scriptTxt(a:b));
        
        % Write file
        scriptFiles{i,j} = strrep(thisFile, 'jobs', sets(i)+"_"+targets(j));
        fileID = fopen(scriptFiles{i,j}, 'w');
        fprintf(fileID, '%s', scriptTxt);
        fclose(fileID);
    end
end

%% Send job commands to console

qList = {'pia-batch.q', 'skull-batch.q', 'spirit-batch'};
qList = qList([1 2 3 1 3 1 3]);
qList = repmat(qList, 1, 100);

for k = 1 : numel(scriptFiles)
    jobCmd = NP.Data.GetJobCmd(scriptFiles{k}, true, 'Queue', qList{k}, 'RAM', 8*8);
    fprintf("\n%s\n", jobCmd);
end
