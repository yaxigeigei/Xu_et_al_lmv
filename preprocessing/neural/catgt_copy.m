%% Paths

add_np_paths;

rawRoot = NP.Data.GetRawRoot;
preprocRoot = NP.Data.GetPreprocRoot;

%% Recording info

if ~exist('recId', 'var')
    subjId = "NP47";
    blockId = "B3";
    recId = subjId + '_' + blockId;
end

%% Run CatGT

isCatGT = true;

% Locate the raw ap.bin file
% e.g. "/datastore_spirit/human/Neuropixels/NP41/NP41_B1_g0/NP41_B1_g0_imec0/NP41_B1_g0_t0.imec0.ap.bin"
apPattern = fullfile(rawRoot, subjId, ...
    strjoin([recId "g*"], '_'), ...
    strjoin([recId "g*_imec*"], '_'), ...
    "*.ap.bin");

apSearch = MBrowse.Dir2Table(apPattern);

if height(apSearch) < 1
    error("No file matches '%s'", apPattern);
% else
%     apFile = fullfile(apSearch.folder, apSearch.name);
end

% Find gates, probes, and associated paths
nApFiles = height(apSearch);
gates = strings(1, nApFiles);
probes = strings(1, nApFiles);
for i = 1 : nApFiles
    apFile = string(fullfile(apSearch.folder{i}, apSearch.name{i}));
    [probeDir, apName] = fileparts(apFile);
    [runDir, probeDirName] = fileparts(probeDir);
    gates(i) = regexp(probeDirName, '(?<=_g)\d+', 'match');
    probes(i) = regexp(probeDirName, '(?<=_imec)\d+', 'match');
end
gates = strjoin(unique(gates), ',');
probes = strjoin(unique(probes), ',');

[subjDir, runDirName] = fileparts(runDir);

if ~isCatGT
    % Simply copy the run folder to userdata
    destDir = fullfile(preprocRoot, recId, runDirName);
    copyfile(runDir, destDir);
else
    % Run CatGT and output result to userdata
    isLinux = strcmp(computer, 'GLNXA64');
    if isLinux
        cgtPath = "/userdata/dxu/project_np/code/third_party_others/CatGT-linux/runit.sh";
    else
        cgtPath = fullfile(pkgDir, "CatGT-win/CatGT.exe");
    end
    
    destDir = fullfile(preprocRoot, recId, "sglx");
    if ~exist(destDir, 'dir')
        mkdir(destDir);
    end
    
    cgtArgs = strjoin([ ...
        "-dir=" + subjDir, ...          % subject directory of source data, e.g. '/datastore_spirit/human/Neuropixels/NP11'
        "-run=" + recId, ...            % run name, e.g. 'NP11_B4'
        "-g=" + gates, ...              % gates
        "-t=" + "0,0", ...              % trigger range
        "-zerofillmax=0", ...           % zero for not leaving any gap when concatenating
        "-prb=" + probes, ...           % probes
        "-ap -lf -ni", ...              % stream
        "-gblcar", ...                  % ap common average referencing
        "-gfix=0.40,0.10,0.02", ...     % ap artifact removal: ||amp(mV)||, ||slope(mV/sample)||, ||noise(mV)||
        "-prb_fld -out_prb_fld", ...    % a folder per probe for both input and output
        "-dest=" + destDir ...          % directory of output data, e.g. '/userdata/dxu/project_np/preproc/NP11/sglx'
        ], ' ');
    
    if isLinux
        cgtCmd = strjoin(["bash", cgtPath, cgtArgs], " ");
    else
        cgtCmd = strjoin([cgtPath, cgtArgs], " ");
    end
    
    disp('CatGT command:')
    disp(cgtCmd);
    
    system(cgtCmd);
    
    % Copy nidq.bin file
    % e.g. "/datastore_spirit/human/Neuropixels/NP42/NP42_B1_g0/NP42_B1_g0_t0.nidq.bin"
    niPattern = fullfile(rawRoot, subjId, recId+"_g*", "*.nidq.bin");
    niSearch = MBrowse.Dir2Table(niPattern);
    niSrcBin = fullfile(niSearch.folder{1}, niSearch.name{1});
    
    niPattern = fullfile(destDir, "catgt_"+recId+"_g*", "*.nidq.meta");
    niSearch = MBrowse.Dir2Table(niPattern);
    niTarBin = fullfile(niSearch.folder{1}, niSearch.name{1});
    niTarBin = strrep(niTarBin, ".meta", ".bin");
    
    if ~exist(niTarBin, 'file')
        copyfile(niSrcBin, niTarBin);
    end
end
