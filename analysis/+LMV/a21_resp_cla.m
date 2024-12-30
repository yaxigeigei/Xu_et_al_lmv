%% Population sentence classification

if ~ispc && ~ismac
    addpath(genpath("/userdata/dxu/project_np/code/third_party_matlab"));
    addpath("/userdata/dxu/project_np/code/babble");
    addpath("/userdata/dxu/project_np/code/preproc");
end

anaDir = LMV.Data.GetAnalysisDir("coding", "pop_sen4_ecoc");

regionName = "all";
groupName = "bridge";

%% Load data

% Single-trial responses
if regionName == "all"
    srcTb = LMV.Data.FindSource([]);
else
    srcTb = LMV.Data.FindSource([], 'Region', regionName);
end
cacheDir = LMV.Data.GetAnalysisDir('data', 'resp_m2_detrend_ex3_sen4_loo');
sData = arrayfun(@load, fullfile(cacheDir, srcTb.recId+".mat"));

% Task phase responsiveness
sTest = LMV.Resp.LoadPhaseResponseTest();

% Linker groups
lkTb = LMV.Linker.LoadClusTb();

%% Prepare data

switch groupName
    case {"mirror", "bridge", "feedback"}
        uid = lkTb.clusId(lkTb.hcGroup==groupName);
        if groupName == "bridge"
            uid = [uid; LMV.Linker.bridgeFN(:)];
        end
    case "responsive"
        uid = sTest.clusTb.clusId(sTest.clusTb.tId1 ~= "none");
    otherwise
        error("'%s' is not a valid group name.", groupName);
end

for i = 1 : numel(sData)
    % Mark units to be included
    sData(i).clusTb.isUnit = ismember(sData(i).clusTb.clusId, uid);
    
    % Downsample responses
    R = sData(i).resp;
    R = smoothdata(R, 3, "movmean", 20);
    R = R(:,:,5:10:end);
    
    % Normalize responses in each time bin for each unit
    for t = 1 : size(R,3)
        R(:,:,t) = MMath.Normalize(R(:,:,t), 'zscore');
    end
    sData(i).X = R;
end

%% Classification

% Compute
nBoot = 250;

fprintf("\nCross-time classification\n");
xrBoot = LMV.Decode.ClassifyXT(sData, 'NBoot', nBoot, 'Shuffle', false);

fprintf("\nCross-time classification with shuffled labels\n");
xrNullBoot = LMV.Decode.ClassifyXT(sData, 'NBoot', nBoot, 'Shuffle', true);

% Save results
fileName = sprintf("rx_%s_%s.mat", regionName, groupName);
save(fullfile(anaDir, fileName), 'regionName', 'groupName', 'xrBoot', 'xrNullBoot');
