%% Extract spike-triggered features

anaDir = LMV.Data.GetAnalysisDir("st_feat");
% srcTb = LMV.Data.FindSource([]);

%% Load ce

ceSearch = MBrowse.Dir2Table(fullfile(anaDir, "ce", "*_ce.mat"));
ceArray = NP.CE.LoadSession(fullfile(ceSearch.folder, ceSearch.name));

%% 

cacheDir = fullfile(anaDir, "extracted_feat");
if ~exist(cacheDir, 'dir')
    mkdir(cacheDir);
end

phases = ["stim", "prod"];

feats = "mic";

feats = [{'env', 'peakEnv', 'peakRate', 'rF0', 'drF0'}, ...
    {'ja', 'la', 'pro', 'ttcc', 'tbcc', 'tdcc', 'ttcl', 'tbcl', 'v_x', 'v_y'}, ...
    NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.consonants, ...
    ];


for i = 4 %1 : numel(ceArray)
    % 
    ce = ceArray(i).Duplicate;
    tSample = ce.userData.ops.rsBinSize;
    
    pn = phases(2);
    
    % Get the mask for the task phase of interest
    %   Use time shifts to includes a period before or after to account for response lag
    if pn == "stim"
        [~, M, mNames] = ce.GetArray('feat', [], pn, 'TimeShifts', 0 : tSample : 0.2);
    elseif pn == "prod"
        [~, M, mNames] = ce.GetArray('feat', [], pn, 'TimeShifts', -0.2 : tSample : 0);
    end
    M = any(M, 2);
    
    % Get features
    [T, F, fNames] = ce.GetArray('feat', [], feats);
    F(~M,:) = NaN;
    
    [~, pitchInd] = MMath.SortLike(string(feats), {'rF0', 'drF0'});
    if ~isempty(pitchInd)
        pitch = F(:,pitchInd);
        isOut = isoutlier(pitch, 1);
        pitch(isOut) = NaN;
        F(:,pitchInd) = pitch;
    end
    
    if ismember(feats, ["mic" "speaker1"])
%         F = MMath.Normalize(F, 'minmax');
        F(isnan(F)) = -80;
    else
        F = MMath.Normalize(F, 'max');
        F(isnan(F)) = 0;
    end
%     imagesc(F');
%     ax = gca;
%     ax.YTick = 1:numel(fNames);
%     ax.YTickLabel = fNames;
    
    tShift = -0.5 : tSample : 0.5;
    indShift = round(tShift ./ tSample);
    FF = MMath.RepShiftMat(F, indShift);
    
    % Get responses
    sTest = LMV.Resp.LoadSentenceSelectTest(NP.SE.GetID(ce));
    uInd = find(sTest.sigTb.(pn));
    uInd = 23; % u245
%     uInd = 8; % u284
    [~, R, uNames] = ce.GetArray('resp', [], uInd+1);
    R(~M,:) = 0;
end

%% 

[stFv, N] = LMV.ST.TriggerFeatures(R, FF);

Z = run_umap(stFv, 'metric', 'euclidean', 'verbose', 'text');

%% 

MPlot.Figure(7456); clf
LMV.ST.PlotFeatureUMAP(Z, stFv, fNames, tShift);



