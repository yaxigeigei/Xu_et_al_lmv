%% Distributions of linking positions

mdlName = "smooth_lm";
srcTb = LMV.Data.FindSource([]);

%% Load data

% Cached profiles
cacheDir = LMV.Data.GetAnalysisDir("linker", mdlName, "computed_profile");
sPf = cell(size(srcTb.recId));
for recIdx = 1 : height(srcTb)
    % Load cache
    recId = NP.SE.GetID(srcTb.name{recIdx});
    disp(recId);
    cacheFile = fullfile(cacheDir, sprintf("%s.mat", recId));
    s = load(cacheFile, 'ce', 'recId');
    
    % Update linking positions
    LMV.Linker.LM.ComputeLinkingScores(s.ce);
    s.posTb = LMV.Linker.LM.FindLinkedPositions(s.ce);
    sPf{recIdx} = s;
end

% Linker clusTb from hierarchical clustering
[clusTb, clusTbPath] = LMV.Linker.LoadClusTb(mdlName);

%% Add linking positions to linker clusTb

lkIds = clusTb.clusId(ismember(clusTb.hcGroup, LMV.Linker.types));

for recIdx = 1 : height(sPf)
    % Get variables
    s = sPf{recIdx};
    recId = s.recId;
    ce = s.ce;
    posTb = s.posTb;
    
    % Add posTb
    cid = intersect(posTb.clusId, lkIds);
    for u = 1 : numel(cid)
        % Find candidate positions for the unit
        isPos = posTb.clusId==cid(u);
        pTb = sortrows(posTb(isPos,:), "score", "descend");
        
        % Keep linked positions with outlier scores
        isOut = isoutlier(pTb.score, "median", ThresholdFactor=3);
        pTb.isOut = isOut;
        if any(isOut)
            pTb = pTb(isOut,:);
        else
            pTb = pTb(1,:);
        end
        
        % Save posTb to clusTb
        isUnit = clusTb.clusId==cid(u);
        clusTb.posTb{isUnit} = pTb;
    end
end

%% Save modified linker clusTb

save(clusTbPath, 'clusTb');
