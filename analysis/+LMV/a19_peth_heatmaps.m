%% 

mdlName = "smooth_lm";
typeName = "bridge";
figDir = fullfile(LMV.Data.GetAnalysisDir, "linker", mdlName);

%% Find linker units

clusTb = LMV.Linker.LoadClusTb(mdlName);
cTb = clusTb(clusTb.hcGroup==typeName, :);

%% Load cached profile data

recIdList = unique(cTb.recId);
cacheDir = LMV.Data.GetAnalysisDir("linker", mdlName, "computed_profile");
sPf = cell(size(recIdList));
for recIdx = 1 : numel(sPf)
    disp(recIdList(recIdx));
    cacheFile = fullfile(cacheDir, sprintf("%s.mat", recIdList(recIdx)));
    sPf{recIdx} = load(cacheFile);
end

%% Find best linked sentences and positions for each unit

for recIdx = 1 : height(sPf)
    % Get variables
    s = sPf{recIdx};
    recId = s.recId;
    ce = s.ce;
    triTb = s.triTb;
    
    % Update linking positions
    LMV.Linker.LM.ComputePredictedPeth(ce);
    pkTb = LMV.Linker.LM.FindPredPeaks(ce);
    
    % Add peak linking info to brTb
    cid = intersect(pkTb.clusId, cTb.clusId);
    for u = 1 : numel(cid)
        isUnit = cTb.clusId==cid(u);
        isPk = pkTb.clusId==cid(u);
        pkTbSub = sortrows(pkTb(isPk,:), "score", "descend");
        cTb.pkTb{isUnit} = pkTbSub;
    end
end

%% Extract response timeseries from M2 unit cache

cTb.uCache = NP.Unit.LoadUnitCache(cTb.clusId, 'DataSource', 'm2');

lmvDur = cTb.uCache{1}.tt.prodMatchOff(1) - cTb.uCache{1}.tt.stimMatchOn(1);
rsPad = [-1.5 .5];
rsBinSize = 0.01;
rsKerSize = 0.1;

t = MMath.BinEdges2Centers(rsPad(1) : rsBinSize : lmvDur+rsPad(2));
isBase = t < 0 | t > lmvDur;

for i = 1 : height(cTb)
    s = cTb.uCache{i};
    pkTb = cTb.pkTb{i};
    
    % Find positions with highest scores
    isOut = isoutlier(pkTb.score, "median", ThresholdFactor=3);
    if any(isOut)
        pkTb = pkTb(isOut,:);
        [~, indPk] = maxk(pkTb.score, 3);
    else
        % continue
        [~, indPk] = max(pkTb.score);
    end
    
    % Compute response timeseries between pairs of positions
    nRep = NaN(size(indPk));
    srCell = cell(size(indPk));
    
    for j = 1 : numel(indPk)
        % Find the sentence
        k = indPk(j);
        isSen = s.tv.stimId == pkTb.stimId(k);
        
        % Find the time window of bridging
        phSeq = Cut(pkTb.seqTge{k});
        phStim = Cut(Cut(s.tt.stim(isSen)));
        phProd = Cut(Cut(s.tt.prod{isSen}));
        
        phSeqStr = erase(phSeq.GetParentLabel, digitsPattern);
        phStimStr = erase(phStim.GetParentLabel, digitsPattern);
        phProdStr = erase(phProd.GetParentLabel, digitsPattern);
        
        nPh = numel(phSeq);
        [~, a] = max(smoothdata(matches(phStimStr, phSeqStr), 'movmean', nPh));
        [~, b] = max(smoothdata(matches(phProdStr, phSeqStr), 'movmean', nPh));
        lkWin = double([phStim(a) phProd(b)]);
        
        % Resample spike times
        rsEdges = lkWin(1)+rsPad(1) : rsBinSize : lkWin(1)+lmvDur+rsPad(2);
        st = s.st{isSen};
        sr = MNeuro.MeanEventRate(st, rsEdges);
        sr = MNeuro.Filter1(sr, 1/rsBinSize, 'gaussian', rsKerSize);
        
        nRep(j) = numel(st);
        srCell{j} = sr;
    end
    
    % Find the strongest bridging
    [~, I] = max(cellfun(@(x) mean(x(~isBase))/mean(x(isBase)), srCell));
    
    % brTb.brWin(i,:) = brWins{I};
    cTb.nRep(i) = nRep(I);
    cTb.spikeRate{i} = srCell{I};
end

%% Plot heatmap of bridging responses

cTb2 = cTb(cTb.nRep>0,:);

% Get responses
R = cat(2, cTb2.spikeRate{:});
R = R - median(R(isBase,:), 1);
R = MMath.Normalize(R, "max");

% Sort units
% D = pdist(R', "seuclidean");
% tree = linkage(D);
% uInd = optimalleaforder(tree, D);
[~, uInd] = LMV.NMF.SortPatterns(R, 'csc');
% [~, uInd] = sort(sum(R), 'descend');

% Make plot
f = MPlot.Figure(783); clf
ax = nexttile;
imagesc(t, [], R(:,uInd)');
colorbar;
ax.CLim(1) = 0;
ax.YTick = 1:height(cTb2);
ax.YTickLabel = "u"+cTb2.clusId(uInd);
ax.XLabel.String = "Aligned time from bridge onsets (s)";
ax.Title.String = "Evoked bridge responses";
MPlot.Axes(ax);
MPlot.Paperize(f, 1, 0.5);
exportgraphics(f, fullfile(figDir, "lmv_heatmap_"+typeName+".png"));

return
%% Case studies of delay responses








