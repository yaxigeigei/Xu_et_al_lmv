%% Population sentence classification

anaDir = LMV.Data.GetAnalysisDir("coding", "pop_sen4_ecoc");

%% Load data

% SCA input for finding the factors that zscored the input for SCA
load(fullfile(LMV.Data.GetAnalysisDir, "data", "ce_m2_ex3_sentence-avg.mat"), 'ce');
respTb = ce.GetTable("resp");
R = cell2mat(respTb{:,2:end});
Rm = mean(R, "omitmissing");
Rsd = std(R, "omitmissing");
Rsd = max(Rsd, eps);

% Single-trial responses
srcTb = LMV.Data.FindSource([]);
cacheDir = LMV.Data.GetAnalysisDir('data', 'resp_m2_ex3_sen4_loo');
sData = arrayfun(@load, fullfile(cacheDir, srcTb.recId+".mat"));

% Task phase responsiveness
sTest = LMV.Resp.LoadPhaseResponseTest();

%% 

for i = 1 : numel(sData)
    % Mark task responsive units
    [~, I] = MMath.SortLike(sTest.clusTb.clusId, sData(i).clusTb.clusId);
    sData(i).clusTb.isResp = sTest.clusTb.tId1(I) ~= "none";
    
    % Downsample responses
    R = sData(i).resp;
    R = smoothdata(R, 3, "movmean", 20);
    R = R(:,:,5:10:end);
    
    % Normalize responses
    [~, I] = MMath.SortLike(ce.clusTb.clusId, sData(i).clusTb.clusId);
    R = (R - Rm(I)) ./ Rsd(I);
    sData(i).respDS = R;
end

%% Load SCA output

regions = [LMV.Param.regions "mPrCG"]';
conds = [LMV.Param.regions "mPrCG_no-bridge"]';
nComp = 12;
scaTbs = LMV.SCA.LoadResults(fullfile(LMV.Data.GetAnalysisDir, "pop_dynamics", "sca", "computed_sca", "sca_"+conds+".mat"), nComp);

regTb = table;
regTb.cond = conds;
regTb.region = regions;
regTb = [regTb vertcat(scaTbs{:})];
regTb = LMV.SCA.EnrichResultTable(regTb);

%% Classification

% Get SCA subspace
indComp = [4 6:12];
indComp = 1:12;
U = regTb.U{1}(:,indComp);
bu = regTb.b_u{1}(indComp);

% Make trial labels
yTrain = repelem(1:4, 7)';
yTest = (1:4)';

% Compute cross-time classification
isReg = srcTb.Region == "mPrCG";
nBoot = 50;
nTime = size(sData(1).respDS, 3);
rxBoot = zeros(nTime, nTime, nBoot);
for n = 1 : nBoot
    % Construct training and testing arrays
    xTrain = arrayfun(@(x) x.respDS(x.indTrain(:,n), x.clusTb.isResp, :), sData(isReg), 'Uni', false);
    xTest = arrayfun(@(x) x.respDS(x.indTest(:,n), x.clusTb.isResp, :), sData(isReg), 'Uni', false);
    xTrain = cat(2, xTrain{:});
    xTest = cat(2, xTest{:});
    
    % Project to sentence-selective subspace
    zTrain = pagemtimes(xTrain, U) + bu;
    zTest = pagemtimes(xTest, U) + bu;
    
    % Normalize predictors in each time bin
    zAll = [zTrain; zTest];
    for t = 1 : size(zAll,3)
        [~, c, k] = MMath.Normalize(zAll(:,:,t), 'zscore');
        zTrain(:,:,t) = (zTrain(:,:,t)-c) ./ k;
        zTest(:,:,t) = (zTest(:,:,t)-c) ./ k;
    end
    
    % 
    str = sprintf('Iteration %i/%i\n', n, nBoot);
    fprintf(str);
    parfor a = 1 : nTime
        for b = 1 : nTime
            mdl = fitcecoc(zTrain(:,:,a), yTrain, 'Coding', 'onevsall');
            rxBoot(a,b,n) = 1 - loss(mdl, zTest(:,:,b), yTest);
        end
    end
    if n < nBoot
        fprintf(repmat('\b', [1 numel(str)]));
    end
end

%% 

f = MPlot.Figure(53201); clf
tl = tiledlayout("flow");
tl.Padding = "compact";
ax = nexttile;

t = ce.GetTable("resp").time{1}(5:10:end);
[rxMean, ~, rxSE] = MMath.MeanStats(rxBoot, 3);
imagesc(t, t, rxMean); hold on

tt = ce.GetTable("taskTime");
wStim = [tt.stimMatchOn(1) tt.stimMatchOff(1)];
wProd = [tt.prodMatchOn(1) tt.prodMatchOff(1)];
plot([wStim; wStim], ax.YLim', 'Color', LMV.Param.GetTaskPhaseColors("stim"));
plot(ax.YLim', [wStim; wStim], 'Color', LMV.Param.GetTaskPhaseColors("stim"));
plot([wProd; wProd], ax.YLim', 'Color', LMV.Param.GetTaskPhaseColors("prod"));
plot(ax.YLim', [wProd; wProd], 'Color', LMV.Param.GetTaskPhaseColors("prod"));

axis(ax, 'square', 'equal', 'tight');
b = colorbar;
b.Label.String = "Accuracy (frac.)";
ax.CLim(1) = 0.25;
ax.CLim(2) = 0.55;
ax.YLabel.String = "Trained from (s)";
ax.XLabel.String = "Tested on (s)";
MPlot.Axes(ax);
MPlot.Paperize(f, .8, .7);



