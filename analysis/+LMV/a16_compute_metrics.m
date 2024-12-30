%% Compute metrics and stats on 

anaDir = LMV.Data.GetAnalysisDir("pop_dynamics", "sca");
srcTb = LMV.Data.FindSource([]);

cachePath = fullfile(anaDir, "regTb.mat");
if isfile(cachePath)
    fprintf("Cached regTb is found at the following location and will not be overwritten.\n%s\n", cachePath)
    return
end

%% Single-trial responses

% Load SCA input to find factors that zscored the input for SCA
load(fullfile(LMV.Data.GetAnalysisDir, "pop_dynamics", "ce_m2_ex3_sentence-avg.mat"), 'ce');
respTb = ce.GetTable("resp");
R = cell2mat(respTb{:,2:end});
Rm = mean(R, "omitmissing");
Rsd = std(R, "omitmissing");

% Load single-trial responses
sTrial = load(fullfile(LMV.Data.GetAnalysisDir, "data", "resp_m2_ex3_rep8_cat.mat"));

% Normalize single-trial responses using the same factors
[~, I] = MMath.SortLike(sTrial.clusTb.clusId, ce.clusTb.clusId);
R = sTrial.resp(:,I,:,:);
R = R - Rm;
R = R ./ Rsd;

%% Load unit info from pickled DataFrame

clusDf = py.pandas.read_pickle(fullfile(anaDir, "df", "clus.pkl"));
clusTb = table(clusDf);
[~, I] = MMath.SortLike(clusTb.clusId, ce.clusTb.clusId);
clusTb = clusTb(I,:);

%% Load SCA output

regions = [LMV.Param.regions "mPrCG"]';
conds = [LMV.Param.regions "mPrCG_no-bridge"]';
nComp = 12;
scaTbs = LMV.SCA.LoadResults(fullfile(anaDir, "computed_sca", "sca_"+conds+".mat"), nComp);

regTb = table;
regTb.cond = conds;
regTb.region = regions;
regTb = [regTb vertcat(scaTbs{:})];
regTb = LMV.SCA.EnrichResultTable(regTb);

%% Project single-trial responses

% Bootstrap options
nShuffle = 50;
rng(61);

% Find responsive units
isResp = any(clusTb{:,4:8}, 2);

for i = 1 : height(regTb)
    disp(regTb.cond(i));
    
    % Get single-trial responses
    isUnit = isResp & clusTb.region==regTb.region(i);
    if regTb.cond(i)=="mPrCG_no-bridge"
        isUnit = isUnit & clusTb.hcGroup~="bridge";
    end
    r4 = R(:,isUnit,:,:);
    
    % Direct projections
    r3 = MMath.CombineDims(r4, [1 3]);
    zzTr = pagemtimes(r3, regTb.U{i}) + regTb.b_u{i};
    
    % Shuffled projections
    sz = size(r4);
    zzNull = zeros(sz(1)*sz(3), size(zzTr,2), nShuffle); % time*sen-by-comp-by-iter
    for n = 1 : nShuffle
        r3 = MMath.CombineDims(r4, [3 4]);          % sen*rep-by-time-by-unit
        ind = randsample(sz(3)*sz(4), sz(3)*sz(4));
        r3 = r3(ind,:,:);
        r3 = reshape(r3, sz([3 4 1 2]));            % sen-by-rep-by-time-by-unit
        r3 = MMath.CombineDims(r3, [3 1]);          % time*sen-by-rep-by-unit
        r2 = squeeze(mean(r3, 2, "omitmissing"));   % time*sen-by-unit
        zzNull(:,:,n) = r2*regTb.U{i} + regTb.b_u{i};
    end
    zzNullCI = prctile(zzNull, [2.5 97.5], 3);
    
    regTb.Ztr{i} = zzTr;
    regTb.ZnullCI{i} = zzNullCI;
end

%% Cache result

save(cachePath, "regTb", "clusTb");

return
%% 

X = ssSCA(2).X;
sComp = ssSCA(2).decomp_12;
Z = X*sComp.U + sComp.b_u;

figure(1); clf
plot([Z(:,1) sComp.Z(:,1)])
plot([Z(:,1) sComp.Z(:,1)])

figure(2); clf
histogram(Z(:,1)-sComp.Z(:,1))
