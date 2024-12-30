%% Determine unit responsiveness using ZETA-test

anaDir = LMV.Data.GetAnalysisDir('phase_resp', 'zeta');
srcTb = LMV.Data.FindSource([]);

%% 

addpath(genpath(fullfile(getenv('NP_ROOT'), "code\third_party_matlab\zetatest")));

% Load M2 se
sePaths = fullfile(NP.Data.GetAnalysisRoot, "data", "se_m2", srcTb.recId + "_se_m2.mat");

% Specify task phases to test
phaseNames = {'atten', 'stim', 'delay', 'init', 'prod', 'iti'};

% Make folder to save cache files
cacheDir = fullfile(anaDir, "computed_tests");
if ~exist(cacheDir, 'dir')
    mkdir(cacheDir);
end

% Compute Zeta responsiveness for each recording
for i = 1 : numel(sePaths)
    % Check cache status
    recId = srcTb.recId{i};
    cachePath = fullfile(cacheDir, recId+"_clusTb.mat");
    if exist(cachePath, 'file')
        fprintf("\nZeta test result of %s has been computed and cached\n", recId);
        continue
    end
    
    % Load M2 se
    se = NP.SE.LoadSession(sePaths{i});
    se.userData.ksMeta.clusTb = NP.Unit.AddRecMeta(se, se.userData.ksMeta.clusTb);
    
    % Compute Zeta responsiveness
    fprintf("\nCompute task phase response for %s\n", recId);
    zTb = NP.Resp.PhaseZeta(se, phaseNames);
    
    % Save result
    clusTb = [NP.Unit.GetClusTb(se) zTb];
    save(cachePath, 'clusTb', 'phaseNames');
end

%% Save concatenated table

% Load test tables
nRec = height(srcTb);
clusTbs = cell(nRec, 1);
for i = 1 : nRec
    load(fullfile(cacheDir, srcTb.recId{i}+"_clusTb.mat"));
    clusTbs{i} = clusTb;
end
clusTb = cat(1, clusTbs{:});

% Remove bulky data
clusTb.spikeTimes = [];
clusTb(:, phaseNames+"Win") = [];

% Get pvals
pTb = clusTb(:, phaseNames);

% Construct task phase index table
phaseInd = [5 1 2 3 4 5]';
iTb = pTb;
for i = 1 : width(iTb)
    iTb.(i)(:) = phaseInd(i);
end

% Find best task phase responsiveness
P = pTb{:,:};
I = iTb{:,:};
[zetaP, k] = min(P, [], 2);
k = sub2ind(size(P), (1:size(P,1))', k);
zetaId = I(k);

% Make output table
clusTb.zetaId = zetaId;
clusTb.zetaP = zetaP;

save(fullfile(anaDir, 'computed_test_clusTb.mat'), 'clusTb');

return
%% Plot speech and rasters for inspection

% Load se
m2epPaths = fullfile(NP.Data.GetAnalysisRoot, "data", "computed_se_m2", srcTb.recId + "_se_m2_ep.mat");
seM2ep = NP.SE.LoadSession(m2epPaths);

% Compute sentence PETH
senTbs = arrayfun(@NP.TaskBaseClass.SplitBySentence, seM2ep, 'Uni', false);
[ceCell, senTbs] = cellfun(@(x) LMV.SE.ComputeSentencePETH(x), senTbs, 'Uni', false);

%% 

groupNames = {'lr-only', 'zeta-only', 'both'};

fSig = @(x) any(x < 0.01, 2);

for i = 1 : numel(uTbs)
    for t = 1 : numel(phaseNames)
        vn = phaseNames{t} + "P";
        if phaseNames{t} == "engage"
            continue
        end
        for g = 1 : numel(groupNames)
            % Add modulation indices
            uSubTb = uTbs{i};
            uSubTb.mi = ceCell{i}.userData.ksMeta.clusTb.mi4;
            
            % Select units
            switch groupNames{g}
                case 'lr-only'
                    m = fSig(rTbs{i}.(vn)) & ~fSig(zTbs{i}.(vn));
                case 'zeta-only'
                    m = ~fSig(rTbs{i}.(vn)) & fSig(zTbs{i}.(vn));
                case 'both'
                    m = fSig(rTbs{i}.(vn)) & fSig(zTbs{i}.(vn));
            end
            uSubTb = uSubTb(m,:);
            
            % Sort units by modulation index
            uSubTb = sortrows(uSubTb, 'mi', 'descend');
            
            % Sort units by depth within each page
            upp = 15;
            nPages = ceil(height(uSubTb) / upp);
            pInd = repelem(1:nPages, upp);
            uSubTb.pInd = pInd(1:height(uSubTb))';
            uSubTb = sortrows(uSubTb, {'pInd', 'depth'}, 'ascend');
            
            % 
            uIdPage = NaN(upp, nPages);
            uIdPage(1:height(uSubTb)) = uSubTb.clusId;
            for p = 1 : min(nPages, 20)
                f = MPlot.Figure(22300); clf
                LMV.Plot.Overview(senTbs{i}, ceCell{i}, 'TaskPhase', 'full', ...
                    'UnitIds', uIdPage(:,p), 'Page', p, 'Folder', fullfile(anaDir, 'compare_nlr-zeta', [phaseNames{t} '_' groupNames{g}]));
            end
        end
    end
end

%% 

groupNames = {'sig0', 'sig1', 'sig2', 'sig3'};

fSig = @(x) max(x, [], 2);

for i = 1 : numel(uTbs)
    pn = 'engage';
    for g = 1 : numel(groupNames)
        % Add modulation indices
        uSubTb = uTbs{i};
        uSubTb.mi = ceCell{i}.userData.ksMeta.clusTb.mi4;

        % Select units
        m = fSig(zTbs{i}.(pn)) == str2double(groupNames{g}(end));
        uSubTb = uSubTb(m,:);
        
        % Sort units by modulation index
        uSubTb = sortrows(uSubTb, 'mi', 'descend');

        % Sort units by depth within each page
        upp = 15;
        nPages = ceil(height(uSubTb) / upp);
        pInd = repelem(1:nPages, upp);
        uSubTb.pInd = pInd(1:height(uSubTb))';
        uSubTb = sortrows(uSubTb, {'pInd', 'depth'}, 'ascend');

        % Make plot
        uIdPage = NaN(upp, nPages);
        uIdPage(1:height(uSubTb)) = uSubTb.clusId;
        for p = 1 : min(nPages, 20)
            f = MPlot.Figure(22300); clf
            LMV.Plot.Overview(senTbs{i}, ceCell{i}, 'TaskPhase', 'full', ...
                'UnitIds', uIdPage(:,p), 'Page', p, 'Folder', fullfile(anaDir, 'engage', groupNames{g}));
        end
    end
end

%%

uId = 410100138; % delay NLR 3*, Zeta 0
uId = 410100451; % 

pn = 'prod';

uIdx = NP.Unit.ClusId2Ind(uId, uTb);
tSpk = zTb.spikeTimes{uIdx};
phaseWin = zTb.(pn+"Win"){uIdx,4};
% [pVal, sZeta, sRate, lat] = zetatest(tSpk, phaseWin(:,1)-1, median(diff(phaseWin,1,2))+2, 500, 3);
[pVal, sZeta, sRate, lat] = zetatest(tSpk, phaseWin, median(diff(phaseWin,1,2)), 500, 3);

fprintf("Unit %i: ZETA pval = %f; NLR sig = %i\n", uId, pVal, rTb.(pn)(uIdx));


