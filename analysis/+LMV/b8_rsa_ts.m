%% 

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'rsa');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end

%% Load data

% Sentence averaged (M1) data
ceSearch = MBrowse.Dir2Table(fullfile(NP.Data.GetAnalysisRoot, 'data', 'ce_m1_sentence-avg', '*_ce_m1.mat'));
ceArray = NP.CE.LoadSession(fullfile(ceSearch.folder, ceSearch.name));

for i = 1 : numel(ceArray)
    ce = ceArray(i);
    
    % Remove sentences that have too few repeats
    nRep = ce.GetColumn('taskValue', 'numTrial');
    ce.RemoveEpochs(nRep < 3);
    
%     % Remove sentences outside of the selected
%     tv = ce.GetTable('taskValue');
%     isSen = ismember(tv.stimId, sList);
%     ce.RemoveEpochs(~isSen);
%     
%     % Sort sentences to standard order
%     tv = ce.GetTable('taskValue');
%     I = zeros(size(sList));
%     for i = 1 : numel(sList)
%         I(i) = find(tv.stimId == sList(i), 1);
%     end
%     ce.SortEpochs(I);
end

% Load unit responsiveness table
sTest = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'ttest', "computed_ttest_clusTb.mat"));

%% Compute with a range of time shifts in spike rates

maskNames = {'stim', 'prod'};

dtLists = struct;
dtLists.stim = -0.4 : 0.02 : 0.1;
dtLists.prod = -0.1 : 0.02 : 0.4;

nRec = numel(ceArray);
nMask = numel(maskNames);
nDT = max(structfun(@numel, dtLists));
sArray = cell(nRec, nMask, nDT);

for i = 1 : nRec
    ce = ceArray(i);
    fprintf("\nRun RSA for %s (%i/%i)\n", NP.SE.GetID(ce), i, nRec);
    
    for j = 1 : nMask
        mn = maskNames{j};
        dt = dtLists.(mn);
        fprintf("%s task phase\n", mn);
        
        % Select speech responsive units
        isRec = sTest.clusTb.recId == NP.SE.GetID(ce);
        pval = sTest.clusTb{isRec,mn};
        uMask = any(pval*size(pval,2) < 0.05, 2);
        
%         % Select all units
%         uMask = true(ce.numResp, 1);
        
        parfor k = 1 : numel(dt)
            fprintf("dt = %g (%i/%i)\n", dt(k), k, numel(dt));
            
            % Configure analysis
            ops = ce.userData.rsOps;
            
            featFam = struct;
            featFam.Phone = [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.consonants];
            featFam.MannerPlace = {'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', 'nasals', 'approximants', 'labial', 'velar', 'coronal', 'glottal', 'dental'};
            featFam.Artic = ops.featVars.artic;
            featFam.Acous = [ops.featVars.inten ops.featVars.pitch];
            switch mn
                case 'stim'
                    featFam.Mel = {'speaker1'};
                case 'prod'
                    featFam.Mel = {'mic'};
            end
            
            ops.featFam = featFam;
            ops.unitInd = find(uMask);
            ops.maskName = mn;
            ops.dtNeural = dt(k);
            ops.downsample = 2;
            ops.nBoot = 0; % not computing pval
            
            % Run RSA
            sRez = NP.RSA.Run(ce, ops);
            sArray{i,j,k} = NP.RSA.LightenResult(sRez);
        end
    end
end

cachePath = fullfile(anaDir, "computed_temporal_psm.mat");
save(cachePath, 'ceSearch', 'maskNames', 'dtLists', 'sArray');

%% Plot similarity across time shifts

recTb = table;
recTb.recId = arrayfun(@(x) string(NP.SE.GetID(x)), ceArray);
recTb.region = arrayfun(@(x) string(NP.SE.GetRegion(x)), ceArray);
regions = unique(recTb.region);

for r = 1 : numel(regions)
    % Get region name
    rn = regions(r);
    rm = rn == recTb.region;
    
    % Set up figure
    f = MPlot.Figure(100+r); clf
    f.Name = rn;
    NP.RSA.PlotSimTimecourse(sArray(rm,:,:), recTb.region(rm));
    MPlot.Paperize(f, 'ColumnsWide', 1.8, 'ColumnsHigh', 0.5);
    figPath = fullfile(anaDir, sprintf("rs_by_time_offsets_%s.png", rn));
    exportgraphics(f, figPath);
end

%% Plot similarity at optimal time shifts

f = MPlot.Figure(200); clf
tl = tiledlayout('flow');
regions = ["mPrCG", "vPrCG", "IFG", "STG"];
NP.RSA.PlotSimOptimal(sArray(:,1,:), recTb, regions);
NP.RSA.PlotSimOptimal(sArray(:,2,:), recTb, regions);
MPlot.Paperize(f, 'ColumnsWide', 2, 'ColumnsHigh', 0.7);
exportgraphics(f, fullfile(anaDir, "rs_max.png"));

%% Standalone RSA

mn = 'prod';
ce = ceArray(4);
disp(NP.SE.GetID(ce));

% Select responsive units
isRec = sTest.clusTb.recId == NP.SE.GetID(ce);
pval = sTest.clusTb{isRec,mn};
uMask = any(pval < 0.05, 2);

% Configure analysis
ops = ce.userData.rsOps;

featFam = struct;
% featFam.Phone = [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.consonants];
featFam.Phone = unique([NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.approximants, NP.Phone.velar, NP.Phone.postalveolar, NP.Phone.alveolar, NP.Phone.dental, NP.Phone.labial], 'stable');
featFam.MannerPlace = {'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', 'nasals', 'approximants', 'labial', 'velar', 'coronal', 'dental'};
featFam.Artic = ops.featVars.artic;
featFam.Acous = [ops.featVars.inten ops.featVars.pitch];
switch mn
    case 'stim'
        featFam.Mel = {'speaker1'};
    case 'prod'
        featFam.Mel = {'mic'};
end

ops.featFam = featFam;
ops.unitInd = find(uMask);
ops.maskName = mn;
ops.dtNeural = 0.17;
ops.downsample = 2;
ops.nBoot = 0; % not computing pval

sRez = NP.RSA.Run(ce, ops);

%% Check input

f = MPlot.Figure(575); clf
tl = tiledlayout('flow');
tl.Padding = 'compact';
ax = nexttile;
NP.RSA.PlotInputFeatures(sRez);
ax.XLim = [0.01 1.57];
MPlot.Paperize(f, 'ColumnsWide', 1, 'ColumnsHigh', 1.3, 'FontSize', 4);
exportgraphics(f, fullfile(anaDir, sprintf("example_feat_%s_%s.png", NP.SE.GetID(ce), mn)));

%% Plot RDMs with sentence grouping

RDM = sRez.RDM;
repNames = sRez.repNames;

% Use sentences for grouping
cond = "s" + sRez.stimIdx;

% Plot RDMs
f = MPlot.Figure(11); clf
tl = tiledlayout('flow');
tl.Padding = 'compact';
for i = 1 : numel(RDM)
    ax = nexttile;
    NP.RSA.PlotRDM(RDM{i}, cond);
    ax.Title.String = repNames(i);
%     colorbar;
end
MPlot.Paperize(f, 'ColumnsWide', 1.4, 'Aspect', 0.667);
% exportgraphics(f, fullfile(anaDir, sprintf("rdms_%s_%s_group-sent.png", NP.SE.GetID(ce), mn)));

%% Plot RDMs with phoneme grouping

RDM = sRez.RDM;
repNames = sRez.repNames;

% Use phonemes for grouping
% evts = {'high', 'mid', 'low', 'plosives', 'fricatives', 'nasals', 'approximants'};
% evts = featFam.Phone;
F = sRez.input{2};
[~, iMaxF] = max(F, [], 2);
evts = string(MLing.ARPA2IPA(featFam.Phone));
cond = categorical(evts(iMaxF)', evts, 'Ordinal', true);

% Plot RDMs
f = MPlot.Figure(11); clf
tl = tiledlayout('flow');
tl.Padding = 'compact';
for i = 1 : numel(RDM)
    ax = nexttile;
    NP.RSA.PlotRDM(RDM{i}, cond);
    ax.Title.String = repNames(i);
%     colorbar;
end
MPlot.Paperize(f, 'ColumnsWide', 1.4, 'Aspect', 0.667);
% exportgraphics(f, fullfile(anaDir, sprintf("rdms_%s_%s_group-ph.png", NP.SE.GetID(ce), mn)));

%% 

f = MPlot.Figure(2); clf
NP.RSA.PlotPSM(sRez);
h = gca;
h.Title = strrep(NP.SE.GetID(ce), '_', '-') + " " + mn;
% h.Interpreter = 'none'; % this option is only available in 2023b
MPlot.Paperize(f, 'ColumnsWide', .6, 'Aspect', .85);
exportgraphics(f, fullfile(anaDir, sprintf("psm_%s_%s.png", NP.SE.GetID(ce), mn)));

return
%% Output data for Python

% save(fullfile(anaDir, NP.SE.GetID(ce)+".mat"), '-struct', 'sRSA');

