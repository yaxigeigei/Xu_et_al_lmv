%% 

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'rsa');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end
% srcTb = NP.Data.FindSource('LMV');

%% Load data

% Load unit responsiveness table
sTest = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'ttest', "computed_ttest_clusTb.mat"));

% Load event ce
ceSearch = MBrowse.Dir2Table(fullfile(NP.Data.GetAnalysisRoot, 'rsa', 'ce', '*_ce.mat'));
ceArray = NP.CE.LoadSession(fullfile(ceSearch.folder, ceSearch.name));

%% 

% Configure analysis
ops = ce.userData.ops;

repFeats = struct;
repFeats.Artic = {'ja', 'la', 'pro', 'ttcc', 'tbcc', 'tdcc', 'ttcl', 'tbcl', 'v_x', 'v_y'};
repFeats.Acous = [{'env', 'peakEnv', 'peakRate'}, {'voicing', 'phrase', 'accent', 'drF0'}];
repFeats.MannerPlace = {'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', 'nasals', 'approximants', 'labial', 'velar', 'coronal', 'dental'};
ops.repFeats = repFeats;

ops.nBoot = 0; % not computing pval

targets = ["stim", "prod"];
conds = {{'P', 'B', 'M', 'F', 'V' ,'TH' ,'DH', 'T', 'D', 'N', 'S', 'Z', 'CH', 'JH', 'SH', 'ZH', 'NG', 'K', 'G', 'Y', 'L', 'R', 'W'}, ...
    {'OY', 'OW', 'AO', 'AA', 'AW', 'AY', 'AE', 'EH', 'EY', 'IY', 'IH', 'AH', 'UW', 'ER', 'UH'}};

for r = 1 : numel(ceArray)
    % Run RSA
    ce = ceArray(r);
    disp(NP.SE.GetID(ce));
    ceSub = ce.Split(ones(ce.numEpochs,1));
    
    ssRez = cell(numel(targets), numel(conds));
    for i = 1 : numel(targets)
        for j = 1 : numel(conds)
            % Select responsive units
            isRec = sTest.clusTb.recId == NP.SE.GetID(ce);
            pval = sTest.clusTb{isRec, targets(i)};
            uMask = any(pval < 0.05, 2);
            ops.unitInd = find(uMask);
            
            % Select phonemes
            ops.eventMask = conds{j};
            
            % Compute
            sRez = NP.RSA.PrepareEventInputs(ceSub(i), ops);
            sRez = NP.RSA.ComputeRDMs(sRez);
%             sRez = NP.RSA.ComputeSim(sRez);
            ssRez{i,j} = sRez;
        end
    end
    
    % Plot RDMs
    f = MPlot.Figure(12); clf
    tl = tiledlayout('flow');
    tl.Padding = 'compact';
    for i = 1 : numel(targets)
        for j = 1 : numel(conds)
            sRez = ssRez{i,j};
            RDM = sRez.RDM;
            repNames = sRez.repNames;
            cond = categorical(MLing.ARPA2IPA(sRez.eventNames), MLing.ARPA2IPA(conds{j}), 'Ordinal', true);
            for k = 1 : numel(RDM)
                ax = nexttile;
                NP.RSA.PlotRDM(RDM{k}, cond);
                ax.Title.String = sprintf("%s (%s, %s)", repNames(k), NP.SE.GetID(sRez.ce), targets(i));
                % colorbar;
            end
        end
    end
    % MPlot.Paperize(f, 'ColumnsWide', 1.4, 'Aspect', 0.667);
%     exportgraphics(f, fullfile(anaDir, sprintf("rdms_by-phone_%s_%s.png", NP.SE.GetRegion(ce), NP.SE.GetID(ce))));
end



%% Standalone RSA

ce = ceArray(1);
disp(NP.SE.GetID(ce));
ce

% Configure analysis
ops = ce.userData.ops;

repFeats = struct;
repFeats.Artic = {'ja', 'la', 'pro', 'ttcc', 'tbcc', 'tdcc', 'ttcl', 'tbcl', 'v_x', 'v_y'};
repFeats.Acous = [{'env', 'peakEnv', 'peakRate'}, {'voicing', 'phrase', 'accent', 'drF0'}];
repFeats.MannerPlace = {'high', 'mid', 'low', 'front', 'back', 'rounded', 'plosives', 'fricatives', 'nasals', 'approximants', 'labial', 'velar', 'coronal', 'dental'};
ops.repFeats = repFeats;

% ph2plot = [NP.Phone.approximants, NP.Phone.velar, NP.Phone.postalveolar, NP.Phone.alveolar, NP.Phone.dental, NP.Phone.labial];
% ph2plot = unique(ph2plot, 'stable');
% ph2plot = flip(ph2plot);
% ph2plot = {'P', 'B', 'M', 'F', 'V' ,'TH' ,'DH', 'T', 'D', 'N', 'S', 'Z', 'CH', 'JH', 'SH', 'ZH', 'NG', 'K', 'G', 'Y', 'L', 'R', 'W'};
% ph2plot = {'OY', 'OW', 'AO', 'AA', 'AW', 'AY', 'AE', 'EH', 'EY', 'IY', 'IH', 'AH', 'UW', 'ER', 'UH'};

ops.nBoot = 0; % not computing pval

% Run RSA
targets = ["stim", "prod"];
conds = { ...
    {'P', 'B', 'M', 'F', 'V' ,'TH' ,'DH', 'T', 'D', 'N', 'S', 'Z', 'CH', 'JH', 'SH', 'ZH', 'NG', 'K', 'G', 'Y', 'L', 'R', 'W'}, ...
    {'OY', 'OW', 'AO', 'AA', 'AW', 'AY', 'AE', 'EH', 'EY', 'IY', 'IH', 'AH', 'UW', 'ER', 'UH'}};

ceSub = ce.Split(ones(ce.numEpochs,1));

ssRez = cell(numel(targets), numel(conds));
for i = 1 : numel(targets)
    for j = 1 : numel(conds)
        % Select responsive units
        isRec = sTest.clusTb.recId == NP.SE.GetID(ce);
        pval = sTest.clusTb{isRec, targets(i)};
        uMask = any(pval < 0.05, 2);
        ops.unitInd = find(uMask);
        
        % 
        ops.eventMask = conds{j};
        
        % 
        sRez = NP.RSA.PrepareEventInputs(ceSub(i), ops);
        sRez = NP.RSA.ComputeRDMs(sRez);
        sRez = NP.RSA.ComputeSim(sRez);
        ssRez{i,j} = sRez;
    end
end

%% Plot RDMs with phoneme grouping

% Use phonemes for grouping
% ph2plot = [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.approximants, NP.Phone.velar, NP.Phone.postalveolar, NP.Phone.alveolar, NP.Phone.dental, NP.Phone.labial];
% ph2plot = [NP.Phone.approximants, NP.Phone.velar, NP.Phone.postalveolar, NP.Phone.alveolar, NP.Phone.dental, NP.Phone.labial];
% ph2plot = unique(ph2plot, 'stable');

% Plot RDMs
f = MPlot.Figure(12); clf
tl = tiledlayout('flow');
tl.Padding = 'compact';
for i = 1 : numel(targets)
    for j = 1 : numel(conds)
        sRez = ssRez{i,j};
        RDM = sRez.RDM;
        repNames = sRez.repNames;
        cond = categorical(MLing.ARPA2IPA(sRez.eventNames), MLing.ARPA2IPA(conds{j}), 'Ordinal', true);
        for k = 1 : numel(RDM)
            ax = nexttile;
            NP.RSA.PlotRDM(RDM{k}, cond);
            ax.Title.String = sprintf("%s (%s, %s)", repNames(k), NP.SE.GetID(sRez.ce), targets(i));
%             colorbar;
        end
    end
end
% MPlot.Paperize(f, 'ColumnsWide', 1.4, 'Aspect', 0.667);
% exportgraphics(f, fullfile(anaDir, sprintf("rdms_%s_group-ph.png", NP.SE.GetID(ce))));

%% 

f = MPlot.Figure(2); clf
NP.RSA.PlotPSM(sRez);
h = gca;
h.Title = strrep(NP.SE.GetID(ce), '_', '-') + " " + mn;
% h.Interpreter = 'none'; % this option is only available in 2023b
MPlot.Paperize(f, 'ColumnsWide', .6, 'Aspect', .85);
% exportgraphics(f, fullfile(anaDir, sprintf("psm_%s_%s.png", NP.SE.GetID(ce), mn)));


return
