%% 

anaDir = LMV.Data.GetAnalysisDir('feat_qc');

%% Load TIMIT corpus data

% Load TIMIT corpus se
s = load(fullfile(anaDir, 'timit_corpus_se.mat'));
seCorp = s.se.Duplicate;

% % Choose AKT version
% seCorp.SetTable('artic', seCorp.GetTable('artic1'), 'timeSeries', seCorp.GetReferenceTime('artic1'));

% Vectorize se
seCorp = seCorp.Duplicate({'taskTime', 'artic1', 'artic2'}, false); % make it lightweight
seCorp.SliceSession(0, 'absolute');

%% Load LMV data

% Load recording
srcTb = LMV.Data.FindSource('phone_stats');
sePaths = [srcTb.path; fullfile(NP.Data.GetProjectRoot, "analysis", "feat_qc", "preproc_output", "NP41_B1_se_new-aai.mat")];
seLMV = NP.SE.LoadSession(sePaths);

% Enrich se
ops = NP.Param.Enrich;
ops.isMel = true;
ops.isTimitGT = true;
seLMV(1) = LMV.SE.Transform(seLMV(1), 'enrich', ops);

% Vectorize se
for i = 1 : numel(seLMV)
    seLMV(i) = seLMV(i).Duplicate({'taskTime', 'artic'}, false); % make it lightweight
    seLMV(i).SliceSession(0, 'absolute');
end

% Concatenate se's
% seLMV = Merge(seLMV);

%% Make phoneme tables

poolList = {'stimGT', 'stim', 'prod'};
poolList = {'stim', 'prod'};

phTbs = cell(numel(poolList)+1, 2);
dVowels = phTbs;
dCons = phTbs;

for k = 1 : 2
    seCorp.SetTable('artic', seCorp.GetTable("artic"+k), 'timeSeries', 0);
    phTbs{1,k} = NP.Phone.PoolPhonemes(seCorp, 'stimGT');
    
    for j = 1 : numel(poolList)
        phTbs{j+1,k} = NP.Phone.PoolPhonemes(seLMV(k), poolList{j});
    end
    
    % Fit LDA models in vowel space
    for j = 1 : size(phTbs,1)
        dVowels{j,k} = NP.Phone.ComputeLDA(phTbs{j,k}, 'vowel');
        dCons{j,k} = NP.Phone.ComputeLDA(phTbs{j,k}, 'consonant');
    end
end

%% Plot LDA vowel space

f = MPlot.Figure(230); clf
tl = tiledlayout('flow');
tl.Padding = 'compact';

titles = ["Corpus audio, GT labels", "LMV stim audio, MFA labels", "LMV prod audio, MFA labels"];

for k = 1 : 2
    for j = 1 : size(phTbs,1)
        ax = nexttile;
        NP.Phone.PlotLDA(ax, phTbs{j,k}, dVowels{1,k}); % always embed in corpus space
        ax.Title.String = titles{j};
    end
end

MPlot.Paperize(f, 'ColumnsWide', 1.6, 'Aspect', .6);
exportgraphics(f, fullfile(anaDir, "phone_in_aai_vowel_space.png"));

%% Plot LDA consonant space

f = MPlot.Figure(231); clf
tl = tiledlayout('flow');
tl.Padding = 'compact';

titles = ["Corpus audio, GT labels", "LMV stim audio, MFA labels", "LMV prod audio, MFA labels"];

for k = 1 : 2
    for j = 1 : size(phTbs,1)
        ax = nexttile;
        NP.Phone.PlotLDA(ax, phTbs{j,k}, dCons{1,k}); % always embed in corpus space
        ax.Title.String = titles{j};
    end
end

MPlot.Paperize(f, 'ColumnsWide', 1.6, 'Aspect', .6);
exportgraphics(f, fullfile(anaDir, "phone_in_aai_consonant_space.png"));


return
%% Plot LDA vowel space

f = MPlot.Figure(230); clf
tl = tiledlayout('flow');
tl.Padding = 'compact';

% Corpus
ax0 = nexttile;
NP.Phone.PlotLDA(ax0, phCorp, dCorp);
ax0.Title.String = 'TIMIT corpus (GT)';
% ax0.XLim = [-.6 .7];
% ax0.YLim = [-.7 .6];
% ax0.XTick = [];
% ax0.YTick = [];
% h = legend(ax0);
% h.Visible = 'on';

% LMV stim (GT)
ax = nexttile;
NP.Phone.PlotLDA(ax, phStimGT, dStim);
ax.Title.String = 'LMV stim (GT)';
% h = legend(ax);
% h.Visible = 'on';

% LMV stim (MFA)
ax = nexttile;
NP.Phone.PlotLDA(ax, phStim, dStim);
ax.Title.String = 'LMV stim (MFA)';
% h = legend(ax);
% h.Visible = 'on';

% LMV prod
ax = nexttile;
NP.Phone.PlotLDA(ax, phProd, dProd);
ax.Title.String = 'LMV production (MFA)';
% h = legend(ax);
% h.Visible = 'on';

MPlot.Paperize(f, 'ColumnsWide', 2.2, 'Aspect', .22);
% exportgraphics(f, fullfile(anaDir, "phone_in_aai-"+k+"_vowel_space.png"));


%% UMAP embedding of artic features

% phList = {'AA', 'AE', 'AH', 'UH', 'UW', 'IH', 'IY'};
% phSelect = setdiff(NP.Phone.vowels, NP.Phone.diphthongs);
% isSelect = ismember(phe.GetParentLabel, phList);
% phe = phe(isSelect).PickupTimeseries('artic', artic.time{1}, cell2mat(artic{:,2:end}));

% Select phonemes
isSelect = ismember(phCorp.arp, phList);
phCSub = phCorp(isSelect, :);
phCSub = sortrows(phCSub, 'arp');

X = double(phCSub.artic(:,1:13));

% Compute embedding
rng(61);
Z = run_umap(X, 'metric', 'correlation', 'verbose', 'text', ...
    'min_dist', 0.2, 'n_neighbors', 15, 'randomize', false);

%% Plot UMAP

f = MPlot.Figure(231); clf
gscatter(Z(:,1), Z(:,2), categorical(phCSub.arp));
h = legend();
h.String = MLing.ARPA2IPA(phList);
ax = gca;
ax.XLabel.String = "UMAP 1";
ax.YLabel.String = "UMAP 2";
ax.XTick = [];
ax.YTick = [];
ax.Box = 'off';
MPlot.Paperize(f, 'ColumnsWide', 0.8, 'ColumnsHigh', 0.6);


%% 

% Options for resampling
ops.rsWin = [tt.stimOn tt.stimOff] + [-1 1]*0.5;
ops = NP.Param.Resample;
ops.rsBinSize = 0.01;
ops.rsArgs = {};
ops.featVars.mel = {'speaker1'};
ops.featVars.artic = {'la', 'pro', 'ttcl', 'tbcl', 'ja', 'ttcc', 'tbcc', 'tdcc'};

fTb = NP.SE.ResampleFeatures(seSen, ops);


