%% 

anaDir = LMV.Data.GetAnalysisDir('feat_qc');

%% Load TIMIT corpus data

% Load TIMIT corpus se
s = load(fullfile(anaDir, 'timit_corpus_se.mat'));
seCorp = s.se.Duplicate;

% Choose AKT version
seCorp.SetTable('artic', seCorp.GetTable('artic1'), 'timeSeries', seCorp.GetReferenceTime('artic1'));

% Vectorize se
seCorp = seCorp.Duplicate({'taskTime', 'artic'}, false); % make it lightweight
seCorp.SliceSession(0, 'absolute');

%% Load LMV data

% Load recording
srcTb = LMV.Data.FindSource('phone_stats');
seLMV = NP.SE.LoadSession(srcTb.path);

for i = 1 : numel(seLMV)
    % Enrich se
    ops = NP.Param.Enrich;
    ops.isMel = true;
    ops.isPitch = false;
    ops.isArtic = false;
    ops.isTimitGT = true;
    seLMV(i) = LMV.SE.Transform(seLMV(i), 'enrich', ops);
    
    % Vectorize se
    seLMV(i) = seLMV(i).Duplicate({'taskTime', 'artic'}, false); % make it lightweight
    seLMV(i).SliceSession(0, 'absolute');
end

% Concatenate se's
% seLMV = Merge(seLMV);

%% Make phoneme tables

phCorp = NP.Phone.PoolPhonemes(seCorp, 'stimGT');
phStimGT = NP.Phone.PoolPhonemes(seLMV, 'stimGT');
phStim = NP.Phone.PoolPhonemes(seLMV, 'stim');
phProd = NP.Phone.PoolPhonemes(seLMV, 'prod');

%% Compute phoneme probabilities

[N, C] = histcounts(phCorp.cat, 'Normalization', 'probability');
histTb = table;
histTb.C = categorical(C)';
histTb.NCorp = N';
histTb.isVowel = ismember(histTb.C, NP.Phone.vowels);
histTb.isNotDiph = ~ismember(histTb.C, NP.Phone.diphthongs);
histTb = sortrows(histTb, {'isVowel', 'isNotDiph', 'NCorp'}, 'descend');

[N, C] = histcounts(phStimGT.cat, 'Normalization', 'probability');
for i = 1 : numel(C)
    k = histTb.C == C{i};
    histTb.NStimGT(k) = N(i);
end

[N, C] = histcounts(phStim.cat, 'Normalization', 'probability');
for i = 1 : numel(C)
    k = histTb.C == C{i};
    histTb.NStim(k) = N(i);
end

[r, pval] = corr(histTb.NCorp, histTb.NStimGT);

%% Plot phoneme probabilities

f = MPlot.Figure(54590); clf
text(histTb.NCorp*100, histTb.NStimGT*100, MLing.ARPA2IPA(histTb.C), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
ax = gca;
ax.XLim = [-.5 10];
ax.YLim = ax.XLim;
ax.XLabel.String = 'TIMIT corpus (%)';
ax.YLabel.String = 'LMV stim (%)';
ax.Title.String = sprintf('r = %.2f, p = %.2e', r, pval);
ax.Box = 'off';
MPlot.Paperize(f, 'ColumnsWide', 0.5, 'ColumnsHigh', 0.4);
exportgraphics(f, fullfile(anaDir, 'phone_corr.png'));
print(f, fullfile(anaDir, 'phone_corr.pdf'), '-dpdf');

%% 

f = MPlot.Figure(54592); clf
bar((1:height(histTb))', [histTb.NCorp histTb.NStimGT histTb.NStim]*100, 1);
legend('TIMIT corpus (GT)', 'LMV stim (GT)', 'LMV stim (MFA)');
ax = gca;
ax.XTick = (1:height(histTb))';
ax.XTickLabel = MLing.ARPA2IPA(histTb.C);
ax.TickLength(1) = 0;
ax.YLabel.String = 'Probability (%)';
ax.YGrid = 'on';
MPlot.Axes(ax);
MPlot.Paperize(f, 'ColumnsWide', 2.5, 'ColumnsHigh', 0.25);
MPlot.SavePNG(f, fullfile(anaDir, 'phone_prob'));
