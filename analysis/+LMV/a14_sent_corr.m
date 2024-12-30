%% Compute response correlations between different task phases

anaDir = LMV.Data.GetAnalysisDir("sent_resp", "corr");

srcTb = LMV.Data.FindSource([]);

%% Load test data

resp = LMV.Resp.LoadPhaseResponse(srcTb.recId);
rTest = LMV.Resp.LoadPhaseResponseTest(srcTb.recId);
sTest = LMV.Resp.LoadSentenceSelectTest(srcTb.recId);

%% Compute response correlation between task phases

phaseNames = ["stim", "delay", "init", "prod"];
stimIdList = LMV.Param.stimIdList12;

nRec = height(srcTb);
nPhase = numel(phaseNames);
nSen = numel(stimIdList);
clusTbs = cell(nRec,1);

for i = 1 : nRec
    fprintf("\nCompute phase response correlations for %s\n", srcTb.recId{i});
    
    % Unpack variables
    respTb = resp(i).respTb;
    dropoutTb = resp(i).dropoutTb;
    phaseTb = resp(i).phaseTb;
    stimIdTb = resp(i).stimIdTb;
    clusTb = resp(i).clusTb;
    rSigTb = rTest(i).sigTb(:, phaseNames);
    sSigTb = sTest(i).sigTb(:, phaseNames);
    
    % Set dropout periods to NaN
    respTb{:,:}(dropoutTb{:,:}) = NaN;
    
    % Exclude sentence non-selective units
    isEx = ~any(sSigTb{:,:}, 2);
    if all(isEx)
        continue
    end
    respTb(:,isEx) = [];
    clusTb(isEx,:) = [];
    rSigTb(isEx,:) = [];
    
    % Add significance to clusTb
    clusTb = [clusTb rSigTb];
    
    nUnit = width(respTb);
    nTrial = sum(phaseTb.(1));
    for u = 1 : nUnit
        % Collect trial-by-phase response matrix for each unit
        rPhases = zeros(nTrial, nPhase);
        for j = 1 : nPhase
            pn = phaseNames(j);
            isPhase = phaseTb.(pn);
            rPhases(:,j) = respTb{phaseTb.(pn), u};
        end
        
        % Set responses in non-responsive phases to NaN
        isResp = rSigTb{u,:};
        rPhases(:, ~isResp) = NaN;
        
        % Compute pairwise correlations and p-values
        [clusTb.phaseCorr{u}, clusTb.phaseCorrPval{u}] = corr(rPhases, 'Rows', 'pairwise');
%         [clusTb.phaseCorr{u}, clusTb.phaseCorrPval{u}] = LMV.Resp.BootCorr(rPhases, 'NBoot', 100);
    end
    
    clusTbs{i} = clusTb;
end

clusTbCat = cat(1, clusTbs{:});

%% Inspect session rasters with interactive scatter violin plots

f = MPlot.Figure(764282); clf
LMV.Resp.PlotPhaseCorrScatter(clusTbCat, phaseNames);
MPlot.Paperize(f, 'ColumnsWide', 2, 'ColumnsHigh', 0.4);
exportgraphics(f, fullfile(anaDir, "phase_corr_violin.png"));

%% Plot CDFs of task phase correlations

f = MPlot.Figure(764285); clf
LMV.Resp.PlotPhaseCorrCDF(clusTbCat, phaseNames);
MPlot.Paperize(f, 'ColumnsWide', 2, 'ColumnsHigh', 0.3);
exportgraphics(f, fullfile(anaDir, "phase_corr_cdf.png"));




