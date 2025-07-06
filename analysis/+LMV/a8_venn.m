%% Save responsiveness table as a parquet file for plotting in Python

method = 'signrank';
% method = 'ttest';
% method = 'zeta';
anaDir = LMV.Data.GetAnalysisDir('phase_resp', method, 'venn');

%% Load responsiveness test results

load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', method, "computed_test_clusTb.mat"), 'clusTb');

%% 

% Exclude inhibited units
phaseNames = intersect(["atten", "stim", "delay", "init", "prod", "iti"], clusTb.Properties.VariableNames);
P = clusTb{:,phaseNames};
if all(ismember(phaseNames+"Resp", clusTb.Properties.VariableNames))
    R = clusTb{:,phaseNames+"Resp"};
    P(R<0) = NaN;
end
clusTb{:,phaseNames} = P;

% Take minimal pvals of each phase
for i = 1 : numel(phaseNames)
    pn = phaseNames(i);
    clusTb.(pn) = min(clusTb.(pn), [], 2);
end

%% Venn diagrams of response types

% Export clusTb for python
colMask = ~varfun(@(x) iscell(x) | isstruct(x) | size(x,2)>1, clusTb, 'OutputFormat', 'uniform');
parquetwrite(fullfile(anaDir, "clusTb.parquet"), clusTb(:,colMask), ...
    'VariableCompression', 'uncompressed', 'VariableEncoding', 'plain');

% Run python scirpt
% \babble\python\venn\venn_task_phase_resp.py

return
%% Plot model weights of selected units for inspection

for i = 1 %: numel(clusTbs)
    % Find example units
    clusTb = clusTbs{i};
    recId = clusTb.recId(1);
    uId = NP.Unit.GetSelectedClusId(recId);
    uInd = find(ismember(clusTb.clusId, uId));
    
    % 
    f = MPlot.Figure(8772); clf
    
    nr = 5;
    nc = 6;
    nameGrid = repmat({'phase'}, nr, nc);
    uIndGrid = NaN(nc, nr);
    uIndGrid(1:numel(uInd)) = uInd;
    uIndGrid = uIndGrid';
    NP.TRF.PlotGrid(clusTb, uIndGrid, nameGrid);
    
    MPlot.Paperize(f, 'ColumnsWide', nc*0.3, 'ColumnsHigh', nr*0.25);
    figName = fullfile(anaDir, 'weights', sprintf("%s_phase_page%i.png", recId, 0));
%     exportgraphics(f, figName);
end

