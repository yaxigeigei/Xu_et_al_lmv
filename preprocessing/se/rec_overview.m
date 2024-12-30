%% This is the standard script for overview of LMV task

% Load se
[se, sePath] = NP.SE.LoadSession([], 'UserFunc', @(x) x.RemoveTable('LFP', 'niTime'));
if isempty(se)
    return
end

% Merge ks metadata (for old se)
NP.Unit.MergeMeta(se);

% Remove tiny clusters (for auto clustering)
clusTb = NP.Unit.GetClusTb(se);
NP.Unit.RemoveUnits(se, clusTb.numSpikes < 100);

% Create figure folder
seDir = fileparts(string(sePath));
recId = string(NP.SE.GetID(se));
ovDir = fullfile(seDir, recId);
if ~exist(ovDir, 'dir')
    mkdir(ovDir);
end

% Source table
srcTb = LMV.Data.FindSource([]);
srcTb = srcTb(srcTb.recId==recId,:);

%% Recon of recording site

addpath(genpath(NP.Recon.iElvisPath));

f = MPlot.Figure(789); clf
NP.Recon.PlotSitesOnAvgBrain(srcTb);
[~, region] = NP.SE.GetRegion(se);
title(sprintf("%s, %s", recId, region), 'Interpreter', 'none');
MPlot.Paperize(f, 'ColumnsWide', 1, 'ColumnsHigh', 0.8);
exportgraphics(f, fullfile(ovDir, recId+"_site_recon.png"));

%% Generate unit quality control plots

% Compute unit quality table
clusTb = NP.Unit.ComputeQualityTable(se);
save(fullfile(ovDir, recId+"_unit_qc_clusTb.mat"), 'clusTb');

% Convert to short IDs for plotting
clusTb.clusId = mod(clusTb.clusId, 1e5);
se.userData.ksMeta.clusTb = clusTb;

% Plot span
spanFigPath = fullfile(ovDir, recId+"_cluster_span");
NP.UnitPlot.Span(se, [], spanFigPath);

% Plot waveform
wfFigBasePath = fullfile(ovDir, recId+"_cluster_waveform");
NP.UnitPlot.WaveformArray(clusTb, wfFigBasePath);

% Plot ISI
isiFigBasePath = fullfile(ovDir, recId+"_cluster_isi");
NP.UnitPlot.ISIArray(clusTb, isiFigBasePath);
