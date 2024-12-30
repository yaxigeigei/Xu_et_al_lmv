%% MNI brain surface recon with NP site overlay

% Add the iElvis code to path
addpath(genpath(NP.Recon.iElvisPath));

% Folder to save plots
anaDir = fullfile(NP.Data.GetAnalysisRoot, "recon");

% Source table
srcTb = LMV.Data.FindSource([]);

%% 

% Plot sites in colors
f = MPlot.Figure(700); clf
NP.Recon.PlotSitesOnAvgBrain(srcTb);
MPlot.Paperize(f, .5, .4);
exportgraphics(f, fullfile(anaDir, "recon_lmv_sites.png"), Resolution=300);
exportgraphics(f, fullfile(anaDir, "recon_lmv_sites.pdf"), Resolution=300);

% Plot sites in black
f = MPlot.Figure(700); clf
NP.Recon.PlotSitesOnAvgBrain(srcTb, 'Color', 'k');
MPlot.Paperize(f, .5, .4);
exportgraphics(f, fullfile(anaDir, "recon_lmv_sites_bk.png"), Resolution=300);
exportgraphics(f, fullfile(anaDir, "recon_lmv_sites_bk.pdf"), Resolution=300);

% Plot sites with colors and labels
f = MPlot.Figure(701); clf
NP.Recon.PlotSitesOnAvgBrain(srcTb, 'ShowNames', true);
MPlot.Paperize(f, 1, .8);
exportgraphics(f, fullfile(anaDir, "recon_lmv_sites_labeled.png"), Resolution=300);
exportgraphics(f, fullfile(anaDir, "recon_lmv_sites_labeled.pdf"), Resolution=300);

return
%% Save an empty recon for schematic drawing

f = MPlot.Figure(703); clf

% Configure plotting
cfg = [];
cfg.fsurfSubDir = fullfile(NP.Recon.iElvisPath, 'freesurfer_brains');

cfg.figId = gcf;
cfg.view = 'l';
cfg.title = "";

cfg.elecShape = 'sphere';
cfg.elecColorScale = [0 1 1];
cfg.showLabels = 'n';
cfg.elecUnits = [];
cfg.elecCbar = 'n';
cfg.elecSize = 2;

plotPialSurf("cvs_avg35_inMNI152", cfg);

% Make the brain brighter than default
lightLevel = 2; % originally 1
hl = light('Color',[1 1 1]*0.15*lightLevel);
lightangle(hl, -90, -20);
hl = light('Color',[1 1 1]*0.2*lightLevel);
lightangle(hl, -30, 70);
hl = light('Color',[1 1 1]*0.2*lightLevel);
lightangle(hl, -150, 70);

exportgraphics(f, fullfile(anaDir, "recon_empty.png"));
