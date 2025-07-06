%% Plot sentence PETH overlay for all linker units

mdlName = LMV.Linker.currentModel;
anaDir = LMV.Data.GetAnalysisDir("linker", mdlName);

%% Load data

% Load linker clusTb with M2 sentence PETHs
load(fullfile(anaDir, "linker_clusTb_peth.mat"), "clusTb");

% Load linked positions
load(fullfile(anaDir, "linked_positions.mat"), "posTb");

%% Plot sentence response overlay of all linker units

types = LMV.Linker.types;
libDir = LMV.Data.GetAnalysisDir("linker", mdlName, "lib_peth");

for i = 1 : numel(types)
    f = MPlot.Figure(660+i); clf
    f.WindowState = "maximized";
    tl = tiledlayout(4,7);
    tl.Padding = "tight";
    
    lkInd = find(clusTb.hcGroup==types(i));
    for j = 1 : min(numel(lkInd), 24)
        k = lkInd(j);
        isPos = posTb.clusId == clusTb.clusId(k);
        ax = nexttile;
        LMV.Fig.SentenceRespOverlay(ax, clusTb(k,:), posTb(isPos,:), "ShowNonLinked", true, "ShowText", true);
        ax.YTickLabel = [];
        ax.XLabel.String = [];
        ax.XTickLabel = [];
        % return
    end
    % return
    % MPlot.Paperize(f, 2, 1.2);
    exportgraphics(f, fullfile(libDir, sprintf("sent_peth_%s.png", types(i))));
end
