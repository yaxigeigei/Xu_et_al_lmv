%% Create a library of plots for inspecting conv xc results

mdlName = "smooth_lm";
anaDir = LMV.Data.GetAnalysisDir("linker", mdlName);
figDir = LMV.Data.GetAnalysisDir("linker", mdlName, "lib_raster");
srcTb = LMV.Data.FindSource([]);

%% Find units

clusTb = LMV.Linker.LoadClusTb(mdlName);

types = LMV.Linker.types;
hcGroups = arrayfun(@(x) clusTb.clusId(clusTb.hcGroup==x), types, 'Uni', false);
maGroups = arrayfun(@(x) LMV.Linker.GetSelectedClusId(x, srcTb.recId)', types, 'Uni', false);
maGroups = cellfun(@(x,y) setdiff(x,y), maGroups, hcGroups, 'Uni', false);

% disp(hcGroups);
% disp(maGroups);

%% Make a library of plots

f = MPlot.Figure(8740);
f.WindowState = "maximized";

for k = 1 : numel(types)
    % Find units
    clusId = [NP.UnitPlot.MakeArrayID(hcGroups{k}, 10), NP.UnitPlot.MakeArrayID(maGroups{k}, 10)];
    
    % Plot summary plots
    for p = 1 : size(clusId,2)
        clf(f);
        LMV.Overview.SessionFromCache(clusId(:,p), 'DataSource', 'm2', 'StimIdList', LMV.Param.stimIdList14);
        figName = sprintf("%s_page%i.png", types(k), p);
        exportgraphics(f, fullfile(figDir, figName));
    end
end
