%% Create a library of plots for inspecting results

mdlName = LMV.Linker.currentModel;
srcTb = LMV.Data.FindSource([]);

%% Load cached data

cacheDir = LMV.Data.GetAnalysisDir("linker", mdlName, "computed_profile");

sPf = cell(size(srcTb.recId));
for i = 1 : height(srcTb)
    disp(srcTb.recId(i));
    cacheFile = fullfile(cacheDir, sprintf("%s.mat", srcTb.recId(i)));
    sPf{i} = load(cacheFile);
end

%% Find units

clusTb = LMV.Linker.LoadClusTb(mdlName);

types = LMV.Linker.types;
hcGroups = arrayfun(@(x) clusTb.clusId(clusTb.hcGroup==x)', types, 'Uni', false);
maGroups = arrayfun(@(x) LMV.Linker.GetSelectedClusId(x, srcTb.recId), types, 'Uni', false);
maGroups = cellfun(@(x,y) setdiff(x,y), maGroups, hcGroups, 'Uni', false);
% cidGroups = cellfun(@(x,y) [x y], hcGroups, maGroups, 'Uni', false);

%% Make a library of plots

f = MPlot.Figure(8740);
f.WindowState = "maximized";

for recIdx = 1 : height(sPf)
    % Get variables
    s = sPf{recIdx};
    recId = s.recId;
    ce = s.ce;
    triTb = s.triTb;
    
    % Update linking positions
    LMV.Linker.LM.PredictProdResp(ce);
    pkTb = LMV.Linker.LM.FindScorePeaks(ce);
    
    % Find units
    clusId = unique(pkTb.clusId);
    
    % Plot summary plots
    for i = 1 : numel(clusId)
        for k = 1 : numel(types)
            if ismember(clusId(i), hcGroups{k})
                figDir = LMV.Data.GetAnalysisDir("linker", mdlName, "lib_profile", types(k));
            elseif ismember(clusId(i), maGroups{k})
                figDir = LMV.Data.GetAnalysisDir("linker", mdlName, "lib_profile", types(k)+"_manual");
            else
                continue
            end
            clf(f);
            LMV.Linker.LM.PlotUnitProfile(ce, pkTb, triTb, clusId(i));
            % return
            exportgraphics(f, fullfile(figDir, "u"+clusId(i)+".png"));
        end
    end
end
