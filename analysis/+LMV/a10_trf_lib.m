%% 

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'trf');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end
srcTb = NP.Data.FindSource('lmv');

% Load PETH data for modulation indices
sNMF = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'nmf', 'computed_nmfc_all.mat'));

%% Plot TRFs for all units

sets = NP.TRF.featSets;
sets = "strf";
targets = NP.TRF.targets;

for f = 1 : numel(sets)
    % Load model tables
    mdlNames = sets(f) + "_" + targets;
    mdlDir = fullfile(anaDir, sets(f), 'mdls');
    uTbs = NP.TRF.LoadModels(mdlNames, mdlDir, srcTb.recId);
    
    % Add modulation indices
    for i = 1 : numel(uTbs)
        uTb = uTbs{i};
        for u = 1 : height(uTb)
            k = find(sNMF.ce.clusTb.clusId == uTb.clusId(u), 1);
            uTb.mi(u) = sNMF.ce.clusTb.mi(k);
        end
        uTbs{i} = uTb;
    end
    
    % Plot TRFs for all units
    for i = 1 : numel(uTbs)
        uTb = uTbs{i};
        
        % Select stim or prod tuned units
        isSpeech = ~cellfun(@isempty, uTb.(mdlNames(1)));
        uTb = uTb(isSpeech,:);
        
        % Sort units by modulation index
        uTb = sortrows(uTb, 'mi', 'descend');
        
        % Sort units by depth within each page
        upp = 20;
        nPages = ceil(height(uTb) / upp);
        pInd = repelem(1:nPages, upp);
        uTb.pInd = pInd(1:height(uTb))';
        uTb = sortrows(uTb, {'pInd', 'depth'}, 'ascend');
        
        % Plot units
        for p = 1 : min(nPages, 3)
            MPlot.Figure(4430); clf
            NP.TRF.PlotWeightsCombo(uTb, mdlNames, 'UnitsPerPage', upp, 'Page', p, ...
                'Folder', fullfile(anaDir, sets(f), 'weights'));
            
            MPlot.Figure(4430); clf
            NP.TRF.PlotWeightsCombo(uTb, mdlNames, 'UnitsPerPage', upp, 'Page', p, ...
                'PanelArgs', {'WeightAlpha', 0.05}, ...
                'Folder', fullfile(anaDir, sets(f), 'weights-p05'));
        end
    end
end

%% Load se

% Load se
seArray = NP.SE.LoadSession(srcTb.path);

% Compute sentence PETH
[~, senTbs] = arrayfun(@(x) LMV.SE.Transform(x), seArray, 'Uni', false);
[ceCell, senTbs] = cellfun(@(x) LMV.SE.ComputeSentencePETH(x), senTbs, 'Uni', false);

%% Plot rasters for all units

for i = 1 : numel(uTbs)
    uTb = uTbs{i};
    
    % Select stim or prod tuned units
    isSpeech = ~cellfun(@isempty, uTb.(targets(1)));
    uTb = uTb(isSpeech,:);
    
    % Sort units by modulation index
    uTb = sortrows(uTb, 'mi', 'descend');
    
    % Sort units by depth within each page
    upp = 3;
    nPages = ceil(height(uTb) / upp);
    pInd = repelem(1:nPages, upp);
    uTb.pInd = pInd(1:height(uTb))';
    uTb = sortrows(uTb, {'pInd', 'depth'}, 'ascend');
    
    %  
    uIdPage = NaN(upp, nPages);
    uIdPage(1:height(uTb)) = uTb.clusId;
    
    % 
    for t = 1 : numel(targets)
        for p = 1 : min(nPages, 20)
            MPlot.Figure(22300); clf
            LMV.Plot.Overview12(senTbs{i}, ceCell{i}, 'TaskPhase', targets{t}, ...
                'UnitIds', uIdPage(:,p), 'Page', p, 'Folder', fullfile(anaDir, ['raster_' targets{t}]));
        end
    end
end
