%% Create a library of plots for each unit

anaDir = LMV.Data.GetAnalysisDir('coding', 'stRF');
srcTb = LMV.Data.FindSource([]);

%% 

% Load segmental units
load(fullfile(anaDir, "segmental_clusTb.mat"), "clusTb");

% Load cached sequence data
cacheDir = LMV.Data.GetAnalysisDir("peri_event", "extracted_seq");
s = cell(height(srcTb), 1);
for i = 1 : numel(s)
    cacheFile = fullfile(cacheDir, srcTb.recId(i)+"_seqData.mat");
    if ~isfile(cacheFile)
        continue
    end
    s{i} = load(cacheFile, "se", "feats", "seqData");
    s{i}.seqData(2:end,:) = [];
end

% Get spectral frequencies
fMel = s{1}.se.userData.melMeta.F;

%% Find valid models

sets = LMV.TRF.featSets;
% targets = LMV.TRF.targets;
targets = ["stim", "prod"];

nTargets = numel(targets);
nSets = numel(sets);
mdlTbs = cell(nSets, nTargets);

for i = 1 : nSets
    for j = 1 : nTargets
        % Find valid models
        tb = clusTb;
        mn = sets(i)+"_"+targets(j);
        isMdl = ~cellfun(@isempty, tb.(mn));
        mdls = tb.(mn)(isMdl);
        
        if sets(i) == "strf"
            for k = 1 : numel(mdls)
                mdls{k}.feats = fMel;
            end
        end
        
        tb = tb(isMdl, 1:end-9);
        tb.mdls = mdls;
        
        mdlTbs{i,j} = tb;
    end
end

%% Extract peri-event responses and features from se

for i = 1 : nSets
    for j = 1 : nTargets
        tb = mdlTbs{i,j};
        if ismember(tb.mdls{1}.name, ["strf_prod", "artic3_stim"])
            continue
        end
        
        for u = 1 : height(tb)
            % Find the source recording
            M = tb.mdls{u};
            uid = str2double(erase(M.resps, "u"));
            isRec = strcmp(NP.SE.GetID(uid), srcTb.recId);
            assert(any(isRec), "No se matches the unit's source recording");
            se = s{isRec}.se;
            phSeqCell = s{isRec}.seqData.(targets(j)){"phone"}.seq1;
            
            % Extract data
            seqTb = LMV.RF.MakeSeqTb(M, phSeqCell, se);
            
            % Cache to tb
            tb.seqTb{u} = seqTb;
        end
        mdlTbs{i,j} = tb;
        % return
    end
end

%% Plot optimal weights of example units

f = MPlot.Figure(711);

for i = 3 %1 : numel(sets)
    % Create figure folder
    libDir = LMV.Data.GetAnalysisDir('coding', 'stRF', sets(i), 'lib');
    
    for j = 1 : numel(targets)
        tb = mdlTbs{i,j};
        if ismember(tb.mdls{1}.name, ["strf_prod", "artic3_stim"])
            continue
        end
        
        for k = 1 : height(tb)
            % Check computed
            M = tb.mdls{k};
            figName = sprintf("%s_%s_%s.png", sets(i), targets(j), M.resps);
            figPath = fullfile(libDir, figName);
            % if isfile(figPath)
            %     continue
            % end
            
            % Set up figure
            clf(f);
            nRows = 3;
            nCols = 3;
            tl = tiledlayout(nRows, nCols);
            tl.Padding = "compact";
            
            % Plot weight matrix
            ntArgs = MPlot.FindTileInd(nRows, nCols, 1:3, 1);
            ax = nexttile(ntArgs{:});
            LMV.TRF.PlotWeights2(M, 'ClusInfo', tb(k,:), 'Parent', ax);
            MPlot.Axes(ax);
            
            % Plot optimal weight vector
            ntArgs = MPlot.FindTileInd(nRows, nCols, 1:3, 2);
            ax = nexttile(ntArgs{:});
            m = M.mdls{M.r2Idx};
            m.feats = M.feats;
            % m.Beta = zscore(m.Beta);
            LMV.RF.PlotWeights(m, 'ClusInfo', tb(k,:), 'Orientation', 'vertical', 'Parent', ax);
            switch sets(i)
                case "phone"
                    ax.YTickLabel = MLing.ARPA2IPA(m.feats);
                case "artic3"
                    ax.YTickLabel = NP.Artic.GetLabels(m.feats);
                case "strf"
                    ax.YDir = "normal";
                    ax.YTick = [1:20:numel(m.feats) numel(m.feats)];
                    ax.YTickLabel = round(m.feats(ax.YTick)/1e3, 1);
                otherwise
                    ax.YTickLabel = feats;
            end
            MPlot.Axes(ax);
            
            % Plot rasters
            seqTb = tb.seqTb{k};
            if isempty(seqTb)
                continue
            end
            switch sets(i)
                case "phone"
                    ntArgs = MPlot.FindTileInd(nRows, nCols, 1:3, 3);
                    ax = nexttile(ntArgs{:});
                    s2plot = struct('seqTb', seqTb, 'positions', 0, 'nElem', 1, 'level', "phone", 'seed', "NA");
                    LMV.Fig.PlotSeqResp(ax, s2plot, 'UnitIdx', 1);
                case "strf"
                    ntArgs = arrayfun(@(x) MPlot.FindTileInd(nRows, nCols, x, 3), 1:3, 'Uni', false);
                    axArray = cellfun(@(x) nexttile(x{:}), ntArgs);
                    I = seqTb.binIdx(1);
                    LMV.Fig.PlotSpectralResp(M, seqTb, I, fMel(I), 1, axArray);
                case "artic3"
                    featNames = string(seqTb.Properties.VariableNames(end-2:end));
                    for n = 1 : numel(featNames)
                        ntArgs = MPlot.FindTileInd(nRows, nCols, n, 3);
                        ax = nexttile(ntArgs{:});
                        LMV.PE.PlotArticResp(M, seqTb, featNames(n), 'Parent', ax);
                    end
            end
            
            MPlot.Paperize(f, 1, 1);
            return
            exportgraphics(f, figPath);
        end
    end
end
