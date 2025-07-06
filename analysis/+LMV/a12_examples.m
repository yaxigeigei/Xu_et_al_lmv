%% Plot sentence rasters and peri-phone responses for example units

mdlName = LMV.Linker.currentModel;
anaDir = LMV.Data.GetAnalysisDir("linker", mdlName, "examples");
% srcTb = LMV.Data.FindSource([]);

%% Load data

% Linker clusTb with sentence PETHs
load(fullfile(LMV.Data.GetAnalysisDir, "linker", mdlName, "linker_clusTb_peth.mat"), "clusTb");

% Linked positions
load(fullfile(LMV.Data.GetAnalysisDir, "linker", mdlName, "linked_positions.mat"), "posTb");

%% Specify units and sequences

% Get example unit info
types = LMV.Linker.types;
uuu = arrayfun(@(x) LMV.Linker.GetExampleUnitInfo(x), types, 'Uni', false);

% Load seqData based on what the examples need
if ~exist("seqData", "var")
    seqData = {};
end
cidAll = cellfun(@(x) x.clusId, [uuu{:}]);
for i = 1 : numel(cidAll)
    recIdList = string(cellfun(@(x) NP.SE.GetID(x.se), seqData, 'Uni', false));
    recId = string(NP.SE.GetID(cidAll(i)));
    if ~any(recId==recIdList)
        seqData{end+1} = LMV.Linker.LoadSeqData(recId);
    end
end
recIdList = cellfun(@(x) string(NP.SE.GetID(x.se)), seqData);

% Add info
for k = 1 : numel(types)
    uu = uuu{k};
    
    % Find stimId based on seqStr
    stimTextList = struct2array(LMV.Param.stimDict)';
    stimIdList = string(fieldnames(LMV.Param.stimDict));
    for i = 1 : numel(uu)
        isStim = contains(stimTextList, uu{i}.seqStr);
        uu{i}.stimId = stimIdList(isStim);
    end
    
    % Add sequence data
    for i = 1 : numel(uu)
        recId = string(NP.SE.GetID(uu{i}.clusId));
        uu{i}.seq = seqData{recIdList==recId};
    end
    
    % Add sentence responses
    for i = 1 : numel(uu)
        isUnit = clusTb.clusId==uu{i}.clusId;
        s = table2struct(clusTb(isUnit,:));
        uu{i}.peth = s;
    end
    
    % Add linking info
    for i = 1 : numel(uu)
        isUnit = posTb.clusId==uu{i}.clusId;
        uu{i}.posTb = posTb(isUnit,:);
    end
    
    uuu{k} = uu;
end

%% Plot rasters and sequence overlays

for k = 1 : numel(types)
    uu = uuu{k};
    nRows = numel(uu);
    colDist = [4 1];
    
    f = MPlot.Figure(120020+k); clf
    tl = tiledlayout(6, sum(colDist));
    tl.Padding = "compact";
    
    for i = 1 : nRows
        % Raster
        ntArgs = MPlot.FindTileInd(nRows, colDist, i, 1);
        ax = nexttile(ntArgs{:});
        LMV.Fig.LinkerRaster(ax, uu{i}.clusId, 'StimIdList', uu{i}.stimId);
        
        % Overlay
        ntArgs = MPlot.FindTileInd(nRows, colDist, i, 2);
        ax = nexttile(ntArgs{:});
        C = uu{i}.seq.triTb{uu{i}.seed, 1:2};
        uIdx = NP.Unit.ClusId2Ind(uu{i}.clusId, NP.Unit.GetClusTb(uu{i}.seq.se));
        LMV.Fig.SeqRespOverlay(ax, C, 'SeqStr', uu{i}.seqStr, 'UnitIdx', uIdx);
    end
    MPlot.Paperize(f, 1.4, 1.5);
    exportgraphics(f, fullfile(anaDir, sprintf("example_%s_units.png", types(k))));
    print(f, fullfile(anaDir, sprintf("example_%s_units.pdf", types(k))), '-dpdf');
end

%% Plot sentence overlay and linking positions

f = MPlot.Figure(120033); clf

types = LMV.Linker.types;
nRows = 2;
tl = tiledlayout(nRows, numel(types));
tl.Padding = "compact";

for k = 1 : numel(types)
    uu = uuu{k};
    for i = 1 : nRows
        ntArgs = MPlot.FindTileInd(nRows, numel(types), i, k);
        ax = nexttile(ntArgs{:});
        if types(k) == "bridge"
            isSen = any(~isnan(uu{i}.peth.bridgeResp));
        else
            isSen = any(~isnan(uu{i}.peth.mirrorResp));
        end
        LMV.Fig.SentenceRespOverlay(ax, uu{i}.peth, uu{i}.posTb, 'SentenceMask', isSen, "ShowNonLinked", false, "ShowScore", false);
    end
end

MPlot.Paperize(f, 1.8, .6);
exportgraphics(f, fullfile(anaDir, "example_sentence_peth_overlay.png"));
print(f, fullfile(anaDir, "example_sentence_peth_overlay.pdf"), '-dpdf');

return
%% Plot sentence overlay and linking positions

f = MPlot.Figure(120042); clf

types = LMV.Linker.types;
nRows = max(cellfun(@numel, uuu));
tl = tiledlayout(nRows, numel(types));
tl.Padding = "compact";

for k = 1 : numel(types)
    uu = uuu{k};
    for i = 1 : numel(uu)
        ntArgs = MPlot.FindTileInd(nRows, numel(types), i, k);
        ax = nexttile(ntArgs{:});
        if types(k) == "bridge"
            isSen = any(~isnan(uu{i}.peth.bridgeResp));
        else
            isSen = any(~isnan(uu{i}.peth.mirrorResp));
        end
        LMV.Fig.SentenceRespOverlay(ax, uu{i}.peth, uu{i}.posTb, 'SentenceMask', isSen, 'MinScore', 0.01);
        ax.XLabel.String = [];
        ax.XTickLabel = [];
    end
end

% MPlot.Paperize(f, 1.2, 1.6);
% exportgraphics(f, fullfile(anaDir, "example_sentence_peth_overlay_all.png"));
% print(f, fullfile(anaDir, sprintf("example_%s_units.pdf", types(k))), '-dpdf');

%% Plot rasters, sequence overlays, and sentence response overlays

f = MPlot.Figure(120021); clf

types = LMV.Linker.types;
for k = 2 %1 : numel(types)
    uu = uuu{k};
    nRows = numel(uu);
    colDist = [4 1 1.5]*2;
    tl = tiledlayout(nRows, sum(colDist));
    tl.Padding = "compact";
    
    for i = 1 : nRows
        % Raster
        ntArgs = MPlot.FindTileInd(nRows, colDist, i, 1);
        ax = nexttile(ntArgs{:});
        LMV.Fig.LinkerRaster(ax, uu{i}.clusId, 'StimIdList', uu{i}.stimId);
        
        % Overlay
        ntArgs = MPlot.FindTileInd(nRows, colDist, i, 2);
        ax = nexttile(ntArgs{:});
        C = uu{i}.data.triTb{uu{i}.seed, 1:2};
        uIdx = NP.Unit.ClusId2Ind(uu{i}.clusId, NP.Unit.GetClusTb(uu{i}.data.se));
        LMV.Fig.SeqRespOverlay(ax, C, 'SeqStr', uu{i}.seqStr, 'UnitIdx', uIdx);
        
        % Sentence response
        s = uu{i}.senResp;
        if types(k)=="feedback"
            s.lkResp = s.("mirrorResp");
        else
            s.lkResp = s.(types(k)+"Resp");
        end
        ntArgs = MPlot.FindTileInd(nRows, colDist, i, 3);
        ax = nexttile(ntArgs{:});
        LMV.Fig.SentenceRespOverlay(ax, s);
    end
    MPlot.Paperize(f, 1.8, 1.3); return
    % exportgraphics(f, fullfile(anaDir, sprintf("example_%s_units.png", types(k))));
    % print(f, fullfile(anaDir, sprintf("example_%s_units.pdf", types(k))), '-dpdf');
end

%% Plot 14 sentences for selection

f = MPlot.Figure(120001); clf
LMV.Overview.SentencesFromCache(expTb.clusId, 'StimIdList', LMV.Param.stimIdList14);

%% 

clear ss

% Mirror examples
s = struct;
s.clusId = 460100003;
s.data = s4601;
s.seed = "DH";
s.seqStr = ["earth's the matter", "the girl"];
s.stimIdList = ["mbbr0_si2315", "fsjk1_si2285"];
ss(1) = s;

s = struct;
s.clusId = 410100454;
s.data = s4101;
s.seed = "UH";
s.seqStr = ["something pulled my", "something's pulled my", "you took me"];
s.stimIdList = ["mewm0_si1978", "msjs1_si1899"];
ss(2) = s;

% Feedback examples
s = struct;
s.clusId = 440300585;
s.data = s4403;
s.seed = "HH";
s.seqStr = ["he may", "this house"];
s.stimIdList = ["madd0_si1295", "fltm0_si2330"];
ss(3) = s;

% s = struct;
% s.clusId = 440300575;
% s.data = s4403;
% s.seed = "OW";
% s.seqStr = ["was nobody's fault", "i'm going to"];
% s.stimIdList = ["fdxw0_si2141", "fltm0_si2330"];
% ss(4) = s;

s = struct;
s.clusId = 440200479;
s.data = s4402;
s.seed = "OW";
s.seqStr = ["was nobody's fault", "i'm going to"];
s.stimIdList = ["fdxw0_si2141", "fltm0_si2330"];
ss(4) = s;

expTb = struct2table(ss);
disp(expTb)

%% 

for i = 1 : height(expTb)
    C = expTb.data(i).triTb{expTb.seed(i), 1:2};
    uIdx = NP.Unit.ClusId2Ind(expTb.clusId(i), NP.Unit.GetClusTb(expTb.data(i).se));
    
    f = MPlot.Figure(120120+i); clf
    ax = gca;
    LMV.Fig.SeqPethOverlay(ax, C, 'SeqStr', expTb.seqStr{i}, 'UnitIdx', uIdx, 'MaxSpikeRate', 80);
    MPlot.Paperize(f, 'ColumnsHigh', 0.5, 'ColumnsWide', 0.5);
    exportgraphics(f, fullfile(anaDir, sprintf("%i_u%i_seq_overlay.png", i, expTb.clusId(i))));
%     return
end
