%% 

figDir = "C:\chang_lab\project_np\misc\20240425_simon_keynote";

%% 

cid = [460100003, 410100454, 440300585];
sid = ["mbbr0_si2315", "mewm0_si1978", "fltm0_si2330"];

%% Load sequence responses

s4101 = LMV.Linker.LoadSeqData('NP41_B1');
s4402 = LMV.Linker.LoadSeqData('NP44_B2');
s4403 = LMV.Linker.LoadSeqData('NP44_B3');
s4601 = LMV.Linker.LoadSeqData('NP46_B1');

%% 

clear ss

% Mirror
s = struct;
s.clusId = 460100003;
s.data = s4601;
s.seed = "DH";
s.seqStr = ["earth's the matter", "the girl"];
s.stimIdList = ["mbbr0_si2315", "fsjk1_si2285"];
ss(1) = s;

% Bridge
s = struct;
s.clusId = 410100454;
s.data = s4101;
s.seed = "UH";
s.seqStr = ["something pulled my", "something's pulled my", "you took me"];
s.stimIdList = ["mewm0_si1978", "msjs1_si1899"];
ss(2) = s;

% Feedback
s = struct;
s.clusId = 440300585;
s.data = s4403;
s.seed = "HH";
s.seqStr = ["this house", "he may"];
s.stimIdList = ["fltm0_si2330", "madd0_si1295"];
ss(3) = s;

expTb = struct2table(ss);
disp(expTb)

%% Plot 14 sentences for selection

f = MPlot.Figure(120001); clf
LMV.Overview.SentencesFromCache(expTb.clusId, 'StimIdList', LMV.Param.stimIdList14);

%% 

for i = 1 : height(expTb)
    f = MPlot.Figure(120020+i); clf
    LMV.Linker.FigTrialRasters(expTb.clusId(i), 'StimIdList', expTb.stimIdList(i,1));
    MPlot.Paperize(f, 'ColumnsHigh', 0.4, 'ColumnsWide', 2);
    % exportgraphics(f, fullfile(figDir, sprintf("%i_u%i_trial_rasters.png", i, expTb.clusId(i))));
    % return
end

%% 

for i = 1 : height(expTb)
    C = expTb.data(i).triTb{expTb.seed(i), 1:2};
    uIdx = NP.Unit.ClusId2Ind(expTb.clusId(i), NP.Unit.GetClusTb(expTb.data(i).se));
    
    f = MPlot.Figure(120120+i); clf
    ax = gca;
    LMV.Linker.FigSeqPethOverlay(ax, C, 'SeqStr', expTb.seqStr{i}(1), 'UnitIdx', uIdx, 'MaxSpikeRate', 80);
    MPlot.Paperize(f, 'ColumnsHigh', 0.4, 'ColumnsWide', 0.5);
    exportgraphics(f, fullfile(figDir, sprintf("%i_u%i_seq_overlay.png", i, expTb.clusId(i))));
    % return
end

%% 










