%% Anatomical and electrophysiological correlations with firing sparsity

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'sparsity');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end
srcTb = LMV.Data.FindSource('lmv');

%% Load data

% Sparsity scores
load(fullfile(anaDir, "computed_sparsity_clusTb.mat"), 'clusTb');

% Tests
% resp = LMV.Resp.LoadPhaseResponse(srcTb.recId);
rTest = LMV.Resp.LoadPhaseResponseTest();
sTest = LMV.Resp.LoadSentenceSelectTest();

%% Sparsity profile acoross 

rTest.sigTb = LMV.Resp.GetSigTable(rTest.clusTb, 'MinRespMag', 1);

regions = LMV.Param.regions;
% regions = "mPrCG";
phases = ["stim", "delay", "init", "prod"];

for i = 1 %: numel(regions)
    rn = regions(i);
    
    f = MPlot.Figure(2001); clf
    tl = tiledlayout(1, numel(phases)+1);
    tl.Padding = 'compact';
    
    for j = 1 : numel(phases)+1
        if j <= numel(phases)
            pn = phases(j);
            isResp = clusTb.region==rn & rTest.sigTb.(pn)>0; % responsive units
            isSelect = clusTb.region==rn & sTest.sigTb.(pn)>0; % sentence selective units
            vn = "sp_"+pn;
        else
            rn = regions(i);
            pn = "all";
            isResp = clusTb.region==rn & any(rTest.sigTb{:,:}, 2); % responsive units
            isSelect = clusTb.region==rn & any(sTest.sigTb{:,:}, 2); % sentence selective units
            vn = "sp";
        end
        isUnit = isSelect; % isResp & 
        
        subTb = clusTb(isUnit,:);
        subTb = sortrows(subTb, 'depth');
        
        F = flip(subTb.F(1,:));
        S = flip(subTb.(vn), 2);
        d = subTb.depth;
        
%         S = MMath.Normalize(S', 'zscore')';
        
        ax = nexttile;
%         pcolor(F, d/1000, S);
%         shading flat
%         ax.Title.String = sprintf("%s %s", rn, pn);
%         ax.XLabel.String = "Frequency (Hz)";
%         ax.XScale = 'log';
%         ax.XTick = [0.5 1 2 4 8 16];
%         ax.YLabel.String = "Depth (mm)";
%         ax.YDir = 'reverse';
        
        h = imagesc(ax, S);
        ax.Title.String = sprintf("%s %s", rn, pn);
        ax.XLabel.String = "Frequency (Hz)";
        ax.XTick = round(linspace(1, numel(F), 5));
        ax.XTickLabel = round(F(ax.XTick), 1);
        ax.YLabel.String = "Units";
        ax.YDir = 'reverse';
        MPlot.Axes(ax);
        
        ax.UserData.clusTb = subTb;
        h.ButtonDownFcn = @LMV.Sparsity.ClickHeatmap;
        h.BusyAction = 'cancel';
        
%         MPlot.PlotTraceLadder(F, -S'*20, d/1000);
%         MPlot.Axes(ax);
%         ax.Title.String = sprintf("%s %s", rn, pn);
%         ax.XLabel.String = "Frequency (Hz)";
%         ax.XScale = 'log';
%         ax.XTick = [0.5 1 2 4 8 16];
%         ax.YLabel.String = "Depth (mm)";
%         ax.YDir = 'reverse';
    end

    MPlot.Paperize(f, 'ColumnsWide', 1.8, 'ColumnsHigh', 1);
%     exportgraphics(f, fullfile(anaDir, sprintf("sparsity_by_depth_%s.png", rn)));
end


%% Relationship between phasic score and cortical depth (interactive)

regions = LMV.Param.regions;
phases = ["stim", "delay", "init", "prod"];

f = MPlot.Figure(1001); clf
tl = tiledlayout(numel(regions), numel(phases));
tl.Padding = 'compact';
for i = 1 : numel(regions)
    for j = 1 : numel(phases)
        rn = regions(i);
        pn = phases(j);
        
        unitGroups = { ...
            clusTb.region==rn & rTest.sigTb.(pn)>0, ... % responsive units
            clusTb.region==rn & sTest.sigTb.(pn)>0, ... % sentence selective units
            };
        
        recordings = unique(clusTb.recId(clusTb.region==rn));
        for k = 1 : numel(recordings)
            unitGroups{end+1} = unitGroups{2} & clusTb.recId==recordings(k);
        end
        
        cc = [.9 .9 .9; 1 1 1; lines];
        
        ax = nexttile;
        for k = 1 : numel(unitGroups)
            isUnit = unitGroups{k};
%             if ~any(isUnit)
%                 continue
%             end
            x = clusTb.phasicScore(isUnit,1+j);
            y = clusTb.depth(isUnit)/1000;
            h = plot(x, y, '.', 'Color', cc(k,:)); hold on
            
            if k == 2
                h.UserData.forBrush = true;
                ax.UserData.clusTb = clusTb(isUnit,:);
                ax.UserData.taskPhase = 'full'; % pn
            else
%                 h.HandleVisibility = 'off';
            end
        end
        
        if j == numel(phases)
            legend(flip(ax.Children(1:end-2)), recordings, 'Location', 'eastoutside', 'Interpreter', 'none', 'Box', 'off');
        end
        
        ax.XLim = [0 .8];
        ax.YLim = [0 6.5];
        ax.XLabel.String = 'Phasic score';
        ax.YLabel.String = 'Depth (um)';
        ax.Title.String = sprintf("%s, %s", rn, pn);
        ax.Title.Interpreter = 'none';
        ax.YDir = 'reverse';
        MPlot.Axes(ax);
    end
end
MPlot.Paperize(f, 'ColumnsWide', 1.2, 'ColumnsHigh', 1.5);
exportgraphics(f, fullfile(anaDir, "phasic-score_vs_depth.png"));

% Initialize brush callback
brushObj = brush(f);
brushObj.ActionPostCallback = @LMV.Overview.SessionOnBrush;
