%% Compute sparsity score and analyze unit compositions

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'sparsity');
if ~exist(anaDir, 'dir')
    mkdir(anaDir);
end
srcTb = LMV.Data.FindSource('lmv');

%% Load data

% Sentence averaged (M1) data
ceSearch = MBrowse.Dir2Table(fullfile(anaDir, 'ce_m1_sentence-avg', '*_ce_m1.mat'));
ceArray = NP.CE.LoadSession(fullfile(ceSearch.folder, ceSearch.name));

% Load unit responsiveness
rTest = LMV.Resp.LoadPhaseResponseTest();
clusTb = rTest.clusTb;

% Load unit sentence selectivity
sTest = LMV.Resp.LoadSentenceSelectTest([], ["sen4", "sen14"]);

%% Compute phasic scores

cachePath = fullfile(anaDir, "computed_sparsity_clusTb.mat");

if exist(cachePath, 'file')
    % Load computed
    load(cachePath, 'clusTb');
else
    % Construct a center-surround filter
    fs = 1 / ceArray(1).userData.rsOps.rsBinSize;
    
    % Define task phases
    phaseNames = ["atten", "stim", "delay", "init", "prod"];
    
    for i = 1 : numel(ceArray)
        ce = ceArray(i);
        fprintf("\n%s\n", NP.SE.GetID(ce));
        
        % Add task phase events
        phases = struct;
        phases.atten = {'cue1On', 'stimOn'};
        phases.delay = {'stimOff', 'cue3On'};
        phases.init = {'cue3On', 'prodOn'};
        NP.TaskBaseClass.AddEventObjects(ce, phases, 'taskTime');
        
        % Make task phase masks
        ce.Column2Cell("taskTime", phaseNames);
        tt = ce.GetTable("taskTime");
        maskTb = table;
        maskTb.time = ce.GetTable('resp').time;
        for j = 1 : numel(phaseNames)
            pn = phaseNames(j);
            maskTb.(pn) = cellfun(@(e,t) e.MaskTimestamps(t), tt.(pn), maskTb.time, 'Uni', false);
        end
        
        % Find the most repeated sentences
        tv = ce.GetTable('taskValue');
        % [~, I] = sort(tv.numTrial, 'descend');
        % ind = I(1:4);
        ind = tv.numTrial > 2;
        
        % Get responses
        [~, y] = ce.GetArray('resp', ind);
        [~, m] = ce.GetArray(maskTb, ind);
        m = logical(m);
        
        % Compute sparsity scores across all phases
        disp("all");
        [scores, F] = LMV.Sparsity.ComputeScores(y, fs);
        ce.clusTb.F = repmat(F, [height(ce.clusTb) 1]);
        ce.clusTb.sp = scores;
        
        % Compute phase specific sparsity scores
        for j = 1 : numel(phaseNames)
            pn = phaseNames(j);
            disp(pn);
            ce.clusTb.("sp_"+pn) = LMV.Sparsity.ComputeScores(y(m(:,j),:), fs);
        end
    end
    
    % Put sparsity scores in clusTb
    for i = 1 : numel(ceArray)
        ceArray(i).clusTb.mi = [];
    end
    pethTb = cat(1, ceArray.clusTb);
    [~, I] = MMath.SortLike(pethTb.clusId, clusTb.clusId);
    pethTb = pethTb(I,:);
    clusTb = [clusTb pethTb(:,17:end)];
    
    % Save table
    save(cachePath, 'clusTb');
end

return
%% Load data

% Sentence averaged (M1) data
ceSearch = MBrowse.Dir2Table(fullfile(anaDir, 'ce_m1_sentence-avg', '*_ce_m1.mat'));
ceArray = NP.CE.LoadSession(fullfile(ceSearch.folder, ceSearch.name));

% Load unit responsiveness
rTest = LMV.Resp.LoadPhaseResponseTest();
clusTb = rTest.clusTb;

% Load unit sentence selectivity
sTest = LMV.Resp.LoadSentenceSelectTest([], ["sen4", "sen14"]);

% Load TRF models
setName = "combo3pm";
mdlNames = setName + "_" + NP.TRF.targets;
mdlDir = fullfile(NP.Data.GetAnalysisRoot, 'trf', setName, "mdls");
mdlTbs = NP.TRF.LoadModels(mdlNames, mdlDir);
mdlTb = cat(1, mdlTbs{:});
[~, I] = MMath.SortLike(mdlTb.clusId, clusTb.clusId);
mdlTb = mdlTb(I,:);
[trfR, trfP] = NP.TRF.GetModelR(mdlTb{:,mdlNames});
trfR(trfP>=0.05) = NaN;

%% Compute phasic scores

cachePath = fullfile(anaDir, "computed_phasic-score_clusTb.mat");

if exist(cachePath, 'file')
    % Load computed
    load(cachePath, 'clusTb');
else
    % Construct a center-surround filter
    fs = 1 / ceArray(1).userData.rsOps.rsBinSize;
    [~, kerS] = MNeuro.Filter1([], fs, 'gaussian', 0.2);
    [~, kerC] = MNeuro.Filter1([], fs, 'gaussian', 0.05, numel(kerS)/fs);
    kerCS = kerC - kerS;
    kerCS = kerCS / sum(kerCS(kerCS>0));
    
    % Define task phases
    phaseNames = ["atten", "stim", "delay", "init", "prod"];
    
    for i = 1 : numel(ceArray)
        % 
        ce = ceArray(i);
        
        % Add task phase events
        phases = struct;
        phases.atten = {'cue1On', 'stimOn'};
        phases.delay = {'stimOff', 'cue3On'};
        phases.init = {'cue3On', 'prodOn'};
        NP.TaskBaseClass.AddEventObjects(ce, phases, 'taskTime');
        
        % Make task phase masks
        ce.Column2Cell("taskTime", phaseNames);
        tt = ce.GetTable("taskTime");
        maskTb = table;
        maskTb.time = ce.GetTable('resp').time;
        for j = 1 : numel(phaseNames)
            pn = phaseNames(j);
            maskTb.(pn) = cellfun(@(e,t) e.MaskTimestamps(t), tt.(pn), maskTb.time, 'Uni', false);
        end
        
        % Find the most repeated sentences
        tv = ce.GetTable('taskValue');
        % [~, I] = sort(tv.numTrial, 'descend');
        % ind = I(1:4);
        ind = tv.numTrial > 2;
        
        % Get responses
        [~, y] = ce.GetArray('resp', ind);
        [~, m] = ce.GetArray(maskTb, ind);
        m = logical(m);
        
        % Compute AUC ratios
        yC = MNeuro.Filter1(y, fs, 'custom', kerC);
        yCS = MNeuro.Filter1(yC, fs, 'custom', kerCS);
        yCS(yCS<0) = 0;
        
        m = [true(size(m,1),1) m];
        scores = NaN(ce.numResp, numel(phaseNames)+1);
        for j = 1 : numel(phaseNames)+1
            scores(:,j) = sum(yCS(m(:,j),:),1) ./ sum(yC(m(:,j),:),1);
        end
        
        % Save result
        ce.clusTb.phasicScore = scores;
    end
    
    % Put phasic scores in clusTb
    for i = 1 : numel(ceArray)
        ceArray(i).clusTb.mi = [];
    end
    pethTb = cat(1, ceArray.clusTb);
    [~, I] = MMath.SortLike(pethTb.clusId, clusTb.clusId);
    pethTb = pethTb(I,:);
    clusTb = [clusTb pethTb(:,17:end)];
    
    % Put significant TRF r-values in clusTb
    clusTb{:,mdlNames} = trfR;
    
    % Save table
    save(cachePath, 'clusTb');
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

%% Relationship between phasic scores and TRF r-values (interactive)

regions = LMV.Param.regions;
phases = ["stim", "prod", "prod"];
recordings = unique(clusTb.recId);

f = MPlot.Figure(1002); clf
tl = tiledlayout('flow');
tl.Padding = 'compact';
for i = 1 : numel(regions)
    rn = regions(i);
    for j = 1 : numel(phases)
        pn = phases(j);
        mn = mdlNames(j);
        
        ax = nexttile;
        
        isUnit = clusTb.region==rn & sTest.sigTb.(pn)>0;
        if ~any(isUnit)
            continue
        end
        
        x = clusTb.(mn)(isUnit);
        y = clusTb.depth(isUnit)/1000;
        h = plot(x, y, 'w.'); hold on
        h.UserData.forBrush = true;
        ax.UserData.clusTb = clusTb(isUnit,:);
        ax.UserData.taskPhase = pn;
        
        for k = 1 : numel(recordings)
            kn = recordings(k);
            isUnit = clusTb.region==rn & sTest.sigTb.(pn)>0 & clusTb.recId==kn;
            if ~any(isUnit)
                continue
            end
            x = clusTb.(mn)(isUnit);
            y = clusTb.depth(isUnit)/1000;
            h = plot(x, y, '.'); hold on
            h.HandleVisibility = 'off';
        end
        
        ax.XLim = [0 .45];
        ax.YLim = [0 6.5];
        ax.XLabel.String = 'TRF r-value';
        ax.YLabel.String = 'Depth (um)';
        ax.Title.String = sprintf("%s, %s", rn, mn);
        ax.Title.Interpreter = 'none';
        MPlot.Axes(ax);
    end
end
MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', 1.5);
% exportgraphics(f, fullfile(anaDir, "trf-r_vs_depth.png"));

% Initialize brush callback
brushObj = brush(f);
brushObj.ActionPostCallback = @NP.UnitPlot.RasterOnBrush;

return
%% Relationship between phasic scores and modulation index (interactive)

regions = unique(clusTb.region);
phases = ["stim", "prod", "prod"];

f = MPlot.Figure(1001); clf
tl = tiledlayout('flow');
tl.Padding = 'compact';
for i = 1 : numel(regions)
    rn = regions(i);
    for j = 1 : numel(phases)
        pn = phases(j);
        mn = mdlNames(j);
        unitGroups = { ...
            clusTb.region==rn, ... % all units from this region
            clusTb.region==rn & any(clusTb.(pn)<0.05,2), ... % responsive units
            clusTb.region==rn & any(clusTb.(pn)<0.05,2) & ~isnan(clusTb.(mn)), ... % TRF significant units
            };
        cc = [.8 .8 .8; .2 .2 .2; lines];
        
        ax = nexttile;
        for k = 1 : numel(unitGroups)
            isUnit = unitGroups{k};
            x = clusTb.phasicScore(isUnit);
%             y = clusTb.mi4(isUnit);
%             y = clusTb.peakSpkRate(isUnit);
            y = clusTb.depth(isUnit)/1000;
%             y = findgroups(clusTb.recId(isUnit));
            h = plot(x, y, '.', 'Color', cc(k,:)); hold on
            if k == 3
                h.UserData.forBrush = true;
                ax.UserData.clusTb = clusTb(isUnit,:);
                ax.UserData.taskPhase = 'full'; % pn
            else
                h.HitTest = 'off';
            end
        end
        ax.XLim = [0 .8];
        ax.YLim = [0 7.660];
        ax.XLabel.String = 'Phasic score';
        ax.YLabel.String = 'Depth (um)';
        ax.Title.String = sprintf("%s, %s", rn, mn);
        ax.Title.Interpreter = 'none';
        ax.YDir = 'reverse';
        MPlot.Axes(ax);
    end
end
MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', 1.5);
exportgraphics(f, fullfile(anaDir, "phasic-score_vs_depth.png"));

% Initialize brush callback
brushObj = brush(f);
brushObj.ActionPostCallback = @NP.UnitPlot.RasterOnBrush;

