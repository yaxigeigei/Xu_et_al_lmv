%% Generate a library of sequence plots

cacheDir = LMV.Data.GetAnalysisDir("peri_event", "extracted_seq");
srcTb = LMV.Data.FindSource([]);

%% Plot sequence tiers as legend

% Load cached sequence data
recId = "NP44_B2";
s = load(fullfile(cacheDir, recId+"_seqData.mat"), "se", "feats", "seqData");
seeds = s.feats;

% Loop through targets
targets = LMV.TRF.targets;
for tIdx = 1 : numel(targets)
    % Get word sequences
    target = targets(tIdx);
    Cw = s.seqData.(target){"word"}{:,:};
    
    % Create folder
    figDir = fullfile(LMV.Data.GetAnalysisDir, "peri_event", "lib_"+target, recId);
    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end
    
    for i = 1 : numel(seeds)
        % 
        figName = sprintf("%s_%i_%s.png", recId, i, seeds(i));
        figPath = fullfile(figDir, figName);
        if exist(figPath, 'file')
            continue
        end
        
        f = MPlot.Figure(234); clf
        LMV.PE.PlotPeriEventSeqArray(Cw(i,:), 'FeatureName', "tgtiers");
        exportgraphics(f, figPath);
%         return
    end
end

%% Plot phoneme tuning, sequence response heatmap, PETHs, feature timeseries for each content unit

for recIdx = 1 : height(srcTb)
    % Load cached sequence data
    recId = NP.SE.GetID(srcTb.name{recIdx});
    s = load(fullfile(cacheDir, recId+"_seqData.mat"), "se", "feats", "seqData");
    seeds = s.feats;
    se = s.se;
    clusTb = NP.Unit.GetClusTb(se);
    
    % Load phone ti-RFs
    targets = LMV.TRF.targets;
    trfTb = LMV.TRF.LoadModels("phone_"+targets, recId);
    [~, trfInd] = MMath.SortLike(trfTb.clusId, clusTb.clusId);
    
    % Loop through targets
    for tIdx = 1 : numel(targets)
        % Extract responses and features from se
        target = targets(tIdx);
        switch target
            case 'stim'
                tWin = [-0.1 0.4];
                phaseName = 'stim';
            case 'prod'
                tWin = [-0.4 0.1];
                phaseName = 'prod';
            case 'feedback'
                tWin = [-0.1 0.4];
                phaseName = 'prod';
            otherwise
                error("'%s' is not a supported target name.", target);
        end
        
        Cp = s.seqData.(target){"phone"}{:,:};
        
        for k = 1 : numel(Cp)
            seqTb = Cp{k}.seqTb;
            if isempty(seqTb)
                continue
            end
            seqTb(seqTb.nSample<3,:) = [];
            seqTb = LMV.PE.AddResp2SeqTb(seqTb, se, tWin);
            seqTb = LMV.PE.AddFeat2SeqTb(seqTb, se, flip(-tWin), "mel", "mic");
            Cp{k}.seqTb = seqTb;
        end
        
        % Mask out sequences with no further branching
        nSeq = cellfun(@(x) height(x.seqTb), Cp);
        for k = size(nSeq,2) : -1 : 2
            m = nSeq(:,k) <= nSeq(:,k-1);
            nSeq(m,k) = NaN;
        end
        
        
        % Loop through units
        for uIdx = 1 : height(clusTb)
            % Get the optimal encoding time
            mdl = trfTb.("phone_"+target){trfInd(uIdx)};
            if isempty(mdl)
                continue
            end
            [mdl.r2, mdl.r2t, mdl.r2Idx] = LMV.RF.FindPeakR2(mdl.r2each, mdl.dt);
            tResp = Cp{1}.seqTb.tResp{1};
            [~, iOpt] = min(abs(tResp - mdl.r2t));
            
            % Use a default +-0.15s if too little variance explained or optimal model is not significant
            if mdl.r2 < 0.01 || mdl.null.r2Pval > 0.05
                iOpt = round(numel(tResp)/2);
            end
            
            % Find the spike rate at optimal encoding time
            R = NaN(size(Cp));
            for i = 1 : numel(Cp)
                seqTb = Cp{i}.seqTb;
                if isempty(seqTb)
                    continue
                end
                R(i) = max(cellfun(@(y) y(iOpt,uIdx), seqTb.hh));
            end
            R(isnan(nSeq)) = NaN;
            
            % Find the max response for each seed and sort seeds in descending order
            maxR = max(R, [] ,2);
            [maxR, iR] = sort(maxR, 'descend');
            
            
            % Create folder
            uidStr = "u" + clusTb.clusId(uIdx);
            figDir = fullfile(LMV.Data.GetAnalysisDir, "peri_event", "lib_"+target, uidStr);
            if exist(figDir, 'dir')
%                 fprintf("\nNot generating figures in existing folder:\n%s\n", figDir);
%                 continue
            else
                mkdir(figDir);
            end
            
            % Plot TRF weights as a heatmap
            f = MPlot.Figure(985); clf
            ax = nexttile;
            LMV.TRF.PlotWeights2(mdl, 'ClusInfo', clusTb(uIdx,:), 'Parent', ax);
            MPlot.Axes(ax);
            MPlot.Paperize(f, 'ColumnsWide', 0.4, 'ColumnsHigh', 1.4);
            exportgraphics(f, fullfile(figDir, uidStr+"_phone_trf.png"));
            
            % Plot RF weights at optimal time
            f = MPlot.Figure(986); clf
            ax = nexttile;
            rfMdl = mdl.mdls{mdl.r2Idx};
            rfMdl.feats = mdl.feats;
            LMV.RF.PlotWeights(rfMdl, 'ClusInfo', clusTb(uIdx,:), 'Parent', ax);
            MPlot.Axes(ax);
            MPlot.Paperize(f, 'ColumnsWide', 1.5, 'ColumnsHigh', 0.4);
            exportgraphics(f, fullfile(figDir, uidStr+"_phone_rf.png"));
            
            % Plot tuning
            f = MPlot.Figure(987); clf
            h = heatmap(1:3, MLing.ARPA2IPA(seeds(iR)), round(R(iR,:)));
            h.Title = uidStr;
            h.Colormap = colormap('copper');
            MPlot.Paperize(f, 'ColumnsWide', 0.4, 'ColumnsHigh', 1.4);
            exportgraphics(f, fullfile(figDir, uidStr+"_seq_tuning.png"));
            
            % Plot sequence responses
            for i = 1 : numel(iR)
                k = iR(i);
                s2plot = Cp(k,:);
                
                figName = sprintf("%s_%i_%s.png", uidStr, i, seeds(k));
                figPath = fullfile(figDir, figName);
                if exist(figPath, 'file')
                    continue
                end
                
                f = MPlot.Figure(234); clf
                LMV.PE.PlotPeriEventSeqArray(s2plot, 'UnitIdx', uIdx);
                exportgraphics(f, figPath);
            end
            
            % Plot sequence feature
            featName = "mel";
            for i = 1 : numel(iR)
                k = iR(i);
                s2plot = Cp(k,:);
                
                figName = sprintf("%s_%i_%s.png", featName, i, seeds(iR(i)));
                figPath = fullfile(figDir, figName);
                if exist(figPath, 'file')
                    continue
                end
                
                f = MPlot.Figure(345); clf
                LMV.PE.PlotPeriEventSeqArray(s2plot, "Feature", featName);
                exportgraphics(f, fullfile(figDir, figName));
            end
        end
    end
end

