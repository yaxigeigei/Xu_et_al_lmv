%% 

anaDir = fullfile(NP.Data.GetAnalysisRoot, 'sent_resp', 'sen4');

recId = "NP41_B1";
recId = "NP44_B3";
% recId = "NP38_B6";

load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', "extracted_resp", recId+"_phase-resp.mat"));
respTb{:,:}(dropoutTb{:,:}) = NaN;

sDec = load(fullfile(anaDir, "computed_mdls", recId+"_clusTb.mat"));
sTest = load(fullfile(anaDir, "computed_tests", recId+"_clusTb.mat"));

%% 

uId = 410100000+[245 224 463 446];

%% 

uId = 440300000+[574 575 585];

%% 

uId = 380600157;

%% 

phaseNames = ["baseline", "atten", "stim", "delay", "init", "prod", "iti"];
% stimIdList = string(stimIdTb.Properties.VariableNames(1:4));
stimIdList = LMV.Param.stimIdList4;
uInd = NP.Unit.ClusId2Ind(uId, clusTb);

f = MPlot.Figure(123); clf
tl = tiledlayout(numel(uInd), 1);
tl.Padding = "compact";

rrr = cell(numel(uInd), numel(phaseNames), numel(stimIdList));
for i = 1 : numel(uInd)
    ax = nexttile;
    
    for k = 1 : numel(phaseNames)
        pn = phaseNames{k};
        phMask = phaseTb.(pn);
        senMask = stimIdTb{:,stimIdList};
        
        rr = cell(size(stimIdList'));
        for n = 1 : numel(stimIdList)
            rr{n} = respTb.(uInd(i))(phMask & senMask(:,n));
        end
        rrr(i,k,:) = rr;
        
        cc = lines;
        for j = 1 : numel(stimIdList)
            y = rr{j};
            dx = repelem(j-1, numel(y))*0.2 + rand(size(y))*0.03;
            x = k-0.3 + dx;
            plot(ax, x, y, 'o', 'Color', cc(j,:));
            hold on
        end
        
        if ~ismember(pn, sTest.clusTb.Properties.VariableNames)
            continue
        end
        
        p = sTest.clusTb.(pn){uInd(i)}.p;
%         if p > 0.05 || isnan(p)
%             continue
%         end
        yText = -7.5;
        text(ax, k, yText, sprintf("p=%.2f", p), 'HorizontalAlignment', 'center');
        
        if isempty(sDec.clusTb.(pn){uInd(i)})
            continue
        end
        a = 1 - sDec.clusTb.(pn){uInd(i)}.loss;
        yText = -15;
        text(ax, k, yText, sprintf("a=%.2f", a), 'HorizontalAlignment', 'center');
    end
    
    ax.XLim = [0.5 numel(phaseNames)+0.5];
    ax.YLim(1) = -20;
    ax.XTick = 1 : numel(phaseNames);
    ax.XTickLabel = phaseNames;
    ax.YLabel.String = "u"+uId(i);
    MPlot.Axes(ax);
end

%% 

rr = squeeze(rrr(1,4,:));
gg = cellfun(@(x,y) repelem(x,numel(y))', num2cell((1:4)'), rr, 'Uni', false);
rr = cat(1, rr{:});
gg = cat(1, gg{:});
p = kruskalwallis(rr, gg, 'off')


%% 

phaseNames = ["atten", "stim", "delay", "init", "prod"];
uInd = NP.Unit.ClusId2Ind(uId, clusTb);

rTb{:,:}(mTb{:,:}) = NaN;

f = MPlot.Figure(456); clf
tl = tiledlayout(numel(uInd), numel(phaseNames));
tl.Padding = "compact";

rrr = cell(numel(uInd), numel(phaseNames));
for i = 1 : numel(uInd)
    for k = 1 : numel(phaseNames)
        pn = phaseNames{k};
        phMask = fTb.(pn);
        stimIdList = LMV.Param.stimIdList4;
        senMask = fTb{:,stimIdList};
        
        rr = cell(size(stimIdList'));
        for n = 1 : numel(stimIdList)
            rr{n} = rTb.(uInd(i))(phMask & senMask(:,n));
        end
        rrr(i,:) = rr(3);
        
        mdl = sDec.clusTb.(pn){uInd(i)};
        if isempty(mdl)
            a = NaN;
        else
            a = 1 - mdl.loss;
        end
        
        s = sTest.clusTb.(pn){uInd(i)};
        if isempty(s)
            p = NaN;
        else
            p = s.KW.p;
        end
        
        ax = nexttile;
        for j = 1 : numel(stimIdList)
            y = rr{j};
            x = repelem(j, numel(y));
            plot(ax, x, y, 'o'); hold on
        end
        ax.XLim = [0.5 4.5];
        ax.YLim(1) = 0;
        ax.XTick = 1 : numel(stimIdList);
        ax.YLabel.String = "u"+uId(i);
        ax.Title.String = sprintf("%s: p=%.2f, a=%.2f", pn, p, a);
        MPlot.Axes(ax);
    end
end

%% 

[p, tbl, stats] = kruskalwallis(cell2mat(rr(2,:)))






%% 

% Load responsiveness
sTTest = load(fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'ttest', "computed_ttest_clusTb.mat"));
sTTest.sigTb = NP.Resp.GetSigTable(sTTest.clusTb);
