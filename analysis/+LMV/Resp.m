classdef Resp
    
    methods(Static)
        function ss = LoadPhaseResponse(recIds)
            % Load extracted task phase responses from files
            % 
            %   ss = LoadPhaseResponse(recIds)
            % 
            recIds = string(recIds);
            ss = cellfun(@(x) load(fullfile(LMV.Data.GetAnalysisDir, "phase_resp", "extracted_resp", x+"_phase-resp.mat")), recIds);
        end
        
        function ss = LoadPhaseResponseTest(recIds)
            % Load results of task phase responsiveness test from files
            % 
            %   ss = LoadPhaseResponseTest()
            %   ss = LoadPhaseResponseTest(recIds)
            % 
            if nargin < 1
                ss = load(fullfile(LMV.Data.GetAnalysisDir, "phase_resp", "signrank", "computed_test_clusTb.mat"));
            else
                ss = cellfun(@(x) load(fullfile(LMV.Data.GetAnalysisDir, "phase_resp", "signrank", "computed_tests", x+"_clusTb.mat")), string(recIds));
            end
            for i = 1 : numel(ss)
                ss(i).sigTb = LMV.Resp.GetSigTable(ss(i).clusTb);
            end
        end
        
        function ss = LoadSentenceSelectTest(recIds, testNames)
            % Load results of sentence selectivity test from files
            % 
            %   ss = LoadSentenceSelectTest()
            %   ss = LoadSentenceSelectTest(recIds)
            %   ss = LoadSentenceSelectTest(recIds, testNames)
            % 
            % Inputs
            %   recIds          1) One or more recording ID that indicate which recording(s) to load.
            %                   2) If not provided or empty [], the function will load from the concatenated cache.
            %   testNames       One or more test names such as "sen4" and "sen14" indicating which test results to 
            %                   include. If more than one tests are included, elements in the sigTb will take the 
            %                   highest value across all tests, but the clusTb is taken from the first specified 
            %                   test as is. Default is ["sen14", "sen4"].
            % Output
            %   ss              An array of struct that contains the clusTb and sigTb.
            % 
            
            if nargin < 2
                testNames = ["sen14", "sen4"];
            end
            testNames = string(testNames);
            
            if nargin < 1
                recIds = [];
            end
            
            if numel(testNames) > 1
                ssCell = arrayfun(@(x) LMV.Resp.LoadSentenceSelectTest(recIds, x), testNames, 'Uni', false);
                ss = ssCell{1};
                for i = 2 : numel(ssCell)
                    for j = 1 : numel(ss)
                        ss(j).sigTb{:,:} = max(ss(j).sigTb{:,:}, ssCell{i}(j).sigTb{:,:});
                    end
                end
                return
            end
            
            if isempty(recIds)
                ss = load(fullfile(LMV.Data.GetAnalysisDir, "sent_resp", testNames, "computed_test_clusTb.mat"));
            else
                recIds = string(recIds);
                ss = cellfun(@(x) load(fullfile(LMV.Data.GetAnalysisDir, "sent_resp", testNames, "computed_tests", x+"_clusTb.mat")), recIds);
            end
            for i = 1 : numel(ss)
                ss(i).sigTb = LMV.Resp.GetSigTable(ss(i).clusTb);
            end
        end
        
        function s = ExtractPhaseResponses(se)
            % Extract mean spike rate from each task phase, in each trial, and for each unit
            % 
            %   s = ExtractPhaseResponses(se)
            % 
            % Input
            %   se              An MSessionExplorer object.
            % Output
            %   s               A struct with the following tables:
            %   s.respTb        A m-by-n table of mean spike rates. m is the number of periods (# trials x # task phases).
            %                   n is the number of units.
            %   s.dropoutTb     A table matching the elements in respTb, where zero and one indicate the presence and 
            %                   absence (aka dropout), respectively, of a unit during the response period.
            %   s.phaseTb       A m-by-p table of binary task phases indicators. m matches the periods. The p columns are 
            %                   "baseline", "atten", "stim", "delay", "init", "prod", "iti".
            %   s.stimIdTb      A m-by-s table of binary sentence indicators. m matches the periods. The s columns includes 
            %                   all sentences (in stimId) sorted by decreasing number of repeats.
            %   s.stimIdInfoTb  Ouput from NP.TaskBaseClass.SplitBySentence. Rows are sorted by decreasing number of repeats.
            %                   The 'se' column is removed to reduce size.
            % 
            
            % Work on a duplicate
            se = se.Duplicate({'taskTime', 'taskValue', 'spikeTime', 'ni'});
            
            % Add unit spans
            NP.Unit.AddSpikeSpanTable(se);
            
            % Add baseline onset events to se
            se.SetColumn('taskTime', 'baselineOn', se.GetTable('taskTime').cue1On - LMV.Param.respBaselineDur);
            
            % Add task phase events
            phases = struct;
            phases.baseline = {'baselineOn', 'cue1On'};
            phases.atten = {'cue1On', 'stimOn'};
            phases.delay = {'stimOff', 'cue3On'};
            phases.init = {'cue3On', 'prodOn'};
            phases.iti = {'prodOff', 'cue1On'};
            
            seVec = se.Duplicate({'taskTime'});
            seVec.SliceSession(0, 'absolute');
            NP.TaskBaseClass.AddEventObjects(seVec, phases);
            seVec.SliceSession(se.GetReferenceTime, 'absolute');
            se.SetTable('taskTime', seVec.GetTable('taskTime'));
            
            % Compute event spike rates
            phaseNames = ["baseline", "atten", "stim", "delay", "init", "prod", "iti"];
            R = cell(size(phaseNames)); % for spike rates
            M = cell(size(phaseNames)); % for task span
            P = cell(size(phaseNames)); % for task phase
            for p = 1 : numel(phaseNames)
                % Compute event spike rates
                [R{p}, M{p}] = NP.Resp.ComputeEventSpkRate(se, phaseNames(p), LMV.Param.respWinOffset, LMV.Param.respWinMinDur);
                
                % Make task phase mask
                isPhase = false(se.numEpochs, numel(phaseNames));
                isPhase(:,p) = true;
                P{p} = isPhase;
            end
            R = cat(1, R{:});
            M = cat(1, M{:});
            P = cat(1, P{:});
            
            uIds = se.GetTable('spikeTime').Properties.VariableNames;
            s = struct;
            s.respTb = array2table(R, 'VariableNames', uIds);
            s.dropoutTb = array2table(M, 'VariableNames', uIds);
            s.phaseTb = array2table(P, 'VariableNames', phaseNames);
            
            % Make sentence labels
            seTb = NP.TaskBaseClass.SplitBySentence(se);
            seTb = sortrows(seTb, 'numTrial', 'descend');
            seTb.se = []; % make lightweight
            tv = se.GetTable('taskValue');
            isSen = false(se.numEpochs, height(seTb));
            for g = 1 : height(seTb)
                isSen(tv.stimId==seTb.stimId(g), g) = true;
            end
            S = repmat(isSen, [numel(phaseNames) 1]);
            s.stimIdTb = array2table(S, 'VariableNames', seTb.stimId);
            s.stimIdInfoTb = seTb;
        end
        
        function sigTb = GetSigTable(tTb, varargin)
            % A wrapper of NP.Resp.GetSigTable with LMV-specific input arguments
            % 
            %   sigTb = GetSigTable(tTb)
            %   sigTb = GetSigTable(tTb, ...)
            %   sigTb = GetSigTable(tTb, ..., 'VariableNames', ["atten", "stim", "delay", "init", "prod"])
            % 
            % See also NP.Resp.GetSigTable
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('VariableNames', ["atten", "stim", "delay", "init", "prod"], @(x) isstring(x) || iscellstr(x) || ischar(x));
            p.parse(varargin{:});
            varList = string(p.Results.VariableNames);
            
            sigTb = NP.Resp.GetSigTable(tTb, 'VariableNames', varList, p.Unmatched);
        end
        
        function lb = FormatSigString(tb, sentIdx)
            % A helper function that construct a string that labels phase responsiveness
            % 
            %   lb = FormatSigString(tb, sentIdx)
            % 
            
            phaseList = ["atten", "stim", "delay", "init", "prod"];
            phaseList = intersect(phaseList, tb.Properties.VariableNames, 'stable');
            
            alphaList = permute([0.05 0.01 0.001], [1 3 2]);
            
            for n = numel(phaseList) : -1 : 1
                % Get p-values
                pn = phaseList(n);
                P = tb.(pn);
                P = P * size(P,2); % bonferroni correction
                
                % Ignore inhibition for two tailed tests
                rn = pn+"Resp";
                if ismember(rn, tb.Properties.VariableNames)
                    R = tb.(rn);
                    P(R<=0) = NaN;
                end
                
                % Compute significance level
                L = sum(P < alphaList, 3);
                lvl0(:,n) = L(:,1); % all trials
                if sentIdx < size(L,2)
                    lvl(:,n) = L(:,1+sentIdx); % trials of a sentence
                else
                    lvl(:,n) = zeros(size(L(:,1)));
                end
            end
            lvl0 = string(lvl0);
            lvl = string(lvl);
            lvl0 = strrep(lvl0, "0", "\_");
            lvl = strrep(lvl, "0", "\_");
            for n = height(tb) : -1 : 1
                lb(n) = strjoin(lvl0(n,:), '') + "/" + strjoin(lvl(n,:), '');
            end
        end
        
        function PlotViolinSets(rCell, cCell, varargin)
            % 
            % 
            %   PlotViolinSets(rCell, cCell)
            %   PlotViolinSets(rCell, cCell, ..., 'Edges', [])
            %   PlotViolinSets(rCell, cCell, ..., 'GroupNames', [])
            %   PlotViolinSets(rCell, cCell, ..., 'ChildrenNames', [])
            %   PlotViolinSets(rCell, cCell, ..., 'ColorFunc', @lines)
            % 
            
            p = inputParser;
            p.addParameter("Edges", [], @isnumeric);
            p.addParameter("GroupNames", [], @(x) isstring(x) || iscellstr(x));
            p.addParameter("ChildrenNames", [], @(x) isstring(x) || iscellstr(x));
            p.addParameter("Color", [], @isnumeric);
            p.parse(varargin{:});
            chNames = p.Results.ChildrenNames;
            gpNames = p.Results.GroupNames;
            rEdges = p.Results.Edges;
            cc = p.Results.Color;
            
            [nG, nC] = size(rCell);
            xc = meshgrid(1:nG, 1:nC)' + linspace(-.25, .25, nC);
            hh = cell(size(rCell));
            xx = cell(size(rCell));
            rr = cell(size(rCell));
            tt = cell(size(rCell));
            hasUnit = cellfun(@(x) any(~isnan(x)), rCell);
            rMax = max(cat(1, rCell{:}));
            if isempty(cc)
                cc = lines(numel(chNames));
            end
            
            ax = gca;
            for i = 1 : nG
                ys = rMax;
                for j = 1 : nC
                    if ~hasUnit(i,j)
                        continue
                    end
                    
                    % Plot all
                    hh{i,j} = MPlot.ViolinScatter(xc(i,j), rCell{i,j}, rEdges, 'Color', cc(j,:), 'Span', 0.12);
                    hold(ax, 'on');
                    xx{i,j} = hh{i,j}.XData';
                    rr{i,j} = hh{i,j}.YData'; % get the actually plotted
                    % [~, I] = MMath.SortLike(rCell{i,j}, rr{i,j});
                    % tt{i,j} = cCell{i,j}(I,:);
                    
                    for k = j+1 : nC
                        ys = ys + rMax/10;
                        
                        % Test significance
                        if all(isnan(rCell{i,j})) || all(isnan(rCell{i,k}))
                            continue
                        end
                        p = ranksum(rCell{i,j}, rCell{i,k});
                        
                        % Plot marker
                        if p < 0.05
                            c = [0 0 0];
                        else
                            c = [0 0 0]+.7;
                        end
                        plot([xc(i,j) xc(i,k)], [ys ys], 'Color', c);
                        text(xc(i,end), ys, sprintf("%.1e", p), 'FontSize', 8, 'Color', c);
                    end
                end
            end
            
            % Set up brushing callback
            X = cat(1, xx{:});
            R = cat(1, rr{:});
            hCB = plot(X, R, 'w.');
            hCB.UserData.forBrush = true;
            ax.Children = circshift(ax.Children, -1); % move callback plot to the bottom
            ax.UserData.clusTb = cat(1, cCell{:});
            ax.UserData.taskPhase = "full";
            brushObj = brush(gcf);
            brushObj.ActionPostCallback = @NP.UnitPlot.RasterOnBrush;
            
            % Formatting
            ldg = legend([hh{end,:}], chNames);
            ldg.Location = 'eastoutside';
            ldg.EdgeColor = 'none';
            
            ax.XLim = [0.5 nG+.5];
            ax.XTick = 1:nG;
            ax.XTickLabel = gpNames;
            
            MPlot.Axes(ax);
        end
        
        function PlotFractions(ff, varargin)
            % 
            % 
            %   PlotFractions(ff)
            %   PlotFractions(ff, ..., 'GroupNames', [])
            %   PlotFractions(ff, ..., 'ChildrenNames', [])
            %   PlotFractions(ff, ..., 'ColorFunc', @lines)
            % 
            
            p = inputParser;
            p.addParameter("GroupNames", [], @(x) isstring(x) || iscellstr(x));
            p.addParameter("ChildrenNames", [], @(x) isstring(x) || iscellstr(x));
            p.addParameter("Color", [], @isnumeric);
            p.parse(varargin{:});
            chNames = p.Results.ChildrenNames;
            gpNames = p.Results.GroupNames;
            cc = p.Results.Color;
            
            [nG, nC] = size(ff);
            [Fm, Fe] = cellfun(@MMath.MeanStats, ff);
            if isempty(cc)
                cc = lines(numel(chNames));
            end
            
            % Plot mean stats
            ax = gca;
            bb = bar(Fm); hold on
            for i = 1 : numel(bb)
                bb(i).FaceColor = "none";
                bb(i).EdgeColor = cc(i,:);
                errorbar(bb(i).XEndPoints, Fm(:,i), Fe(:,i), 'LineStyle', 'none', 'Marker', 'none', 'Color', cc(i,:));
                plot(bb(i).XEndPoints, cat(2, ff{:,i})', 'x', 'Color', cc(i,:));
            end
            
            % Mark significance
            for i = 1 : nG
                y = 0.4;
                for j = 1 : nC
                    for k = j+1 : nC
                        % Test difference
                        p = ranksum(ff{i,j}, ff{i,k});
                        
                        % Plot indicators
                        y = y + 0.05;
                        xx = [bb(j).XEndPoints(i) bb(k).XEndPoints(i)];
                        if p < 0.05
                            cSig = [0 0 0];
                        else
                            cSig = [0 0 0]+.7;
                        end
                        plot(xx, [y y], 'Color', cSig);
                        text(bb(end).XEndPoints(i), y, sprintf("%.1e", p), 'Color', cSig);
                    end
                end
            end
            
            ax.XTickLabel = gpNames;
            ax.YLabel.String = "Fraction";
            MPlot.Axes(ax);
        end
        
        % Correlation analysis
        function [r, pval] = BootCorr(X, varargin)
            % Compute pairwise correlation coefficients and estimate p-values by bootstrap
            % 
            %   [r, pval] = BootCorr(X)
            %   [r, pval] = BootCorr(X, ..., 'NBoot', 100)
            % 
            % Inputs
            %   X               A m-by-n matrix of task phase responses. m is the number of trials. 
            %                   n is the number of task phases.
            %   'NBoot'         The number of bootstrap iteration to compute the null distributions.
            % Outputs
            %   r               A n-by-n matrix of test correlation coefficients.
            %   pval            A n-by-n matrix of bootstrap p-values for r.
            % 
            
            p = inputParser;
            p.addParameter('NBoot', 100, @isscalar);
            p.parse(varargin{:});
            nBoot = p.Results.NBoot;
            
            r = corr(X, 'Rows', 'pairwise');
            
            rNull = NaN([size(r) nBoot]);
            for i = nBoot : -1 : 1
                Xrand = X;
                for j = 1 : size(X,2)
                    Xrand(:,j) = randsample(X(:,j), size(X,1));
                end
                rNull(:,:,i) = corr(Xrand, 'Rows', 'pairwise');
            end
            
            rFlat = r(:);
            rNullFlat = reshape(rNull, numel(r), nBoot);
            pval = MMath.EstimatePval(rFlat', rNullFlat');
            pval = reshape(pval, size(r));
        end
        
        function PlotPhaseCorrCDF(clusTb, phaseNames)
            % Plot correlation coefficients of responses between pairs of task phases
            % 
            %   PlotPhaseCorrCDF(clusTb, phaseNames)
            % 
            
            % Find regions and task phases
            regions = LMV.Param.regions;
            nRegion = numel(regions);
            nPhase = numel(phaseNames);
            
            % Put coefficients in phase-by-phase-by-unit array
            R = cat(3, clusTb.phaseCorr{:});
            P = cat(3, clusTb.phaseCorrPval{:});
            
            % Plot
            tl = tiledlayout('flow');
            tl.Padding = 'compact';
            
            for i = 1 : nPhase
                for j = i+1 : nPhase
                    ax = nexttile;
                    for k = 1 : nRegion
                        % Select units
                        isReg = clusTb.region==regions(k);
                        r = squeeze(R(i,j,isReg));
                        
                        % Set non-significant coefficients to NaN
                        isSig = squeeze(P(i,j,isReg)) < 0.05;
                        r(~isSig) = NaN;
                        
                        % Compute CDF
                        rEdges = 0:0.005:1;
                        rCenters = MMath.BinEdges2Centers(rEdges);
                        N = histcounts(r, rEdges, 'Normalization', 'cdf');
                        
                        % Plot CDF
                        stairs(rCenters, N, 'Color', NP.Param.GetRegionColors(regions(k)));
                        hold(ax, 'on');
                    end
                    
                    %                 if j == 2
                    %                     ldg = legend(ax, regions);
                    %                     ldg.Location = 'southeast';
                    %                     ldg.EdgeColor = 'none';
                    %                 end
                    
                    ax.XLabel.String = "Pearson's r";
                    ax.YLabel.String = "Frac. of selective";
                    ax.Title.String = strjoin(phaseNames([i j]), '-');
                    MPlot.Axes(ax);
                end
            end
        end
        
        function PlotPhaseCorrScatter(clusTb, phaseNames)
            % Plot correlation coefficients of responses between pairs of task phases
            % 
            %   PlotPhaseCorrScatter(clusTb, phaseNames)
            % 
            
            % Find regions and task phases
            regions = LMV.Param.regions;
            nRegion = numel(regions);
            nPhase = numel(phaseNames);
            
            % Put coefficients in phase-by-phase-by-unit array
            R = cat(3, clusTb.phaseCorr{:});
            P = cat(3, clusTb.phaseCorrPval{:});
            
            % Set up plot
            tl = tiledlayout("flow");
            tl.Padding = 'compact';
            ax = nexttile;
            
            % Preallocate variables to collect the plotted data for labeling and callback
            hh = cell(nPhase, nPhase, nRegion);
            tt = hh;
            xx = hh;
            yy = hh;
            pairLb = hh;
            regLb = hh;
            
            x = 0;
            for k = 1 : nRegion
                rn = regions(k);
                for i = 1 : nPhase
                    for j = i+1 : nPhase
                        % Choose units with significant coefficients and in the region of interest
                        m = squeeze(P(i,j,:)) < 0.05 & clusTb.region==rn;
                        r = squeeze(R(i,j,m));
                        
                        % Plot one violin scatter
                        x = x + 1;
                        hh{i,j,k} = MPlot.ViolinScatter(x, r, 'Span', 0.6, 'Color', NP.Param.GetRegionColors(rn), 'MarkerSize', 8);
                        hold(ax, 'on')
                        
                        % Collect plotting data
                        pairLb{j,i,k} = strjoin(phaseNames([i j]), "-");
                        regLb{i,j,k} = rn;
                        if isempty(hh{i,j,k})
                            continue
                        end
                        xx{i,j,k} = hh{i,j,k}.XData';
                        yy{i,j,k} = r;
                        tt{i,j,k} = clusTb(m,:);
                    end
                end
            end
            
            ax.XLim = [0 x+1];
            ax.YLim = [0 1];
            ax.XTick = 1 : x;
            ax.XTickLabel = [pairLb{:}];
            ax.YLabel.String = "Pearson's r";
            MPlot.Axes(ax);
            
            % Set up brushing callback
            X = cat(1, xx{:});
            Y = cat(1, yy{:});
            hCB = plot(X, Y, 'w.');
            hCB.UserData.forBrush = true;
            ax.Children = circshift(ax.Children, -1); % move callback plot to the bottom
            ax.UserData.clusTb = cat(1, tt{:});
            ax.UserData.taskPhase = "full";
            brushObj = brush(gcf);
            brushObj.ActionPostCallback = @LMV.Overview.SessionOnBrush;
            
            % Formatting
            ldg = legend([hh{1,2,:}], regions);
            ldg.Location = 'eastoutside';
            ldg.EdgeColor = 'none';
        end
        
        % Misc
        function [tW, W] = ComputeModulatedTone(tR, R, freq)
            % 
            
            % Construct sine waves at requested freuencies
            fs = 44100;
            tW = (tR(1) : 1/fs : tR(end))';
            W = arrayfun(@(x) sin(2*pi*x*tW), freq, 'Uni', false);
            W = cat(2, W{:});
            
            % Scale modulation
            % R = max(R-0.1, 0);
            R = R.^2;
            
            % Apply modulation
            R = interp1(tR, R, tW);
            W = W.*R;
            
            % Mix tones
            W = sum(W, 2);
            
            % Add ramps at the beginning and the end
            nSpPad = 0.1*fs;
            W(1:nSpPad) = W(1:nSpPad) .* linspace(0, 1, nSpPad)';
            W(end-nSpPad:end) = W(end-nSpPad:end) .* linspace(1, 0, nSpPad+1)';
        end
        
        % not in use
        function PlotPhaseRespViolin(clusTb, phaseNames)
            % Plot task phase response scores
            % 
            %   PlotPhaseRespViolin(clusTb, phaseNames)
            % 
            
            % Find regions and task phases
            regions = LMV.Param.regions;
            nRegion = numel(regions);
            nPhase = numel(phaseNames);
            
            % Set up plot
            tl = tiledlayout("flow");
            tl.Padding = 'compact';
            ax = nexttile;
            
            % Preallocate variables to collect the plotted data for labeling and callback
            hh = cell(nPhase, nRegion);
            tt = hh;
            xx = hh;
            yy = hh;
            
            x = 0;
            dx = 1;
            for i = 1 : nPhase
                for k = 1 : nRegion
                    rn = regions(k);
                    
                    % Choose units in the region of interest
                    m = clusTb.region==rn & clusTb.tuningValid(:,i);
                    r = clusTb.tuningScore(m,i);
                    
                    % Plot one violin scatter
                    x = x + dx;
                    hh{i,k} = MPlot.ViolinScatter(x, r, 'Span', 0.6, 'Color', NP.Param.GetRegionColors(rn), 'MarkerSize', 8);
                    hold(ax, 'on')
                    
                    % Collect plotting data
                    if isempty(hh{i,k})
                        continue
                    end
                    xx{i,k} = hh{i,k}.XData';
                    yy{i,k} = r;
                    tt{i,k} = clusTb(m,:);
                end
            end
            
            xticks = (1:x)*nRegion - (nRegion-1)/2;
            
            ax.XLim = [0 x+1];
            ax.YLim = [0 1];
            ax.XTick = xticks;
            ax.XTickLabel = phaseNames;
            ax.YLabel.String = "Response score";
            MPlot.Axes(ax);
            
            % Set up brushing callback
            X = cat(1, xx{:});
            Y = cat(1, yy{:});
            hCB = plot(X, Y, 'w.');
            hCB.UserData.forBrush = true;
            ax.Children = circshift(ax.Children, -1); % move callback plot to the bottom
            ax.UserData.clusTb = cat(1, tt{:});
            ax.UserData.taskPhase = "full";
            brushObj = brush(gcf);
            brushObj.ActionPostCallback = @LMV.Overview.SessionOnBrush;
            
            % Formatting
            ldg = legend([hh{1,:}], regions);
            ldg.Location = 'eastoutside';
            ldg.EdgeColor = 'none';
        end
        
        function PlotPairedPhaseResp(clusTb, regions, phasePair)
            % Plot response scores of one phase against another for each unit in a scatter plot
            % 
            %   PlotPairedPhaseResp(clusTb, regions, phasePair)
            % 
            
            % Get phase responses
            [~, I] = MMath.SortLike(clusTb.tuningPhase(1,:), phasePair);
            rr = clusTb.tuningScore(:,I);
            
            % Only exclude units that are not responsive in all phases of the pair
            vv = clusTb.tuningValid(:,I);
            rr(~any(vv,2), :) = NaN;
            
%             % Find units that are responsive to at least one phase
%             isResp = ~all(isnan(rr), 2);
%             rr(isnan(rr) & isResp) = 0; % set non-responsive score to 0
            
            % Plot all points for brushing
            ax = gca;
            range = [-0.1 1];
            plot(range, range, '-', 'Color', [0 0 0]+0.8); hold(ax, 'on');
            h = plot(ax, rr(:,1), rr(:,2), '.', 'Color', 'w');
            h.UserData.forBrush = true;
            
            % Plot points for each region
            hh = cell(size(regions));
            for r = 1 : numel(regions)
                isRegion = clusTb.region==regions(r);
                hh{r} = plot(ax, rr(isRegion,1), rr(isRegion,2), '.', 'Color', NP.Param.GetRegionColors(regions(r)));
            end
            hh = [hh{:}];
            legend(ax, hh, regions, 'Location', 'eastoutside');
            
            axis(ax, 'equal');
            ax.XLim = range;
            ax.YLim = range;
            % ax.XTick = [];
            % ax.YTick = [0.1 1 10];
            ax.XLabel.String = "Response in " + phasePair(1);
            ax.YLabel.String = "Response in " + phasePair(2);
            ax.Title.String = sprintf("%s-%s", phasePair(2), phasePair(1));
            ax.Box = 'off';
            
            % Initialize brush callback
            ax.UserData.clusTb = clusTb;
            brushObj = brush(gcf);
            brushObj.ActionPostCallback = @LMV.Overview.SessionOnBrush;
        end
        
        function PlotDiffPhaseResp(clusTb, regions, phasePair)
            % Plot the distributions of differences between two task phase responses as histograms
            % 
            %   PlotDiffPhaseResp(clusTb, regions, phasePair)
            % 
            
            % Get phase responses
            [~, I] = MMath.SortLike(clusTb.tuningPhase(1,:), phasePair);
            rr = clusTb.tuningScore(:,I);
            
            % Get response mask
            m = false(size(clusTb.tuningScore));
            m = m | clusTb.isResp == 0; % not responsive in any phase of the pair
            m = m | abs(clusTb.tuningMean) < 1; % small effect size
            
            % Set responses to NaN if a unit fails in both phases
            rr(all(m(:,I),2), :) = NaN;
            
            % Get the mask of units responsive to any task phase
            isAnyResp = any(~m, 2);
            
            % Compute difference
            d = rr * [1 -1]';
            
            % Plot all points for brushing
            ax = gca;
            edges = -1:0.1:1;
            for i = 1 : numel(regions)
                m = clusTb.region==regions(i) & isAnyResp; % use # of responsive as total
                hh(i) = histogram(ax, d(m), edges, 'Normalization', 'probability');
                hh(i).FaceColor = NP.Param.GetRegionColors(regions(i));
                hh(i).EdgeColor = 'none';
                hold(ax, 'on');
            end
            
            ax.XLabel.String = "Diff. response";
            ax.XTick = [-1 0 1];
%             ax.XTickLabel = [phasePair(2), "eq.", phasePair(1)];
            ax.YLabel.String = "Frac. units";
            ax.Title.String = sprintf("%s - %s", phasePair(2), phasePair(1));
            ax.Box = 'off';
        end
        
        function PlotSelectivity(tTb, regions, phasePair)
            % Plot spike rates of one phase against another for each unit as a scatter plot
            % 
            %   PlotSelectivity(tTb, regions, phasePair)
            % 
            
            % Find units responsive to at least one of the phase pair
            tId = tTb{:, "tId"+(1:5)};
            tR = tTb{:, "tR"+(1:5)};
            nUnits = height(tTb);
            rr = NaN(nUnits, 2);
            for i = 1 : nUnits
                for p = 1 : 2
                    j = find(tId(i,:)==phasePair(p), 1);
                    if ~isempty(j)
                        rr(i,p) = tR(i,j);
                    end
                end
            end
            isResp = any(~isnan(rr), 2);
            isJitter = isnan(rr);
            nJitter = sum(isJitter, 'all');
            rr(isnan(rr)) = 10.^(rand(nJitter,1)) / 100;
            
            % Plot all points for brushing
            ax = nexttile;
            range = [0.05 50];
            plot(range, range, '-', 'Color', [0 0 0]+0.8); hold(ax, 'on');
            plot([0.1 0.1], range, '-', 'Color', [0 0 0]+0.8);
            plot(range, [0.1 0.1], '-', 'Color', [0 0 0]+0.8);
            h = plot(ax, rr(isResp,1), rr(isResp,2), '.', 'Color', 'w');
            h.UserData.forBrush = true;
            
            for r = 1 : numel(regions)
                isRegion = isResp & tTb.region==regions(r);
                hh(r) = plot(ax, rr(isRegion,1), rr(isRegion,2), '.', 'Color', NP.Param.GetRegionColors(regions(r)));
            end
            legend(ax, hh, regions, 'Location', 'eastoutside');
            
            axis(ax, 'equal');
            ax.XScale = 'log';
            ax.YScale = 'log';
            ax.XLim = range;
            ax.YLim = range;
            ax.XTick = [];
            ax.YTick = [0.1 1 10];
            ax.XLabel.String = "\Delta spike rate by " + phasePair(1);
            ax.YLabel.String = "\Delta spike rate by " + phasePair(2);
            ax.Title.String = sprintf("%s-%s", phasePair(2), phasePair(1));
            ax.Box = 'off';
            
            % Initialize brush callback
            ax.UserData.clusTb = tTb(isResp,:);
            ax.UserData.taskPhase = "full";
            brushObj = brush(gcf);
            brushObj.ActionPostCallback = @NP.UnitPlot.RasterOnBrush;
        end
        
        function PlotDiffSelectivity(tTb, regions, phasePair)
            % Plot the distributions of differences between two task phase responses as histograms
            % 
            %   PlotDiffSelectivity(tTb, regions, phasePair)
            % 
            
            % Find units responsive to at least one of the phase pair
            tId = tTb{:, "tId"+(1:5)};
            tR = tTb{:, "tR"+(1:5)};
            nUnits = height(tTb);
            rr = NaN(nUnits, 2);
            for i = 1 : nUnits
                for p = 1 : 2
                    j = find(tId(i,:)==phasePair(p), 1);
                    if ~isempty(j)
                        rr(i,p) = tR(i,j);
                    end
                end
            end
            isResp = any(~isnan(rr), 2);
            rr(isnan(rr)) = 0;
            d = rr * [1 -1]';%./norm([1 -1]);
            
            % Plot all points for brushing
            ax = nexttile;
            
            edges = -15:1:15;
            for r = 1 : numel(regions)
                isRegion = isResp & tTb.region==regions(r);
                hh(r) = histogram(ax, d(isRegion), edges, 'Normalization', 'probability');
                hh(r).FaceColor = NP.Param.GetRegionColors(regions(r));
                hh(r).EdgeColor = 'none';
                hold(ax, 'on');
            end
            legend(ax, hh, regions, 'Location', 'northeast');
            
            ax.XLabel.String = "Difference in spike rates";
            ax.YLabel.String = "# of units";
            ax.Title.String = sprintf("%s-%s", phasePair(2), phasePair(1));
            ax.Box = 'off';
        end
        
    end
    
end
