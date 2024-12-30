classdef Resp
    
    methods(Static)
        % Task phase responsiveness
        function [r, isOut] = ComputeEventSpkRate(se, eventName, tOffset, minDur)
            % Compute mean spike rates during events
            % 
            %   r = ComputeEventSpkRate(se, eventName, tOffset, minDur)
            % 
            
            if nargin < 4
                minDur = 0;
            end
            if nargin < 3
                tOffset = 0;
            end
            
            % Extract time windows
            evts = se.GetColumn('taskTime', eventName);
            tWin = NaN(se.numEpochs, 2);
            for n = 1 : se.numEpochs
                if iscell(evts)
                    evt = evts{n};
                else
                    evt = evts(n);
                end
                if isnan(evt)
                    continue
                end
                if class(evt) == "MSessionExplorer.Event"
                    tWin(n,:) = [evt(1).GetTfield('tOn') evt(end).GetTfield('tOff')];
                else
                    tWin(n,:) = [evt(1).GetTfield('tmin') evt(end).GetTfield('tmax')];
                end
            end
            
            % Refine time windows
            switch eventName
                case 'delay'
                    tWin = tWin + [tOffset, 0];
                case 'init'
                    tWin = tWin + [0, -tOffset];
                case 'retrieval'
                    tWin = tWin + [tOffset, -tOffset];
            end
            
            % Report status
            isValid = diff(tWin,1,2) > minDur; % this also excludes NaN time windows
            tWin(~isValid,:) = NaN;
            fprintf("%s: %i/%i valid time windows.\n", eventName, sum(isValid), numel(isValid));
            
            % Compute spike rates over windows
            st = se.SliceEventTimes('spikeTime', tWin, isValid);
            dur = diff(tWin(isValid,:), 1, 2);
            dur = repmat(dur, [1 width(st)]);
            r = NaN(se.numEpochs, width(st));
            r(isValid,:) = cellfun(@(k,d) sum(~isnan(k))/d, st{:,:}, num2cell(dur));
            
            % Set responses when units are dropped out to NaN
            if ~ismember('spikeSpan', se.tableNames)
                return
            end
            ss = se.SliceTimeSeries('spikeSpan', tWin, isValid);
            isOut = false(size(r));
            isOut(isValid,:) = cellfun(@(x) ~any(x), ss{:,2:end});
        end
        
        function [dataTb, I] = Circshift(dataTb, groupTb)
            % Randomly shift the spike rates to maintain the temporal dynamics but breaking relationship to trials
            % 
            %   dataTb = Circshift(dataTb)
            %   dataTb = Circshift(dataTb, groupTb)
            % 
            
            if nargin < 2
                % Treat all rows as one group
                groupTb = table(true(height(dataTb), 1));
            end
            
            % Circshift row indices by group
            I = zeros(height(dataTb), 1);
            for p = 1 : width(groupTb)
                ind0 = find(groupTb.(p));
                I(ind0) = circshift(ind0, randi(numel(ind0)-1));
            end
            
            % Keep the order of out-of-group rows
            I(I==0) = find(I==0);
            
            % Apply shifts
            dataTb = dataTb(I,:);
        end
        
        function tTb = PhaseTTest(rTb, fTb, phaseNames, stimIdList)
            % Run t-test for each unit, each task phase, across all trials, and across trials of each stentence.
            % 
            %   tTb = PhaseTTest(se, phaseNames)
            % 
            % Inputs
            %   rTb         A m-by-n table of mean spike rates. m is the number of periods (# trials x # task phases).
            %               n is the number of units.
            %   pTb         A m-by-p table of binary task phases indicators. m matches the periods of rTb. p is 
            %               the task phases.
            %   sTb         A m-by-s table of binary sentence indicators. m matches the periods of rTb. s is the 
            %               stim or, more generally, groups.
            % Output
            %   tTb         A table of test reuslts. Rows are units. Columns are task phases. Five p-values are 
            %               computed for each unit and each task phase. The first is across all trials in se. The 
            %               rest are for each of the sentence in LMV.Param.stimIdList4.
            %               Positive p-values indicate activation, negative ones indicate suppression. 
            %               Absence or insufficient trials for a given sentence will leave the correspoding columns 
            %               of p-values to be NaN. If there is no trial where test spike rate differ from baseline, 
            %               ttest will return NaN.
            % 
            
            if ~exist('stimIdList', 'var')
                stimIdList = LMV.Param.stimIdList4;
            end
            
            % Get baseline responses
            m = logical(fTb.("baseline"));
            rBase = rTb{m,:};
            
            % Compute mean firing rate for baseline and each task phase
            phaseNames = setdiff(phaseNames, "baseline");
            nPhase = numel(phaseNames);
            for p = nPhase : -1 : 1
                m = logical(fTb.(phaseNames{p}));
                rPhase(:,:,p) = rTb{m,:};
            end
            rDiff = rPhase - rBase;
            
            % Test with all trials
            nUnit = width(rTb);
            nSent = numel(stimIdList);
            tTb = table;
            for p = 1 : nPhase
                pn = phaseNames{p};
                tTb.(pn) = NaN(nUnit, 1+nSent);
                [~, tTb.(pn)(:,1)] = ttest(rDiff(:,:,p));
                mu = mean(rDiff(:,:,p), 1, 'omitnan');
                tTb.(pn+"Resp")(:,1) = mu;
            end
            
            % Test by sentences
            for s = 1 : nSent
                sid = stimIdList(s);
                if ~ismember(sid, fTb.Properties.VariableNames)
                    fprintf("Stim #%i %s is not found in this recording\n", s, sid);
                    continue
                end
                
                % Test with sentence repeats
                for p = 1 : nPhase
                    pn = phaseNames{p};
                    
                    m = fTb.(sid) & fTb.("baseline");
                    rBase = rTb{m,:};
                    
                    m = fTb.(sid) & fTb.(pn);
                    rPhase = rTb{m,:};
                    
                    rDiff = rPhase - rBase;
                    
                    for u = 1 : size(rDiff,2)
                        val = rDiff(:,u);
                        if sum(~isnan(val)) < 3
                            continue
                        end
                        [~, pval] = ttest(val);
                        tTb.(pn)(u,1+s) = pval;
                    end
                    mu = mean(rDiff, 1, 'omitnan');
                    tTb.(pn+"Resp")(:,1+s) = mu;
                end
            end
        end
        
        function tTb = PhaseSignrank(respTb, phaseTb, groupTb, groupNames)
            % Run Wilcoxon signed rank test on each unit and each task phase, across 1) all trials, and 2) repeats of 
            % each stim or group.
            % 
            %   tTb = PhaseSignrank(respTb, phaseTb, groupTb)
            %   tTb = PhaseSignrank(respTb, phaseTb, groupTb, groupVars)
            % 
            % Inputs
            %   respTb      A m-by-n table of mean spike rates. m is the number of periods (# trials x # task phases).
            %               n is the number of units.
            %   phaseTb     A m-by-p table of binary task phases indicators. m matches the periods of rTb. p is 
            %               the number of task phases.
            %   groupTb     A m-by-s table of binary sentence indicators. m matches the periods of rTb. s is the 
            %               number of stim or, more generally, groups.
            %   groupVars   
            % Output
            %   tTb         A table of test reuslts. Rows are units. Columns are task phases. Five p-values are 
            %               computed for each unit and each task phase. The first is across all trials in se. The 
            %               rest are for each of the sentence in LMV.Param.stimIdList4. 
            %               Positive p-values indicate activation, negative ones indicate suppression. 
            %               Absence or insufficient trials for a given sentence will leave the correspoding columns 
            %               of p-values to be NaN. If there is no trial where test spike rate differ from baseline, 
            %               ttest will return NaN.
            % 
            
            % Get baseline responses
            m = logical(phaseTb.baseline);
            rBase = respTb{m,:};
            
            % Compute mean firing rate for baseline and each task phase
            phaseNames = phaseTb.Properties.VariableNames;
            phaseNames = setdiff(phaseNames, "baseline");
            nPhase = numel(phaseNames);
            for p = nPhase : -1 : 1
                m = logical(phaseTb.(phaseNames{p}));
                rPhase(:,:,p) = respTb{m,:};
            end
            rDiff = rPhase - rBase;
            
            % Test with all trials
            if ~exist('groupNames', 'var') || isempty(groupNames)
                groupNames = string(groupTb.Properties.VariableNames);
            end
            nGroups = numel(groupNames);
            nUnit = width(respTb);
            tTb = table;
            for p = 1 : nPhase
                pn = phaseNames{p};
                tTb.(pn) = NaN(nUnit, 1+nGroups);
                for u = 1 : nUnit
                    dr = rDiff(:,u,p);
                    if all(isnan(dr))
                        continue
                    end
                    tTb.(pn)(u,1) = signrank(dr);
                end
                mu = mean(rDiff(:,:,p), 1, 'omitnan');
                tTb.(pn+"Resp")(:,1) = mu;
            end
            
            % Test by sentences
            for g = 1 : nGroups
                gn = groupNames(g);
                if ~ismember(gn, groupTb.Properties.VariableNames)
                    fprintf("Group #%i %s is not found in this recording\n", g, gn);
                    continue
                end
                
                % Test with sentence repeats
                for p = 1 : nPhase
                    pn = phaseNames{p};
                    
                    m = groupTb.(gn) & phaseTb.baseline;
                    rBase = respTb{m,:};
                    
                    m = groupTb.(gn) & phaseTb.(pn);
                    rPhase = respTb{m,:};
                    
                    rDiff = rPhase - rBase;
                    
                    for u = 1 : size(rDiff,2)
                        dr = rDiff(:,u);
                        if all(isnan(dr))
                            continue
                        end
                        tTb.(pn)(u,1+g) = signrank(dr);
                    end
                    mu = mean(rDiff, 1, 'omitnan');
                    tTb.(pn+"Resp")(:,1+g) = mu;
                end
            end
        end
        
        function zTb = PhaseZeta(se, phaseNames)
            % Run ZETA test for each unit, on each task phase, across all trials and across trials of each stentence
            % 
            %   zTb = PhaseZeta(se, phaseNames)
            % 
            
            warning('off', 'zetatest:InsufficientSamples');

            % Remove trials with bad performance
            se = se.Duplicate({'taskTime', 'taskValue', 'spikeTime'});
            isBad = LMV.SE.IsBadTrials(se, -Inf, 0.5);
            se.RemoveEpochs(isBad);
            
            % Vectorize recording
            se0 = se.Duplicate;
            se0.SliceSession(0, 'absolute');
            
            % Add task phase events to vectorized se
            phases = struct;
            phases.atten = {'cue1On', 'stimOn'};
            phases.delay = {'stimOff', 'cue3On'};
            phases.init = {'cue3On', 'prodOn'};
            phases.iti = {'prodOff', 'cue1On'};
            phases.engage = {'cue1On', 'prodOff'};
            NP.TaskBaseClass.AddEventObjects(se0, phases, 'taskTime');
            
            % Extract data
            [st, tt] = se0.GetTable('spikeTime', 'taskTime');
            for i = 1 : width(st)
                if ~iscell(st.(i))
                    st.(i) = num2cell(st.(i));
                end
            end
            st = st{:,:}';
            
            % Preallocation
            zTb = table;
            zTb.spikeTimes = st;
            nUnit = numel(st);
            nPhase = numel(phaseNames);
            nSent = numel(LMV.Param.stimIdList4);
            for p = 1 : nPhase
                pn = phaseNames{p};
                zTb.(pn+"Win") = cell(nUnit, 1+nSent);
                zTb.(pn) = NaN(nUnit, 1+nSent);
            end
            
            % Compute zeta with all trials
            for p = 1 : nPhase
                pn = phaseNames{p};
                evt = tt.(pn){1};
                zTb = testPhase(zTb, pn, evt, 1); % 1 means storing overall results in the first array column
            end
            
            function tb = testPhase(tb, vn, evt, s)
                if class(evt) == "MSessionExplorer.Event"
                    w = [evt.GetTfield('tOn') evt.GetTfield('tOff')];
                else
                    w = [evt.GetTfield('tmin') evt.GetTfield('tmax')];
                end
                for u = 1 : height(tb)
                    pv = zetatest(tb.spikeTimes{u,1}, w, median(diff(w,1,2)));
                    tb.(vn+"Win"){u,s} = w;
                    tb.(vn)(u,s) = pv;
                end
            end
            
            % Slice recording back to epochs
            tSlice = se.GetReferenceTime;
            tSlice(1) = 0;
            seEp = se0.Duplicate;
            seEp.SliceSession(tSlice, 'absolute');
            seEp.AlignTime('cue1On', 'taskTime');
            
            % Split recording by sentences
            tv = se.GetTable('taskValue');
            seEp.SetTable('taskValue', tv, 'eventValues');
            senTb = NP.SE.SplitConditions(seEp, 'conditionVars', {'taskName', 'stimText', 'stimId'});
            
            % Compute zeta by sentences
            for s = 1 : nSent
                isSen = senTb.stimId == LMV.Param.stimIdList4(s);
                if ~any(isSen)
                    fprintf("Stim #%i %s is not found in this recording\n", s, LMV.Param.stimIdList4(s));
                    continue
                end
                seSen = senTb.se(isSen);
                if seSen.numEpochs < 3
                    warning("Skip stim #%i %s due to low number (%i) of trials", s, LMV.Param.stimIdList4(s), seSen.numEpochs);
                    continue
                end
                seSen.SliceSession(0, 'absolute');
                ttSen = seSen.GetTable('taskTime');
                
                for p = 1 : nPhase
                    pn = phaseNames{p};
                    evt = ttSen.(pn){1};
                    zTb = testPhase(zTb, pn, evt, 1+s); % store result in (1+s)-th column for the specific sentence
                end
            end
            
            warning('off', 'zetatest:InsufficientSamples');
        end
        
        function uTb = PhaseEncoding(ce, uInd, maskVars)
            % Fit task phase encoding models
            % 
            %   uTb = PhaseEncoding(ce)
            %   uTb = PhaseEncoding(ce, uInd)
            %   uTb = PhaseEncoding(ce, uInd, maskVars)
            % 
            % Inputs
            %   ce          An NP.CodingExplorer object including a clusTb in its metadata.
            %   uInd        The indices of units to fit. Default include all units.
            %   maskVars    Features used as mask on top of the default masking.
            % Output
            %   uTb         
            % 
            
            % Get unit table
            NP.Unit.MergeMeta(ce);
            NP.Unit.SetUniqueClusId(ce);
            uTb = NP.Unit.GetClusTb(ce);
            uTb = NP.Unit.AddRecMeta(ce, uTb);
            
            % Configure model fitting
            mdlName = 'phase';
            if ~exist('uInd', 'var') || isempty(uInd)
                uInd = (1 : ce.numResp)';
            end
            feats = {'atten', 'stim', 'delay', 'init', 'prod', 'iti'};
            
            % Get predictors
            [~, X] = ce.GetArray('feat', [], feats);
            X(isnan(X)) = 0;
            
            % Get responses
            [~, Y] = ce.GetArray('resp', [], uInd+1);
            Y(isnan(Y)) = 0;
            
            % Masking
            m = any(X > 0, 2);
            if exist('maskVars', 'var') && ~isempty(maskVars)
                [~, M] = ce.GetArray('feat', [], maskVars);
                m = all([m M], 2);
            end
            if ~any(m)
                warning("'%s' model is not fitted because all samples are masked out", mdlName);
                uTb.(mdlName)(uInd) = cell(size(uInd(:)));
                return
            end
            X = X(m,:);
            Y = Y(m,:);
            
            % Normalize predictors and responses to be all positive and equal variance
            X = MMath.Normalize(X, 'zscore');
            X = X - min(X, [], 1);
            Y = MMath.Normalize(Y, 'zscore');
            Y = Y - min(Y, [], 1);
            
            % Print input info
            fprintf("\n%s\n", NP.SE.GetID(ce));
            fprintf("%i predictors, %i responses\n", size(X,2), size(Y,2));
            fprintf("%i raw samples, %i masked, %i valid\n", numel(m), sum(m), sum(~any(isnan([X Y]) | isinf([X Y]), 2)));
            
            % Model fitting
            fprintf("Fit %s encoding models.\n", mdlName);
            mdls = NP.Fit.BatchLsqnonneg(X, Y, 'PredictorNames', feats, 'NumShifts', 500, 'MinShift', round(0.1*sum(m)));
            uTb.(mdlName)(uInd) = mdls;
        end
        
        function sigTb = GetSigTable(tTb, varargin)
            % Convert p-values to the level of significance
            % 
            %   sigTb = GetSigTable(tTb)
            %   sigTb = GetSigTable(tTb, ..., 'VariableNames', tTb.Properties.VariableNames)
            %   sigTb = GetSigTable(tTb, ..., 'SubMask', [])
            %   sigTb = GetSigTable(tTb, ..., 'Side', 'right')
            %   sigTb = GetSigTable(tTb, ..., 'MinRespMag', 0)
            %   sigTb = GetSigTable(tTb, ..., 'Minimum', true)
            %   sigTb = GetSigTable(tTb, ..., 'Correction', false)
            %   sigTb = GetSigTable(tTb, ..., 'AlphaList', [0.05 0.01 0.001])
            % 
            % Inputs
            %   tTb             A table of p-values for different task phases.
            %   'VariableNames' Only process and return the specified columns.
            %   'SubMask'       An index or logical mask to indicate which sub-tests to include. Default is [] which 
            %                   includes all sub-tests under each variable.
            %   'Side'          Whether to look for responses that are activation ('right'), inhibition ('left'), 
            %                   or both ('both'). Default is 'right'.
            %   'MinRespMag'    Minimal response magnitude required. Default is 0 for no requirement.
            %   'Minimum'       Whether or not to take the minimum p-val over the included sub-tests. Default is true.
            %   'Correction'    Whether or not to perform Bonferroni correction after taking the minimal pval.
            %   'AlphaList'     A vector of significance levels. Default is [0.05 0.01 0.001] which correspond to 
            %                   level [1 2 3] respectively.
            % Output
            %   sigTb           A table of significance levels.
            % 
            
            p = inputParser;
            p.addParameter('VariableNames', tTb.Properties.VariableNames, @(x) isstring(x) || iscellstr(x) || ischar(x));
            p.addParameter('SubMask', [], @(x) isnumeric(x) || islogical(x));
            p.addParameter('Side', 'right', @(x) any(x == ["left", "right", "both"]));
            p.addParameter('MinRespMag', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Minimum', true, @islogical);
            p.addParameter('Correction', false, @islogical);
            p.addParameter('AlphaList', [0.05 0.01 0.001], @isnumeric);
            p.parse(varargin{:});
            varList = string(p.Results.VariableNames);
            subMask = p.Results.SubMask;
            minMag = p.Results.MinRespMag;
            side = p.Results.Side;
            isMin = p.Results.Minimum;
            isBonferroni = p.Results.Correction;
            alphas = permute(p.Results.AlphaList(:), [2 3 1]);
            
            varList = intersect(varList, tTb.Properties.VariableNames, 'stable');
            
            sigTb = table;
            for n = 1 : numel(varList)
                % Get p-values
                pn = varList(n);
                P = tTb.(pn);
                if iscell(P)
                    P = cellfun(@(x) x.p, P);
                end
                
                % Take a subset of p-values
                if isempty(subMask)
                    subMask = true(1, size(P,2));
                end
                P = P(:,subMask);
                
                % Get response magnitudes
                rn = pn+"Resp";
                if ismember(rn, tTb.Properties.VariableNames)
                    R = tTb.(rn)(:,subMask);
                else
                    R = NaN(size(P));
                end
                
                % Exclude p-values associated with subthreshold response magnitude
                P(abs(R) < minMag) = NaN;
                
                % Ignore p-values that are not for the tail of interest
                if side == "left"
                    P(R>0) = NaN;
                elseif side == "right"
                    P(R<0) = NaN;
                end
                
                % Take the smallest p-values
                if isMin
                    if isBonferroni
                        P = P * size(P,2); % bonferroni correction
                    end
                    P = min(P, [], 2);
                end
                
                % Compute significance level
                sigTb.(pn) = sum(P < alphas, 3);
            end
        end
        
        % Selectivity
        function outCell = StimSelectivityKW(respTb, phaseTb, groupTb, varargin)
            % Run Kruskal-Wallis tests to see whether units have selective responses
            % 
            %   outCell = StimSelectivityKW(respTb, phaseTb, groupTb)
            %   outCell = StimSelectivityKW(respTb, phaseTb, groupTb, ..., 'NBoot', 0)
            %   outCell = StimSelectivityKW(respTb, phaseTb, groupTb, ..., 'Mask', true(width(respTb), width(phaseTb)))
            % 
            % Inputs
            %   'NBoot'         The number of iterations to bootstrap the null distribution of the Chi-squared 
            %                   statictic. In each iteration, respnses of a unit will be translated using the 
            %                   MATLAB built-in circshift function, thus retaining the unit's stability profile 
            %                   while breaking the relationship with stimulus labels. Using circshift instead 
            %                   of permutation is robust to temporal correlation of similar consecutive trials.
            % Output
            %   outCell         A u-by-p cell array where u is the number of responses (e.g. units), p is the 
            %                   number task phases. Each element is a struct containing all three outputs (p,
            %                   tbl, stats) of the MATLAB built-in kruskalwallis function. If 'NBoot' is 
            %                   non-zero, bootstrap p-value will replace the original p-value in the 'p' field. 
            %                   One can always retrieve the original p-value from tbl{2,6}.
            % 
            % See also kruskalwallis, circshift
            
            p = inputParser;
            p.addParameter('NBoot', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Mask', true(width(respTb), width(phaseTb)), @(x) all(size(x) == [width(respTb) width(phaseTb)]));
            p.parse(varargin{:});
            nBoot = p.Results.NBoot;
            mask = logical(p.Results.Mask);
            
            nUnits = width(respTb);
            nPhases = width(phaseTb);
            nGroups = width(groupTb);
            
            % Remove instances that do not belong to any group
            if islogical(groupTb{:,:}) || width(groupTb) > 1
                isGroup = any(groupTb{:,:}, 2);
            else
                isGroup = ~ismissing(groupTb.(1));
            end
            respTb(~isGroup,:) = [];
            phaseTb(~isGroup,:) = [];
            groupTb(~isGroup,:) = [];
            
            % Convert logical group labels to string labels
            if islogical(groupTb{:,:}) || width(groupTb) > 1
                groupList = string(groupTb.Properties.VariableNames);
                gInd = sum(groupTb{:,:}.*(1:nGroups), 2);
                groupLabels = groupList(gInd);
            else
                groupLabels = groupTb.(1);
                groupList = string(unique(groupLabels));
            end
            
            % Go through task phases
            outCell = cell(nUnits, nPhases);
            for p = 1 : nPhases
                pn = phaseTb.Properties.VariableNames{p};
                isPhase = logical(phaseTb.(pn));
                fprintf("\nTest selectivity during '%s'.\n", pn);
                
                % Go through units
                for u = nUnits : -1 : 1
                    % Check if the test needs to be computed
                    s = struct;
                    if ~mask(u,p)
                        s.p = NaN;
                        s.tbl = num2cell(NaN(4,6));
                        s.stats = [];
                        outCell{u,p} = s;
                        continue
                    end
                    
                    % Collect spike rates
                    r = respTb.(u)(isPhase);
                    g = groupLabels(isPhase);
                    
                    % Remove NaNs
                    isRm = isnan(r);
                    r(isRm) = [];
                    g(isRm) = [];
                    
                    % Remove groups that have too few samples
                    gc = categorical(g);
                    [N, cats] = histcounts(gc);
                    fewTh = 3;
                    isFew = N < fewTh;
                    if any(isFew)
                        fprintf("%s: include %i/%i groups with samples > %i.\n", respTb.Properties.VariableNames{u}, ...
                            sum(~isFew), numel(groupList), fewTh-1);
                    end
                    isRm = ismember(gc, cats(isFew));
                    r(isRm) = [];
                    g(isRm) = [];
                    
                    % Run test
                    [s.p, s.tbl, s.stats] = kruskalwallis(r, g, 'off');
                    
                    % Compute null
                    if nBoot > 0
                        chiSqNull = zeros(1, nBoot);
                        parfor n = 1 : nBoot
                            rShift = circshift(r, randi(numel(r)-1));
                            [~, tbl] = kruskalwallis(rShift, g, 'off');
                            chiSqNull(n) = tbl{2,5};
                        end
                        s.p = MMath.EstimatePval(s.tbl{2,5}, chiSqNull, 'Tail', 'right');
                    end
                    
                    outCell{u,p} = s;
                end
            end
        end
        
        % Plotting
        function PlotNumberOfSigUnits(sigTbs, varargin)
            % Make bar plots for the number or fraction of statistically significant units for each task phase
            % 
            %     PlotNumberOfSigUnits(sigTbs)
            %     PlotNumberOfSigUnits(sigTbs, BarArgs)
            % 
            
            if istable(sigTbs)
                sigTbs = {sigTbs};
            end
            
            phaseNames = cellfun(@(x) x.Properties.VariableNames, sigTbs, 'Uni', false);
            phaseNames = unique(cat(2, phaseNames{:}), 'stable');
            
            nPhase = numel(phaseNames);
            nTest = numel(sigTbs);
            N = zeros(nPhase, nTest);
            for i = 1 : nTest
                isSig = sigTbs{i}{:,:} > 0;
                n = sum(isSig, 1);
                isPhase = ismember(phaseNames, sigTbs{i}.Properties.VariableNames);
                N(isPhase,i) = n;
            end
            
            ax = gca;
            bar(ax, N, varargin{:});
            ax.YLabel.String = "# of units";
            ax.XTickLabel = phaseNames;
        end
        
        function PlotWeights(ax, clusTb, uIdx, mdlName)
            % Plot weight vectors for models fitted by lasso
            % 
            %   PlotWeights(ax, clusTb, uIdx, mdlName)
            % 
            
            s = clusTb.(mdlName){uIdx};
            
            % Make weight vector
            b = s.Beta;
            x = (1:numel(b))';
            
            % Get error
            [sig, ciLow, ciHigh] = NP.Fit.GetSigLevel(s.Beta, s.null.Beta, 0.05, 'right');
            sig = logical(sig);
            
            % Plot
            plot(ax, [0.5 numel(b)+0.5], [0 0]', 'Color', [0 0 0 .5]); hold(ax, 'on')
            MPlot.ErrorShade(x, b, ciLow, ciHigh, 'IsRelative', false, 'Parent', ax);
            stem(ax, x(~sig), b(~sig), '.', 'Color', [0 0 0 .1]);
            stem(ax, x(sig), b(sig), 'o', 'Color', [0 0 1]);
            
            ax.Title.String = sprintf("u%i, %ium", clusTb.clusId(uIdx), clusTb.depth(uIdx));
%             ax.YLabel.String = "Weights";
            ax.XLim = [0 numel(b)+1];
            ax.YLim = [0 1] * max(mean(ciHigh)*3, max(b));
            ax.XTick = 1 : numel(b);
            ax.XTickLabel = s.mdl.PredictorNames;
            ax.Box = 'off';
        end
        
        function PlotBenchmark(scoreTb, testTbs, taskPhases)
            % Bar plots comparing the flase positive and false negative rates of different testing methods
            % 
            %   PlotBenchmark(scoreTb, testTbs, taskPhases)
            % 
            
            err1 = zeros(numel(taskPhases), numel(testTbs));
            err2 = err1;
            
            for j = 1 : numel(testTbs)
                for i = 1 : numel(taskPhases)
                    pn = taskPhases(i);
                    y = scoreTb.(pn) > 0;
                    yHat = testTbs{j}.(pn) > 0;
                    cm = confusionmat(y, yHat);
                    cm = cm ./ sum(cm,2);
                    err1(i,j) = cm(1,2);
                    err2(i,j) = cm(2,1);
                end
            end
            
            cc = flag;
            cc = cc([1:3 5:7 9:11],:);
            
            ax = nexttile;
            b = bar(err2);
            for i = 1 : numel(b)
                b(i).FaceColor = cc(i,:);
            end
            ax.XTickLabel = taskPhases;
            ax.YLabel.String = "Frac. of GT positive";
            ax.Title.String = "False negative";
            MPlot.Axes(ax);
            
            ax = nexttile;
            b = bar(err1);
            for i = 1 : numel(b)
                b(i).FaceColor = cc(i,:);
            end
            ax.XTickLabel = taskPhases;
            ax.YLabel.String = "Frac. of GT negative";
            ax.Title.String = "False positive";
            MPlot.Axes(ax);
        end
        
    end
    
end
