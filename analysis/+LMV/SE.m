classdef SE
    
    methods(Static)
        % General
        function varargout = Transform(se, stages, ops)
            % Transform se via desired processing stages
            % 
            %   [se, senTb] = Transform(se)
            %   [se, senTb] = Transform(se, stages)
            %   [se, senTb] = Transform(se, stages, ops)
            % 
            % Inputs
            %   se          MSessionExplorer object.
            %   stages      Any combination of {'enrich', 'morph', 'sentence'}.
            %   ops         A struct customized from ops = NP.Param.Enrich
            % Outputs
            %   se          Transformed MSessionExplorer object.
            %   senTb       A table of split se objects by sentences.
            % 
            
            if nargin < 2 || isempty(stages)
                stages = '';
            end
            stages = cellstr(stages);
            
            se = se.Duplicate;
            
            % Some standardization
            NP.Unit.MergeMeta(se); % combine ksMeta
            NP.Unit.SetUniqueClusId(se);
            LMV.SE.AddTaskValueTable(se); % recompute taskValue table to apply latest update
            
            % Remove outdated data
            var2rm = "pitch" + ["Min", "Max", "Up", "Down"];
            tt = se.GetTable('taskTime');
            isRm = ismember(tt.Properties.VariableNames, var2rm);
            tt(:,isRm) = [];
            se.SetTable('taskTime', tt);
            
            % Enrich se
            if ismember('enrich', stages)
                if nargin < 3
                    ops = NP.Param.Enrich;
                    ops.isFiltSpeech = true;
                    ops.isMel = true;
                    ops.isPitch = true;
                    ops.isArtic = true;
                    ops.isSpkRate = true;
                end
                NP.SE.Enrich(se, ops);
            end
            
            % Get LMV epochs
            taskTb = se.SplitConditions('taskName', 'taskValue');
            se = taskTb.se(taskTb.taskName=="lmv");
            varargout{1} = se;
            
            if ismember('morph', stages)
                % Time morphing
                NP.SE.SetMorphTimes(se);
                se = NP.SE.MorphSession(se);
                
                % Make sure all spike times are in cell array
                st = se.GetTable('spikeTime');
                for i = 1 : width(st)
                    if ~iscell(st.(i))
                        st.(i) = num2cell(st.(i));
                    end
                end
                se.SetTable('spikeTime', st);
                
                varargout{1} = se;
            end
            
            if ismember('sentence', stages)
                % Reslice epochs and split se by sentence
                senTb = LMV.SE.SplitBySentence(se);
                varargout{1} = se;
                varargout{2} = senTb;
            end
        end
        
        % Features
        function tv = AddTaskValueTable(se)
            % Make taskValue table
            
            % Continue with existing table or create new
            if ismember('taskValue', se.tableNames)
                tv = se.GetTable('taskValue');
            else
                tv = table;
            end
            tt = se.GetTable('taskTime');
            tv.trialNum = (1:height(tt))';
            tv.recId(:) = string(se.userData.expInfo.recId);
            se.SetTable('taskValue', tv, 'eventValues');
            
            % Set stim ID in taskValue
            NP.TaskBaseClass.MigrateStimId(se, '^.{5}_si.{3,}');
            NP.TaskBaseClass.MigrateStimId(se, '^BR_');
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            
            % Check required columns in taskTime table
            requiredCols = {'stim', 'prod'};
            assert(all(ismember(requiredCols, tt.Properties.VariableNames)), ...
                "LMV task requires the '%s' and '%s' columns in the taskTime table.", requiredCols{:});
            
            numTrials = 0;
            for i = 1 : height(tt)
                % Get task name
                eStim = tt.stim(i);
                eProd = tt.prod{i};
                if isnan(eStim)
                    continue
                end
                taskName = lower(eStim(1).GetVfield('task'));
                if taskName ~= "lmv"
                    continue
                end
                numTrials = numTrials + 1;
                
                % Get stim and prod texts
                dtype = class(eStim);
                if dtype == "NP.TGEvent"
                    stimText = eStim(1).GetParentLabel;
                    prodText = eProd.GetAllParentLabel();
                elseif dtype == "MSessionExplorer.Event"
                    stimText = eStim(1).GetVfield('type');
                    prodText = string(NaN);
                else
                    error("Unrecognized event data type '%s'", dtype);
                end
                
                tv.taskName(i) = taskName;
                tv.stimText(i) = stimText;
                tv.prodText(i) = prodText;
            end
            
            if ~numTrials
                fprintf("This se does not have any LMV trial\n");
                return
            end
            
            % Make unique numeric stim IDs
            m = tv.taskName == "lmv";
            [~, ~, tv.stimNumId(m)] = unique(tv.stimText(m), 'stable');
            
            % Add taskValues table to se
            se.SetTable('taskValue', tv);
        end
        
        function AddTrialEvents(se)
            % Add a eventTimes table named 'trial' to se. Each column is a trial type named after the sentence ID. 
            % Each event object is located at trialOn and spans from this trialOn (tOn) to the next trialOn (tOff).
            % 
            %   AddTrialEvents(se)
            % 
            
            % 
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            rt = se.GetReferenceTime('taskTime');
            trialOn = tt.trialOn + rt;
            trialOn = [trialOn; trialOn(end)+10];
            
            % Make a single epoch table of sentence events
            tr = table;
            for k = 1 : numel(LMV.Param.stimIdList)
                id = LMV.Param.stimIdList(k);
                ind = find(tv.stimId == id);
                evts = MSessionExplorer.Event(trialOn(ind));
                if ~isempty(ind) % recording may not have this sentence
                    evts = evts.SetTfield('tOn', trialOn(ind));
                    evts = evts.SetTfield('tOff', trialOn(ind+1));
                end
                tr.(id){1} = evts;
            end
            
            % Slice to trials
            seTemp = MSessionExplorer;
            seTemp.isVerbose = false;
            seTemp.SetTable('trial', tr, 'eventTimes', 0);
            seTemp.SliceSession(rt, 'absolute');
            tr = seTemp.GetTable('trial');
            
            se.SetTable('trial', tr, 'eventTimes', rt);
            
        end
        
        function AddSentenceOnsets(se)
            % Add sentOn including both stimOn and prodOn to taskTime table
            % 
            %   AddSentenceOnsets(se)
            % 
            tt = se.GetTable('taskTime');
            stimOn = tt.stimOn;
            if ~iscell(stimOn)
                stimOn = num2cell(stimOn);
            end
            prodOn = tt.prodOn;
            if ~iscell(prodOn)
                prodOn = num2cell(prodOn);
            end
            sentOn = cellfun(@(x,y) sort([x;y]), stimOn, prodOn, 'Uni', false);
            sentOn = cellfun(@(x) x(~isnan(x)), sentOn, 'Uni', false);
            se.SetColumn('taskTime', 'sentOn', sentOn);
        end
        
        % Morphing
        function [tt, tv] = MatchProd2Stim(tt, tv, rt)
            % Match the transcript of speech production to stimulus
            % 
            %   [tt, tv] = MatchProd2Stim(tt, tv)
            %   [tt, tv] = MatchProd2Stim(tt, tv, rt)
            % 
            
            % Check required columns in taskTime table
            requiredCols = {'cue3On', 'prod'};
            assert(all(ismember(requiredCols, tt.Properties.VariableNames)), ...
                "MatchProd2Stim requires the '%s' and '%s' column in the taskTime table.", requiredCols{:});
            
            requiredCols = {'stimText'};
            assert(all(ismember(requiredCols, tv.Properties.VariableNames)), ...
                "MatchProd2Stim requires the '%s' column in the taskValue table.", requiredCols{:});
            
            for i = 1 : height(tt)
                % Get stim phone labels
                stimPh = Cut(Cut(tt.stim(i)));
                s1 = stimPh.GetParentLabel();
                
                % Get production phone objects (those within each trial)
                prodPh = Cut(Cut(tt.prod{i}));
                prodPh = prodPh(prodPh > tt.cue3On(i)); % remove events before cue3On or are NaN
                if exist('rt', 'var') && ~isempty(rt) && i < height(tt)
                    prodPh = prodPh(prodPh.GetTfield('tmax')+rt(i) < tt.cue1On(i+1)+rt(i+1)); % remove events that bleed into the next trial
                end
                
                if isempty(prodPh)
                    % Fill with placeholders
                    tv.alignment{i} = [];
                    tv.alignScore(i) = -Inf;
                    tv.tReact(i) = NaN;
                    tt.stimMatchOn(i) = NaN;
                    tt.stimMatchOff(i) = NaN;
                    tt.prodMatchOn(i) = NaN;
                    tt.prodMatchOff(i) = NaN;
                    tt.matchOn(i) = NaN; % kept for backward compatibility
                    tt.matchOff(i) = NaN; % kept for backward compatibility
                    continue
                end
                
                % Align sequences
                s2 = prodPh.GetParentLabel();
                [k1, k2, info] = MLing.FindAlignedTokens(s1, s2);
                aa = info.aa;
                sc = info.tscore / numel(s1); % normalize by maximal edit distance
                
                tv.alignment{i} = string(aa);
                tv.alignScore(i) = sc;
                tv.tReact(i) = prodPh(k2(1)).t - tt.cue3On(i); % get reaction time using the first prod event that matches stim
                tt.stimMatchOn(i) = stimPh(k1(1)).T.tmin;
                tt.stimMatchOff(i) = stimPh(k1(end)).T.tmax;
                tt.prodMatchOn(i) = prodPh(k2(1)).T.tmin;
                tt.prodMatchOff(i) = prodPh(k2(end)).T.tmax;
                tt.matchOn(i) = tt.prodMatchOn(i); % kept for backward compatibility
                tt.matchOff(i) = tt.prodMatchOff(i); % kept for backward compatibility
            end
        end
        
        function [tt, tv] = FindMorphTimes(tt, tv, rt, targetType)
            % Find matching key times between source and temaplte trials
            % 
            %   [tt, tv] = FindMorphTimes(tt, tv)
            %   [tt, tv] = FindMorphTimes(tt, tv, rt)
            %   [tt, tv] = FindMorphTimes(tt, tv, rt, targetType)
            % 
            
            if nargin < 4
                targetType = 'median';
            end
            targetType = lower(targetType);
            
            % Preallocate columns
            tt.morphFrom = num2cell(tt.cue1On);
            tt.morphTo = num2cell(tt.cue1On);
            tv.tempTrialNum = NaN(size(tt.cue1On));
            
            % Cut sentences to phones
            for i = height(tt) : -1 : 1
                % Cut stim
                stimPh{i,1} = Cut(Cut(tt.stim(i)));
                
                % Cut prod
                tge = Cut(Cut(tt.prod{i}));
                tge = tge(tge > tt.cue3On(i)); % remove events before cue3On or are NaN
                if exist('rt', 'var') && ~isempty(rt) && i < height(tt)
                    tge = tge(tge+rt(i) < tt.cue1On(i+1)+rt(i+1)); % remove events that bleed into the next trial
                end
                prodPh{i,1} = tge;
            end
            
            % Find unique stimuli
            [stimList, ~, trialStimId] = unique(tv.stimText);
            
            % Find morph times stim by stim
            for s = 1 : numel(stimList)
                % Find non-empty trials of the same sentence
                trialInd = find(trialStimId==s);
                
                % Find the trial(s) that has the highest alignment score
                sc = tv.alignScore(trialInd);
                bestTrials = trialInd(sc == max(sc));
                
                % Break ties by choosing the trial with reaction time closest to median
                tReact = tv.tReact(bestTrials);
                [~, I] = min(abs(tReact - median(tReact)));
                bestTrial = bestTrials(I);
                
                % Get the template speech events
                stimTp = stimPh{bestTrial};
                prodTp = prodPh{bestTrial};
                assert(~isempty(prodTp), "Template for '%s' is empty", stimList(s));
                
                % Construct time matrix
                preStimEvts = {'cue1On', 'cue1Off'}; % events preceeding stim
                if ismember('cue2On', tt.Properties.VariableNames) && ~any(isnan(tt.cue2On))
                    periDelayEvts = {'stimOff', 'cue2On', 'cue2Off', 'cue3On'}; % events preceeding repetition
                else
                    periDelayEvts = {'stimOff', 'cue3On'}; % events preceeding repetition
                end
                nPreProd = numel(preStimEvts) + numel(stimTp) + numel(periDelayEvts);
                nProdOff = 1;
                T = NaN(numel(trialInd), nPreProd+numel(prodTp)+nProdOff);
                
                % Collect source event times to align from
                for j = 1 : numel(trialInd)
                    % Get the source event times or objects
                    k = trialInd(j);
                    tPreProd = [tt{k,preStimEvts}'; double(stimPh{k}); tt{k,periDelayEvts}'];
                    prodSrc = prodPh{k};
                    
                    % Store source event times
                    if isempty(prodSrc)
                        % No prod, so only use stim and task times
                        tSrc = tPreProd;
                        kTp = 1:nPreProd;
                        info = [];
                    else
                        % Find timestamps of matching prod words
                        [~, tProdSrc, info] = NP.TaskBaseClass.FindMatchedSpeechTimes(prodTp, prodSrc);
                        info.nscore = info.tscore / numel(prodTp); % normalize alignment score to max score
                        
                        % Append prodOff time if the template is being matched to the end
                        if info.kk(1,end) == numel(prodTp)
                            prodOff = prodSrc(info.kk(2,end)).T.tmax;
                        else
                            prodOff = [];
                        end
                        
                        % Store full source event times
                        tSrc = [tPreProd; tProdSrc; prodOff];
                        kTp = [1:nPreProd nPreProd+info.kk(1,:) nPreProd+info.kk(1,end)+ones(size(prodOff))];
                    end
                    T(j,kTp) = tSrc;
                    
                    tv.alignInfo{k} = info;
                    tv.tempTrialNum(k) = tv.trialNum(bestTrial);
                    tt.morphFrom{k} = tSrc;
                end
                
                % Determine template times
                itvl = diff(T, 1, 2);
                if targetType == "median"
                    % Use median event intervals
                    medItvl = median(itvl, 1, 'omitnan');
                    medItvl(nPreProd) = 1; % force the mental holding interval to the designed mean of 1 sec
                    tTp = cumsum([tt.cue1On(bestTrial) medItvl])';
                elseif targetType == "best"
                    % Use the best aligned trial
                    bestRep = trialInd == bestTrial;
                    tTp = cumsum([tt.cue1On(bestTrial) itvl(bestRep,:)])';
                end
                
                % Find and save matching times from the template
                for j = 1 : numel(trialInd)
                    k = trialInd(j);
%                     info = tv.alignInfo{k};
%                     if isempty(info)
%                         continue
%                     end
%                     kTp = [1:nPreProd nPreProd+info.kk(1,:) nPreProd+info.kk(1,end)+1];
                    kTp = ~isnan(T(j,:));
                    tt.morphTo{k} = tTp(kTp);
                end
            end
        end
        
        function isBad = IsBadTrials(se, stimTh, repTh)
            % Return logical indices of trials with bad performance
            % Performance is determined by the two alignment scores in the taskValue table
            % 
            %   isBad = IsBadTrials(se, stimTh, repTh)
            % 
            tv = se.GetTable('taskValue');
            for i = height(tv) : -1 : 1
                stimScore(i,1) = tv.alignScore(i);
                if ~isempty(tv.alignInfo{i})
                    repScore(i,1) = tv.alignInfo{i}.nscore;
                else
                    repScore(i,1) = -Inf;
                end
            end
            isBad = stimScore < stimTh | repScore < repTh;
        end
        
        % PETH
        function senTb = SplitBySentence(se)
            % Split se by unique sentences after reslicing and screening trials
            
            % Slice trials to include a period before and after
            tOffset = [-1 1]; % in sec
            se = NP.SE.BleedTrials(se, tOffset);
            
            % Remove trials with bad performance
            isBad = LMV.SE.IsBadTrials(se, -Inf, 0.5);
            se.RemoveEpochs(isBad);
            
            % Split
            senTb = NP.TaskBaseClass.SplitBySentence(se);
        end
        
        function ce = ComputeSessionPETH(se, varargin)
            % Compute PETHs from the given se
            % 
            %   ce = ComputeSessionPETH(se)
            %   ce = ComputeSessionPETH(..., 'BoundaryEvents', ["trialOn", "prodMatchOff"])
            %   ce = ComputeSessionPETH(..., 'StatFun', @MMath.MeanStats)
            %   ce = ComputeSessionPETH(..., 'Modification', "none")
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('BoundaryEvents', ["trialOn", "prodMatchOff"], @(x) iscellstr(x) || isstring(x));
            p.addParameter('Modification', "none", @(x) ismember(lower(x), ["sem_discounted", "none"]));
            p.parse(varargin{:});
            bEvts = string(p.Results.BoundaryEvents);
            mOpt = string(p.Results.Modification);
            argStruct = p.Unmatched;
            
            % Derive resampling parameters
            tt = se.GetTable('taskTime');
            tOn = median(tt.(bEvts(1)), 'omitnan');
            tOff = median(tt.(bEvts(2)), 'omitnan');
            ops = NP.Param.Resample();
            ops.rsWin = round([tOn tOff] + [-1 1]*0.5, 6); % round off eps introduced during morphing
            ops.rsBinSize = 0.01;
            
            % Resample response timeseries
            rTb = NP.SE.ResampleResponses(se, ops);
            
            % Add recording info to clusTb
            ud = se.userData;
            for r = 1 : numel(ud)
                for p = 1 : numel(ud(r).ksMeta)
                    ud(r).ksMeta(p).clusTb = NP.Unit.AddRecMeta(ud(r), ud(r).ksMeta(p).clusTb);
                end
            end
            se.userData = ud;
            
            % Compute intra-recording average with trial masks
            cTb = NP.Unit.GetClusTb(se);
%             tv = se.GetTable('taskValue');
%             trialMask = tv.recId == cTb.recId';
%             mrTb = NP.SE.MeanTimeseries(rTb, 'TrialMask', trialMask, argStruct);
            mrTb = NP.SE.MeanTimeseries(rTb, argStruct);
            
            % Modify PETHs
            switch mOpt
                case 'sem_discounted'
                    mrTb{1,2:end} = cellfun(@(x,y) max(x-y, 0), mrTb{1,2:end}, mrTb{3,2:end}, 'Uni', false);
            end
            
            % Compute PETHs
            ce = NP.CodingExplorer();
            ce.SetTable('resp', mrTb(1,:), 'timeSeries');
            ce.SetTable('sd', mrTb(2,:), 'timeSeries');
            ce.SetTable('sem', mrTb(3,:), 'timeSeries');
            
            % Add PETH features
            stats = NP.PETH.ComputePatternStats(ce, [], true);
            cTb.peakSpkCount = stats.pkVal * ops.rsBinSize;
            cTb.peakSpkRate = stats.pkVal;
            cTb.mi = stats.mi;
            
            % Add metadata
            ce.userData = se.userData;
            ce.userData.rsOps = ops;
            ce.clusTb = cTb;
            
            % Add tables of an example trial to pe
            LMV.SE.AddTemplateTrial(ce, se);
        end
        
        function [ce, senTb] = ComputeSentencePETH(se, varargin)
            % Compute PETHs for each sentence
            % 
            %   [ce, senTb] = ComputeSentencePETH(se)
            %   [ce, senTb] = ComputeSentencePETH(senTb)
            %   [ce, senTb] = ComputeSentencePETH(..., 'BoundaryEvents', ["trialOn", "prodMatchOff"])
            %   [ce, senTb] = ComputeSentencePETH(..., 'StatFun', @MMath.MeanStats)
            %   [ce, senTb] = ComputeSentencePETH(..., 'Modification', "none")
            % 
            
            if istable(se)
                senTb = se;
            else
                senTb = LMV.SE.SplitBySentence(se);
            end
            
            % Run through sentences
            for s = 1 : height(senTb)
                fprintf("\n%s\n", senTb.stimText(s));
                senTb.ce(s) = LMV.SE.ComputeSessionPETH(senTb.se(s), varargin{:});
            end
            
            % Combine ce of different sentecnes
            ce = Merge(senTb.ce);
            ce.userData = senTb.ce(1).userData;
            ce.SetColumn('taskValue', 'numTrial', senTb.numTrial);
            
            % Compute modulation index (MI)
            clusTb = ce.clusTb;
            [resp, sem] = ce.GetTable('resp', 'sem');
            resp.time = [];
            sem.time = [];
            clusTb.mi = NP.PETH.ComputeModulationIndex(resp, sem)';
            
            % Use MI as a heuristic to find units with salient and reliable patterns
            [~, ind] = sort(senTb.numTrial, 'descend');
            ind = ind(1:4); % the four most repeated
            clusTb.mi4 = mean(clusTb.mi(:,ind), 2);
            ce.clusTb = clusTb;
        end
        
        function AddTemplateTrial(seTo, seFrom, varargin)
            % Add single-trial data tables from one se to another
            % 
            %   AddTemplateTrial(seTo, seFrom)
            % 
            
            p = inputParser;
            p.addParameter('IncludeTables', [], @(x) all(ismember(x, seFrom.tableNames)));
            p.parse(varargin{:});
            tbInclude = p.Results.IncludeTables;
            
            % Find mask for template trials
            tv = seFrom.GetTable('taskValue');
            mTemp = ismember(tv.trialNum, tv.tempTrialNum);
            mTemp = find(mTemp, seTo.numEpochs);
            if isempty(mTemp)
                mTemp = 1;
                fprintf("No template trial found. Use the first trial.\n");
            end
            
            % Set tables
            for i = 1 : numel(seFrom.tableNames)
                tn = seFrom.tableNames{i};
                if isempty(tbInclude) || ismember(tn, tbInclude)
                    seTo.SetTable(tn, seFrom.tot.tableData{i}(mTemp,:), seFrom.tot.tableType{i});
                end
            end
        end
        
        function StandardizeDataType(seArray)
            % Make sure all variables are in standard data containers
            % 'prod', 'prodOn', 'prodOff' from 'taskTime' table should be cell arrays.
            % 
            %   StandardizeDataType(seArray)
            % 
            cellVars = {'prod', 'prodOn', 'prodOff'};
            for i = 1 : numel(seArray)
                tt = seArray(i).GetTable('taskTime');
                for j = 1 : numel(cellVars)
                    vn = cellVars{j};
                    if ~iscell(tt.(vn))
                        fprintf("%i - %s: convert %s to cell array\n", i, NP.SE.GetID(seArray(i)), vn);
                        tt.(vn) = num2cell(tt.(vn));
                    end
                end
                seArray(i).SetTable('taskTime', tt);
            end
        end
        
        % 
        function se = UnitCache2SE(sUnit, varargin)
            % Construct a se object from unit cache
            % 
            %   se = UnitCache2SE(sUnit)
            %   se = UnitCache2SE(sUnit, ..., 'StimIdList', LMV.Param.stimIdList14)
            %   se = UnitCache2SE(sUnit, ..., 'NumTrials', "auto")
            % 
            
            % Process inputs
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('StimIdList', LMV.Param.stimIdList14, @isstring);
            p.addParameter('NumTrials', "auto", @(x) any(x==["max", "min", "median", "mode"]));
            p.parse(varargin{:});
            stimIdList = p.Results.StimIdList;
            nTrialMode = p.Results.NumTrials;
            
            % Find the maximal numbers of repeat for each stim
            nUnit = numel(sUnit);
            nStim = numel(stimIdList);
            nRep = zeros(nStim, nUnit);
            for i = 1 : nUnit
                s = sUnit{i};
                for j = 1 : nStim
                    m = s.tv.stimId == stimIdList(j);
                    if ~any(m)
                        continue
                    end
                    nRep(j,i) = numel(s.st{m});
                end
            end
            switch nTrialMode
                case 'max'
                    nRep = max(nRep, [], 2);
                case 'min'
                    nRep = min(nRep, [], 2);
                case 'median'
                    nRep = median(nRep, 2);
                case 'mode'
                    nRep = mode(nRep, 2);
                case 'auto'
                    [M, F] = mode(nRep, 2);
                    M(F==1) = median(nRep(F==1,:), 2);
                    nRep = M;
            end
            nRep = round(nRep); % avoid 0.5 from median with even samples
            isMissing = nRep==0;
            nRep(isMissing) = [];
            stimIdList(isMissing) = [];
            nStim = numel(stimIdList);
            
            % Make taskValue and taskTime table
            tvCat = CatTables(cellfun(@(x) x.tv, sUnit, 'Uni', false));
            ttCat = CatTables(cellfun(@(x) x.tt, sUnit, 'Uni', false));
            
            function tbCat = CatTables(tbs)
                % Concatenate tables with common variables (currently not handling mismatched data types)
                vars = cellfun(@(x) x.Properties.VariableNames, tbs, 'Uni', false);
                vars = cat(2, vars{:});
                vars = unique(vars, 'stable');
                for n = 1 : numel(tbs)
                    vars(~ismember(vars, tbs{n}.Properties.VariableNames)) = [];
                end
                tbCat = cellfun(@(x) x(:,vars), tbs, 'Uni', false);
                tbCat = cat(1, tbCat{:});
            end
            
            [~, I] = MMath.SortLike(tvCat.stimId, stimIdList, false);
            I = repelem(I, nRep);
            tv = tvCat(I,:);
            tt = ttCat(I,:);
            
            % Make spike times table
            stCell = num2cell(NaN(height(tt), nUnit));
            for i = 1 : nUnit
                s = sUnit{i};
                for j = 1 : nStim
                    % Check if the stim is present
                    m = s.tv.stimId == stimIdList(j);
                    if ~any(m)
                        continue
                    end
                    
                    % Determine the number of trials to include
                    t = s.st{m};
                    nTrials = min(nRep(j), numel(t));
                    
                    % Find the range of rows
                    a = sum(nRep(1:j-1)) + 1;
                    b = a-1 + nTrials;
                    
                    % Evenly sample trials (when nTrials < numel(t))
                    ind = round(linspace(1, numel(t), nTrials));
                    
                    stCell(a:b,i) = t(ind);
                end
            end
            uNames = "u" + cellfun(@(x) x.unitInfo.clusId, sUnit);
            st = cell2table(stCell, 'VariableNames', uNames);
            
            % Make clusTb
            clusTb = CatTables(cellfun(@(x) struct2table(x.unitInfo, 'AsArray', true), sUnit, 'Uni', false));
            
            % Construct se object
            se = MSessionExplorer;
            se.userData.ksMeta.clusTb = clusTb;
            se.SetTable('taskTime', tt, 'eventTimes');
            se.SetTable('taskValue', tv, 'eventValues');
            se.SetTable('spikeTime', st, 'eventTimes');
        end
        
        % Misc
        function out = ToOutStruct(se)
            % Convert se to out
            %
            %   out = ToOutStruct(se)
            % 
            
            % Combine ksMeta
            NP.Unit.MergeMeta(se);
            NP.Unit.SetUniqueClusId(se);
            
            % Standard enrichment
            ops = NP.Param.Enrich;
            ops.isSpkRate = true;
            ops.isMel = true;
            ops.isPitch = true;
            ops.isArtic = true;
            ops.spkBinSize = 0.0025; % set this to ≤0.0015 to have a maximum of one spike (i.e. 1/spkBinSize Hz) per bin
            ops.spkKerSize = 0.015; % set this to zero for no Gaussian smoothing
            NP.SE.Enrich(se, ops);
            
            % Vectorize se
            se0 = se;
            se = se0.Duplicate;
            se.SliceSession(0, 'absolute');
            
            % Add phoneme events
            phones = [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.consonants]';
            phGroups = { ...
                'high',         NP.Phone.high; ...
                'mid',          NP.Phone.mid; ...
                'low',          NP.Phone.low; ...
                'front',        NP.Phone.front; ...
                'back',         NP.Phone.back; ...
                'rounded',      NP.Phone.rounded; ...
                'plosives',     NP.Phone.plosives; ...
                'fricatives',   NP.Phone.fricatives; ...
                'nasals',       NP.Phone.nasals; ...
                'approximants', NP.Phone.approximants; ...
                'labial',       NP.Phone.labial; ...
                'velar',        NP.Phone.velar; ...
                'coronal',      NP.Phone.coronal; ...
                'glottal',      NP.Phone.glottal; ...
                'dental',       NP.Phone.dental; ...
                };
            phEvents = [repmat(phones, [1 2]); phGroups];
            NP.Phone.AddPhoneEvents(se, phEvents);
            
            
            % Initialize option struct from the enrichment option struct
            ops = NP.Param.Resample(ops);
            
            % Set the resmapling window and sampling rate
            tRef = se0.GetReferenceTime;
            ops.rsWin = tRef([1 end])' + [-5 15];
            ops.rsBinSize = 1e-3; % resample at 1000Hz
            
            % Specify the features to resample
            sCell = { ...
                'taskTime', {'cue1On', 'cue1Off', 'cue3On', 'cue3Off', 'stimOn', 'prodOn'}, []; ...
                'mel',      {'mic', 'speaker1'}, []; ...
                'ni',       {'mic', 'speaker1'}, []; ...
                'inten',    {'env', 'peakEnv', 'peakRate'}, []; ...
                'pitch',    {'rF0', 'drF0', 'brF0', 'voicing', 'phrase', 'accent'}, []; ...
                'artic',    {'ja', 'la', 'pro', 'ttcc', 'tbcc', 'tdcc', 'ttcl', 'tbcl', 'v_x', 'v_y'}, []; ...
                'phone',    phones', []; ...
                };
            ops.rsVars = struct('tableName', sCell(:,1), 'varNames', sCell(:,2), 'readout', sCell(:,3));
            
            % Resample features and responses
            fTb = NP.SE.ResampleFeatures(se, ops);
            rTb = NP.SE.ResampleResponses(se, ops);
            
            % Construct ce
            ce = NP.CodingExplorer;
            ce.userData = se.userData;
            ce.userData.ops = ops;
            ce.SetTable('feat', fTb, 'timeSeries', 0);
            ce.SetTable('resp', rTb, 'timeSeries', 0);
            
            % Slice ce back to trials
            tSlice = tRef;
            tSlice(1) = 0;
            ce.SliceSession(tSlice, 'absolute');
            
            % Align time to trial onsets
            tAlign = se0.GetColumn('taskTime', 'stimOn');
            tAlign(1) = tRef(1) + tAlign(1);
            ce.AlignTime(tAlign);
            ce.SetTable('taskValue', se0.GetTable('taskValue'), 'eventValues');
            
            
            % Unload data
            [fTb, ~, tv] = ce.GetTable(ce.tableNames{:});
            tv.stimId(ismissing(tv.stimId)) = "";
            ops = ce.userData.ops;
            
            % Find time windows
            tmin = -0.5;
            tmax = 1;
            epochOn = cellfun(@(x) x(1), fTb.time) + tmin;
            epochOff = cellfun(@(x) x(end), fTb.time) + tmax;
            epochWin = [epochOn epochOff];
            
            % Reslice tables
            fTb = ce.SliceTimeSeries('feat', epochWin, 'Fill', 'bleed');
            rTb = ce.SliceTimeSeries('resp', epochWin, 'Fill', 'bleed');
            
            % Select units
%             recId = NP.SE.GetID(ce);
%             [~, uId] = LMV.Param.GetSelectedClusId(recId);
%             uInd = ismember(ce.clusTb.clusId, uId);
%             uInd = find(uInd);
            uInd = 1:height(ce.clusTb);
            
            % Construct out struct from ce
            for i = ce.numEpochs : -1 : 1
                s = struct;
                s.name = char(tv.stimId(i));
                s.tmin = tmin;
                s.tmax = tmax;
                s.spike_fs = 1/ops.rsBinSize;
                s.chs_original = [];
                s.chs_physical = [];
                s.duration = diff(fTb.time{i}([1 end]));
                s.soundf = 1/ops.rsBinSize;
                s.sound = fTb.mic_ni{i};
                s.aud = fTb.mic{i}';
                s.spikes_sua = cell2mat(rTb{i,uInd+1})';
                s.spike_wave = [];
                s.anin = fTb.speaker1_ni{i}';
                
                varNames = [{'cue1On', 'cue1Off', 'cue3On', 'cue3Off', 'stimOn', 'prodOn'} phones'];
                for n = 1 : numel(varNames)
                    vn = varNames{n};
                    s.(vn) = fTb.(vn){i};
                end
                
                if i == 1
                    sClus = table2struct(ce.clusTb(uInd,:));
                    for u = 1 : numel(sClus)
                        sClus(u).Amplitude = NaN;
                        sClus(u).ContamPct = NaN;
                        sClus(u).KSLabel = NaN;
                        sClus(u).amp = NaN;
                        sClus(u).ch = NaN;
                        sClus(u).id = sClus(u).clusId;
                        sClus(u).fr = NaN;
                        sClus(u).n_spikes = sClus(u).numSpikes;
                        sClus(u).sua_id = NaN;
                    end
                    s.sua_info = sClus;
                else
                    s.sua_info = [];
                end
                
                out(i) = s;
            end
        end
        
        function out = ToOutStruct2(se)
            % Convert se to out
            %
            %   out = ToOutStruct(se)
            % 
            
            % Vectorize se
            se0 = se;
            se = se0.Duplicate;
            se.SliceSession(0, 'absolute');
            
            % Add phoneme events
            phones = [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.consonants]';
            phGroups = { ...
                'high',         NP.Phone.high; ...
                'mid',          NP.Phone.mid; ...
                'low',          NP.Phone.low; ...
                'front',        NP.Phone.front; ...
                'back',         NP.Phone.back; ...
                'rounded',      NP.Phone.rounded; ...
                'plosives',     NP.Phone.plosives; ...
                'fricatives',   NP.Phone.fricatives; ...
                'nasals',       NP.Phone.nasals; ...
                'approximants', NP.Phone.approximants; ...
                'labial',       NP.Phone.labial; ...
                'velar',        NP.Phone.velar; ...
                'coronal',      NP.Phone.coronal; ...
                'glottal',      NP.Phone.glottal; ...
                'dental',       NP.Phone.dental; ...
                };
            phEvents = [repmat(phones, [1 2]); phGroups];
            NP.Phone.AddPhoneEvents(se, phEvents);
            
            
            % Initialize option struct from the enrichment option struct
            ops = NP.Param.Resample();
            
            % Set the resmapling window and sampling rate
            tRef = se0.GetReferenceTime;
            ops.rsWin = tRef([1 end])' + [-5 15];
            ops.rsBinSize = 1e-3; % resample at 1000Hz
            
            % Specify the features to resample
            sCell = { ...
                'taskTime', {'cue1On', 'cue1Off', 'cue3On', 'cue3Off', 'stimOn', 'prodOn'}, []; ...
                'mel',      {'mic', 'speaker1'}, []; ...
                'ni',       {'mic', 'speaker1'}, []; ...
                'inten',    {'env', 'peakEnv', 'peakRate'}, []; ...
                'pitch',    {'rF0', 'drF0'}, []; ...
                'artic',    {'ja', 'la', 'pro', 'ttcc', 'tbcc', 'tdcc', 'ttcl', 'tbcl', 'v_x', 'v_y'}, []; ...
                'phone',    phones', []; ...
                };
            ops.rsVars = struct('tableName', sCell(:,1), 'varNames', sCell(:,2), 'readout', sCell(:,3));
            
            % Resample features and responses
            fTb = NP.SE.ResampleFeatures(se, ops);
            rTb = NP.SE.ResampleResponses(se, ops);
            
            % Construct ce
            ce = NP.CodingExplorer;
            ce.userData = se.userData;
            ce.userData.ops = ops;
            ce.SetTable('feat', fTb, 'timeSeries', 0);
            ce.SetTable('resp', rTb, 'timeSeries', 0);
            
            % Slice ce back to trials
            tSlice = tRef;
            tSlice(1) = 0;
            ce.SliceSession(tSlice, 'absolute');
            
            % Align time to trial onsets
            tAlign = se0.GetColumn('taskTime', 'stimOn');
            tAlign(1) = tRef(1) + tAlign(1);
            ce.AlignTime(tAlign);
            ce.SetTable('taskValue', se0.GetTable('taskValue'), 'eventValues');
            
            
            % Unload data
            [fTb, ~, tv] = ce.GetTable(ce.tableNames{:});
            tv.stimId(ismissing(tv.stimId)) = "";
            ops = ce.userData.ops;
            
            % Find time windows
            tmin = -0.5;
            tmax = 1;
            epochOn = cellfun(@(x) x(1), fTb.time) + tmin;
            epochOff = cellfun(@(x) x(end), fTb.time) + tmax;
            epochWin = [epochOn epochOff];
            
            % Reslice tables
            fTb = ce.SliceTimeSeries('feat', epochWin, 'Fill', 'bleed');
            rTb = ce.SliceTimeSeries('resp', epochWin, 'Fill', 'bleed');
            
            % Select units
%             recId = NP.SE.GetID(ce);
%             [~, uId] = LMV.Param.GetSelectedClusId(recId);
%             uInd = ismember(ce.clusTb.clusId, uId);
%             uInd = find(uInd);
            uInd = 1:height(ce.clusTb);
            
            % Construct out struct from ce
            for i = ce.numEpochs : -1 : 1
                s = struct;
                s.name = char(tv.stimId(i));
                s.tmin = tmin;
                s.tmax = tmax;
                s.spike_fs = 1/ops.rsBinSize;
                s.chs_original = [];
                s.chs_physical = [];
                s.duration = diff(fTb.time{i}([1 end]));
                s.soundf = 1/ops.rsBinSize;
                s.sound = fTb.mic_ni{i};
                s.aud = fTb.mic{i}';
                s.spikes_sua = cell2mat(rTb{i,uInd+1})';
                s.spike_wave = [];
                s.anin = fTb.speaker1_ni{i}';
                
                varNames = [{'cue1On', 'cue1Off', 'cue3On', 'cue3Off', 'stimOn', 'prodOn'} phones'];
                for n = 1 : numel(varNames)
                    vn = varNames{n};
                    s.(vn) = fTb.(vn){i};
                end
                
                if i == 1
                    sClus = table2struct(ce.clusTb(uInd,:));
                    for u = 1 : numel(sClus)
                        sClus(u).Amplitude = NaN;
                        sClus(u).ContamPct = NaN;
                        sClus(u).KSLabel = NaN;
                        sClus(u).amp = NaN;
                        sClus(u).ch = NaN;
                        sClus(u).id = sClus(u).clusId;
                        sClus(u).fr = NaN;
                        sClus(u).n_spikes = sClus(u).numSpikes;
                        sClus(u).sua_id = NaN;
                    end
                    s.sua_info = sClus;
                else
                    s.sua_info = [];
                end
                
                out(i) = s;
            end
        end
        
    end
    
end
