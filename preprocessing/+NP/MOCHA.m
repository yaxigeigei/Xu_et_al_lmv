classdef MOCHA
    
    methods(Static)
        function tv = AddTaskValueTablePTB(se) 
            % Make and add taskValue table to se
            
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
            
            NP.TaskBaseClass.MigrateStimId(se, '^sent_.{10,}$');
            [tt, tv] = se.GetTable('taskTime', 'taskValue');

            % Check required columns in taskTime table
            requiredCols = {'stim', 'green_cue', 'prod'};
            assert(all(ismember(requiredCols, tt.Properties.VariableNames)), ...
                "MOCHA task requires '%s' and '%s' columns in the taskTime table.", requiredCols{:});
            
            numTrials = 0;
            for i = 1 : height(tt)
                % Get task name
                eStim = tt.stim(i);
                if isnan(eStim) 
                    continue
                end
                eProd = tt.prod{i};
                
                taskName = lower(eStim(1).GetVfield('task'));
                if taskName ~= "mocha"
                    continue
                end
                numTrials = numTrials + 1;
                
                % Get stim and prod texts
                stimText = eStim(1).GetParentLabel;
                prodText = eProd.GetAllParentLabel();
                
                tv.taskName(i) = taskName;
                tv.stimText(i) = stimText;
                tv.prodText(i) = prodText;
            end

            if ~numTrials
                fprintf("This se does not have any MOCHA trial\n");
                return
            end
                
            % Make unique numeric stim IDs
            m = tv.taskName == "mocha";
            [~, ~, tv.stimNumId(m)] = unique(tv.stimText(m), 'stable');
            
            % Add taskValues table to se
            se.SetTable('taskValue', tv);
        end

        function tv = AddTaskValueTable(se)
            % Make and add taskValue table to se
            
            % Continue with existing table or create new
            if ismember('taskValue', se.tableNames)
                tv = se.GetTable('taskValue');
            else
                tv = table;
            end
            tt = se.GetTable('taskTime');
            tv.trialNum = (1:height(tt))';
            
            % Check required columns in taskTime table
            requiredCols = {'text', 'prod'};
            assert(all(ismember(requiredCols, tt.Properties.VariableNames)), ...
                "MOCHA task requires '%s' and '%s' columns in the taskTime table.", requiredCols{:});
            
            for i = 1 : height(tt)
                % Get task name
                eText = tt.text(i);
                eProd = tt.prod{i};
                if isnan(eText)
                    continue
                end
                taskName = lower(eText(1).GetVfield('task'));
                if taskName ~= "mocha"
                    continue
                end
                
                % Get stim and prod texts
                stimText = eText.V.type;
                prodText = eProd.GetAllParentLabel();
                
                tv.taskName(i) = taskName;
                tv.stimText(i) = stimText;
                tv.prodText(i) = prodText;
            end
            
            % Make unique numeric stim IDs
            m = tv.taskName == "mocha";
            [~, ~, tv.stimNumId(m)] = unique(tv.stimText(m), 'stable');
            
            % Add taskValue table to se
            se.SetTable('taskValue', tv, 'eventValues');
        end
        
        function [tt, tv] = MatchProd2Stim(tt, tv)
            % Match the transcript of speech producition to stimulus
            
            % Check required columns in taskTime table
            requiredCols = {'prod'};
            assert(all(ismember(requiredCols, tt.Properties.VariableNames)), ...
                "MatchProd2Stim requires the '%s' column in the taskTime table.", requiredCols{:});
            
            requiredCols = {'stimText'};
            assert(all(ismember(requiredCols, tv.Properties.VariableNames)), ...
                "MatchProd2Stim requires the '%s' column in the taskValue table.", requiredCols{:});
            
            % Helper function to preprocess strings
            sFunc = @(x) lower(erasePunctuation(x));
            
            for i = 1 : height(tt)               
                % Get stim label
                s1 = sFunc(tv.stimText(i));
                
                % Get production label
                prod = tt.prod{i};
                if isnan(prod)
                    tv.alignScore(i) = -Inf;
                    tv.tReact(i) = NaN;
                    continue
                end
                prod = Cut(prod); % cut into words
                s2 = sFunc(prod.GetParentLabel());
                
                % Align sequences
                [k1, k2, info] = MLing.FindAlignedWords(s1, s2);
                aa = info.aa;
                sc = info.tscore;
                
                % Normalize score to max score
                sc = sc / strlength(s1);
                
                % Find reaction time by the first prod event that matches stim
                tReact = prod(k2(1)).T.tmin;
                
                tv.alignment{i} = string(aa);
                tv.alignScore(i) = sc;
                tv.tReact(i) = tReact;
                
                tt.matchOn(i) = tReact;
                tt.matchOff(i) = prod(k2(end)).T.tmax;
            end
        end
        
        function [tt, tv] = FindMorphTimes(tt, tv)
            % Find matching key times between source and temaplte trials
            
            % Preallocate columns
            tt.morphFrom = num2cell(tt.trialOn);
            tt.morphTo = num2cell(tt.trialOn);
            tv.tempTrialNum = NaN(size(tt.trialOn));
            
            % Compute overall median reaction time
            rtMed = median(tv.tReact, 'omitnan');
            
            % Find unique stimuli
            [stimList, ~, ic] = unique(tv.stimText);
            
            for i = 1 : numel(stimList)
                % Find trials with the same stimulus
                ind = find(ic==i);
                
                % Find the trial(s) that has the highest alignment score
                sc = tv.alignScore(ind);
                isBest = sc==max(sc);
                indBest = ind(isBest);
                
                % Break ties by choosing the trial with production time closest to median
                dur = tt.matchOff(indBest) - tt.matchOn(indBest);
                [~, I] = min(abs(dur - median(dur)));
                indBest = indBest(I);
                
                % Get the template speech events
                tp = tt.prod{indBest};
                if all(isnan(tp))
                    warning("Production template for '%s' is empty.", stimList(i));
                    continue
                end
                tp = Cut(Cut(tp));
                
                % Align repeats and collect all the aligned times
                T = NaN(numel(ind), numel(tp));
                for j = 1 : numel(ind)
                    % Get the source speech events to be morphed
                    k = ind(j);
                    src = tt.prod{k};
                    if all(isnan(src))
                        continue % if there's no production event
                    end
                    src = Cut(Cut(src));
                    
                    % Find timestamps of matched words
                    [tTp, tSrc, info] = NP.MOCHA.FindMatchedSpeechTimes(tp, src);
                    
                    kTp = info.kk(1,:);
                    T(j,kTp) = tSrc;
                    
                    % Normalize alignment score to max score
                    info.nscore = info.tscore / numel(tp);
                    
                    % Add trial onset time (i.e. text onset in MOCHA, pic onset in SenGen, etc)
                    tTrialOn = tt.trialOn(k);
                    tSrc = [tTrialOn; tSrc];
                    
                    tv.alignInfo{k} = info;
                    tv.tempTrialNum(k) = tv.trialNum(indBest);
                    tt.morphFrom{k} = tSrc;
                end
                
                % Use median reaction time and median event intervals to construct template
                itvl = diff(T, 1, 2);
                medItvl = median(itvl, 1, 'omitnan');
                tTp = cumsum([rtMed medItvl])';
                
                % Find and save matching times from the template
                tpTxt = tt.trialOn(indBest);
                for j = 1 : numel(ind)
                    k = ind(j);
                    info = tv.alignInfo{k};
                    if isempty(info)
                        continue
                    end
                    kTp = info.kk(1,:);
                    tt.morphTo{k} = [tpTxt; tTp(kTp)];
                end
            end
        end
        
    end

end