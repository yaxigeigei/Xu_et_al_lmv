classdef SenGen < NP.TaskBaseClass
    
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
            
            NP.TaskBaseClass.MigrateStimId(se, '^pic\d+_.+_.+_.+$');
            [tt, tv] = se.GetTable('taskTime', 'taskValue');

            % Check required columns in taskTime table
            requiredCols = {'stim', 'prod'};
            assert(all(ismember(requiredCols, tt.Properties.VariableNames)), ...
                "Sentence Generation task requires '%s' and '%s' columns in the taskTime table.", requiredCols{:});
            
            numTrials = 0;
            for i = 1 : height(tt)
                % Get task name
                eStim = tt.stim(i);
                if isnan(eStim) 
                    continue
                end
                eProd = tt.prod{i};

                taskName = lower(eStim(1).GetVfield('task'));
                if taskName ~= "sentgen"
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
                fprintf("This se does not have any Sentence Gen trial\n");
                return
            end
                
            % Make unique numeric stim IDs
            m = tv.taskName == "sentgen";
            [~, ~, tv.stimNumId(m)] = unique(tv.stimText(m), 'stable');
            
            % Add taskValues table to se
            se.SetTable('taskValue', tv);
        end

        function tv = AddTaskValueTable(se, stimTb)
            % Make taskValue table
            
            % Continue with existing table or create new
            if ismember('taskValue', se.tableNames)
                tv = se.GetTable('taskValue');
            else
                tv = table;
            end
            tt = se.GetTable('taskTime');
            tv.trialNum = (1:height(tt))';
            
            % Check required columns in taskTime table
            requiredCols = {'cue', 'prod'};
            assert(all(ismember(requiredCols, tt.Properties.VariableNames)), ...
                "Sentence generation task requires the '%s' and '%s' columns in the taskTime table.", requiredCols{:});
            
            % Ground truth description of stimulus
            stimTb.description = string(stimTb.description);
            k = 0; % row counter
            
            for i = 1 : height(tt)
                % Get task name
                eCue = tt.cue(i);
                eProd = tt.prod{i};
                if isnan(eCue)
                    continue
                end
                taskName = lower(eCue(1).GetVfield('task'));
                if taskName ~= "sentgen"
                    continue
                end
                
                % Get stim and prod texts
                k = k + 1;
                stimText = stimTb.description(k);
                prodText = eProd.GetAllParentLabel();
                
                tv.taskName(i) = taskName;
                tv.stimText(i) = stimText;
                tv.prodText(i) = prodText;
            end
            
            % Make unique numeric stim IDs
            m = tv.taskName == "sentgen";
            [~, ~, tv.stimNumId(m)] = unique(tv.stimText(m), 'stable');
            
            % Add taskValues table to se
            se.SetTable('taskValue', tv, 'eventValues');
        end
        
        function [stimText, I] = SortStimByPOS(stimText, POS)
            % Sort stimulus by parts of speech
            %
            %   [stimText, I] = SortStimByPOS(stimText, POS)
            % 
            
            POS = upper(char(POS));
            assert(all(ismember(POS, 'SVO')), "POS must be a combination of 'S', 'V', and/or 'O'");
            
            parts = arrayfun(@(x) strsplit(x), stimText, 'Uni', false);
            parts = cat(1, parts{:});
            parts = parts(:, [2 4 6]); % take out S, V and O
            
            priority = arrayfun(@(x) find(x=='SVO',1), POS);
            [~, I] = sortrows(parts(:,priority));
            stimText = stimText(I);
        end

        % Morphing
        function [tt, tv] = FindVerbTimes(tt, tv)
            % Find the onset times of the first semantically relevant word in questions
            %
            %   [tt, tv] = FindVerbTime(tt, tv)
            % 

            % 
            % Verb can be start of "pushes" or start of "is" in "is pushed"
            tt.sentOn = NaN(height(tt), 1);
            tt.sentOff = NaN(height(tt), 1);

            tt.verbOn = NaN(height(tt), 1);
            tt.verbOff = NaN(height(tt), 1);

            % Action onsent will be start of pushes or start of "pushed" in "is pushed"
            tt.actionOn = NaN(height(tt), 1);
            tt.actionOff = NaN(height(tt), 1);
            
            for i = 1 : height(tt)
                eStim = tt.stim(i);
                taskName = lower(eStim(1).GetVfield('task'));
                if taskName ~= "sentgen"
                    continue
                end
                
                stimId = tv.stimId(i);
                parts = split(stimId, '_');
                verb = parts{2};
                type = parts{3};
                eProd = tt.prod{i};
                tier = eProd.tier{1};
                
                % If type is passive then the keyword is "is" else it is verb
                if type == "passive"
                    search = "is";
                else
                    search = verb;
                end
                tt.sentOn(i) = eProd.T.tmin;
                tt.sentOff(i) = eProd.T.tmax;

                searchIdx = find(contains(tier.Label, search), 1, 'first');
                tt.verbOn(i) = eProd.T.T1(searchIdx);
                tt.verbOff(i) = eProd.T.T2(searchIdx);

                verbIdx = find(contains(tier.Label, verb), 1, 'first');
                tt.actionOn(i) = eProd.T.T1(verbIdx);
                tt.actionOff(i) = eProd.T.T1(verbIdx);
            end
        end
        
    end
end

