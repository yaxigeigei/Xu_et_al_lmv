classdef CV
    
    properties(Constant)
        stimIdList = ["ba", "da", "ga", "dee", "dae", "do"];
    end
    
    methods(Static)
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
            
            % Set stim ID in taskValue
            tv = NP.CV.MigrateStimId(se, tv);
            
            % Check required columns in taskTime table
            requiredCols = {'stim', 'prod'};
            assert(all(ismember(requiredCols, tt.Properties.VariableNames)), ...
                "CV task requires the '%s' and '%s' columns in the taskTime table.", requiredCols{:});
            
            numTrials = 0;
            for i = 1 : height(tt)
                % Get task name
                eStim = tt.stim(i);
                eProd = tt.prod{i};
                if isnan(eStim)
                    continue
                end
                taskName = lower(eStim(1).GetVfield('task'));
                if taskName ~= "cv"
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
            m = tv.taskName == "cv";
            [~, ~, tv.stimNumId(m)] = unique(tv.stimText(m), 'stable');
            
            % Add taskValues table to se
            se.SetTable('taskValue', tv, 'eventValues');
        end
        
        function tv = AddTaskValueTableCustom(se)
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
            NP.TaskBaseClass.MigrateStimId(se, '^(ANT_|CAT_|SEMSR\d*_)');
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            
            % Check required columns in taskTime table
            requiredCols = {'prod'};
            assert(all(ismember(requiredCols, tt.Properties.VariableNames)), ...
                "CV task requires the '%s' columns in the taskTime table.", requiredCols{:});
            
            % Extract info from speech events
            numTrials = 0;
            for i = 1 : height(tt)
                % Get task name
                eProd = tt.prod{i};
                if isnan(eProd)
                    continue
                end
                taskName = lower(eProd(1).V.task);
                duration = eProd(1).T.tmax;
                if ~startsWith(taskName, "custom_cv")
                    continue
                end
                numTrials = numTrials + 1;
                
                % Get stim and prod texts
                dtype = class(eProd);
                if dtype == "NP.TGEvent"
                    prodText = eProd.GetAllParentLabel();
                elseif dtype == "MSessionExplorer.Event"
                    prodText = eProd(1).V.type;
                else
                    error("Unrecognized event data type '%s'", dtype);
                end
                
                tv.taskName(i) = taskName;
                tv.prodText(i) = prodText;
                tv.duration(i) = duration;
            end
            
            if ~numTrials
                fprintf("This se does not have any custom_cv trial\n");
                return
            end
            
            % Add taskValue table to se
            se.SetTable('taskValue', tv, 'eventValues');
        end
        
        function tv = MigrateStimId(se, tv)
            % Remove columns of NP.TGEvent objects from taskTime table whose names match the CV stimulus IDs 
            % and add the IDs to the taskValue table.
            % Note that because the ID events often preceed the trimmed stim events. This function therefore 
            % matches ID events with stim events by finding the pair closest in time.
            % 
            %   tv = MigrateStimId(se, tv)
            % 
            
            % Collect all valid events
            seLite = se.Duplicate({'taskTime'}, false);
            seLite.SliceSession(0, 'absolute');
            tt = seLite.GetTable('taskTime');
            varNames = string(tt.Properties.VariableNames);
            stimVars = [NP.CV.stimIdList, NP.CV.stimIdList+"On", NP.CV.stimIdList+"Off"];
            sid = cell(size(varNames));
            for i = numel(varNames) : -1 : 1
                vn = varNames(i);
                isRm(i) = ~isempty(find(vn==stimVars, 1));
                if isRm(i)
                    if iscell(tt.(vn))
                        v = tt.(vn){1};
                    else
                        v = tt.(vn);
                    end
                    if class(v) == "MSessionExplorer.Event"
                        sid{i} = v;
                    end
                end
            end
            sid = cat(1, sid{:});
            sid(isnan(sid)) = [];
            
            % Remove ID related columns
            se.SetColumn('taskTime', isRm, []);
            
            % Check the presence of valid IDs
            if isempty(sid)
                return
            end
            
            % Match sid events to the closest stim events
            stim = se.GetColumn('taskTime', 'stim');
            if iscell(stim)
                stim = cellfun(@(x) x(1), stim);
            end
            stim = stim + se.GetReferenceTime('taskTime');
            tStim = double(stim);
            tSID = double(sid);
            for i = 1 : numel(sid)
                [~, I] = min(abs(tStim - tSID(i)));
                tv.stimId(I) = sid(i).GetVfield('type');
            end
        end
        
        % not in use
        function tv = AddSentTaskValueTable(se)
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
            requiredCols = {'exp1', 'prod'};
            assert(all(ismember(requiredCols, tt.Properties.VariableNames)), ...
                "Sentence repetition task requires 'exp1' and 'prod' columns in the taskTime table.");
            
            for i = 1 : height(tt)
                % Get task name
                eExp1 = tt.exp1{i};
                eProd = tt.prod{i};
                if isnan(eExp1)
                    continue
                end
                taskName = lower(eExp1(1).GetVfield('task'));
                if taskName ~= "sentence_repeat"
                    continue
                end
                
                % Get stim and prod texts
                stimText = eExp1(1).GetParentLabel;
                prodText = eProd.GetAllParentLabel();
                
                tv.taskName(i) = taskName;
                tv.stimText(i) = stimText;
                tv.prodText(i) = prodText;
            end
            
            % Make unique numeric stim IDs
            m = tv.taskName == "sentence_repeat";
            [~, ~, tv.stimNumId(m)] = unique(tv.stimText(m), 'stable');
            
            % Add taskValues table to se
            se.SetTable('taskValue', tv, 'eventValues');
        end
        
    end
end

