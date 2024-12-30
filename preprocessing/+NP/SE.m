classdef SE
    
    methods(Static)
        % IO
        function [seArray, filePaths] = LoadSession(varargin)
            % Load session files as MSessionExplorer
            % 
            %   [seArray, filePaths] = NP.SE.LoadSession()
            %   [seArray, filePaths] = NP.SE.LoadSession(filePaths)
            %   [seArray, filePaths] = NP.SE.LoadSession(..., 'Enrich', false)
            %   [seArray, filePaths] = NP.SE.LoadSession(..., 'UserFunc', @(se) se)
            % 
            
            % Handle user inputs
            p = inputParser();
            p.addOptional('filePaths', '', @(x) ischar(x) || iscellstr(x) || isstring(x) || isempty(x));
            p.addParameter('Enrich', false, @islogical);
            p.addParameter('UserFunc', @(x) x, @(x) isa(x, 'function_handle'));
            p.parse(varargin{:});
            filePaths = p.Results.filePaths;
            isEnrich = p.Results.Enrich;
            userFunc = p.Results.UserFunc;
            
            if isempty(filePaths)
                filePaths = MBrowse.Files();
            else
                filePaths = cellstr(filePaths);
            end
            if isempty(filePaths)
                seArray = [];
                return;
            end
            
            % Preallocation
            seArray(numel(filePaths),1) = MSessionExplorer();
            
            for i = 1 : numel(filePaths)
                load(filePaths{i}, 'se');
                fprintf("\nLoading se %i - %s\n", i, NP.SE.GetID(se));
                if isEnrich
                    NP.SE.Enrich(se);
                end
                userFunc(se);
                seArray(i) = se;
            end
        end
        
        % Getters
        function [recId, subjectId, blockId] = GetID(arg)
            % Extract recording ID and subject ID, e.g. "NP00_B0", "NP00"
            % 
            %   [recId, subjectId, blockId] = NP.SE.GetID(se)
            %   [recId, subjectId, blockId] = NP.SE.GetID(se.userData)
            %   [recId, subjectId, blockId] = NP.SE.GetID(se.userData.expInfo)
            %   [recId, subjectId, blockId] = NP.SE.GetID(fileName)
            %   [recId, subjectId, blockId] = NP.SE.GetID(clusId)
            % 
            
            % Get userdata struct(s) if input is a scalar MSessionExplorer
            if isa(arg, 'MSessionExplorer') && isscalar(arg)
                arg = arg.userData;
            end
            
            % Recursively extract IDs for array input
            if ~ischar(arg) && numel(arg) > 1
%                 fprintf("Found metadata of more than one recordings. Combine IDs.\n");
                [recId, subjectId, blockId] = arrayfun(@NP.SE.GetID, arg, 'Uni', false);
                recId = strjoin(unique(recId), ',');
                subjectId = strjoin(unique(subjectId), ',');
                blockId = strjoin(unique(blockId), ',');
                return
            end
            
            % Denesting
            if iscell(arg)
                arg = arg{1};
            end
            
            if isstruct(arg) && isfield(arg, 'expInfo')
                % When input is a se.userData struct
                s = arg.expInfo;
            elseif isstruct(arg) && isfield(arg, 'subjectId')
                % When input is a se.userData.expInfo struct
                s = arg;
            elseif ischar(arg) || isstring(arg)
                % When input is a string that contains recording IDs
                nameParts = strsplit(char(arg), '_');
                s.subjectId = nameParts{1};
                s.blockId = nameParts{2};
            elseif isnumeric(arg)
                % When input is a cluster ID
                idStr = num2str(arg);
                s.subjectId = "NP"+idStr(1:end-7);
                s.blockId = "B"+str2double(idStr(end-6:end-5));
            else
                error('Cannot resolve ID info from the input variable.');
            end
            
            subjectId = char(s.subjectId);
            blockId = char(s.blockId);
            recId = [subjectId '_' blockId];
        end
        
        function [name, region, subregion] = GetRegion(arg)
            % 
            
            if isa(arg, 'MSessionExplorer')
                % When input is a MSessionExplorer object
                s = arg.userData.recMeta;
            elseif isstruct(arg) && isfield(arg, 'recMeta')
                % When input is a se.userData struct
                s = arg.recMeta;
            elseif isstruct(arg) && isfield(arg, 'Region')
                % When input is a se.userData.recMeta struct
                s = arg;
            else
                error('The input data is incorrect');
            end
            
            region = s.Region;
            subregion = s.Subregion;
            
            if contains(region, 'STG')
                name = 'STG';
            elseif contains(region, 'MTG')
                name = 'MTG';
            else
                name = region;
            end
        end
        
        function tRange = GetRecTimeRange(se)
            % Find the absolute time of the first and last sample
            
            tRange = [NaN NaN];
            
            % Iterate through all tables
            for i = 1 : numel(se.tableNames)
                if se.isEventValuesTable(i)
                    continue
                end
                tb = se.GetTable(se.tableNames{i});
                tRef = se.GetReferenceTime(se.tableNames{i});
                
                if se.isTimesSeriesTable(i)
                    colInd = 1; % only check timestamps
                elseif se.isEventTimesTable(i)
                    colInd = 1 : width(tb); % check all event times
                end
                
                % Iterate through columns of interest
                for c = colInd
                    if iscell(tb.(c))
                        val = tb.(c){1};
                        if ~isempty(val)
                            t1 = double(tb.(c){1}(1)) + tRef(1);
                        end
                        val = tb.(c){end};
                        if ~isempty(val)
                            tn = double(tb.(c){end}(end)) + tRef(end);
                        end
                    else
                        t1 = double(tb.(c)(1)) + tRef(1);
                        tn = double(tb.(c)(end)) + tRef(end);
                    end
                    tRange(1) = min(tRange(1), t1);
                    tRange(2) = max(tRange(2), tn);
                end
            end
        end
        
        % Enrich
        function se = Enrich(se, ops)
            % Enrich featuers in se
            
            % Speech
            if ops.isFiltSpeech
                NP.Audio.FilterSpeech(se);
            end
            if ops.isMel
                NP.Audio.AddMelTable(se);
            end
            if ops.isPitch
                NP.Pitch.EnrichPitchTable(se);
            end
            if ops.isArtic
                NP.Artic.EnrichArticTable(se);
            end
            if ops.isTimitGT
                NP.Phone.AddTimitGroundTruth(se);
            end
            
            % Neural
            if ops.isMergeMeta
                NP.Unit.MergeMeta(se);
            end
            if ops.isSpkRate
                NP.Unit.AddSpikeRateTable(se, ops);
            end
            if ops.isSpkSpan
                NP.Unit.AddSpikeSpanTable(se);
            end
        end
        
        function AddSubjectMeta(seArray, metaTb)
            % Find subject metadata in the input table and add them to the userData of the se objects
            % 
            %   NP.SE.AddSubjectMeta(seArray)
            %   NP.SE.AddSubjectMeta(seArray, metaTb)
            % 
            
            if nargin < 2
                metaTb = NP.Preproc.ImportMetaSheet([], "Subjects");
            end
            
            for i = 1 : numel(seArray)
                [~, subjectId] = NP.SE.GetID(seArray(i));
                m = strcmp(subjectId, metaTb.(1));
                if any(m)
                    seArray(i).userData.subjectMeta = table2struct(metaTb(m,:));
                else
                    error("Cannot find the metadata for '%s'", subjectId);
                end
            end
        end
        
        function AddRecMeta(seArray, metaTb)
            % Find recording metadata in the input table and add them to the userData of the se objects
            % 
            %   NP.SE.AddRecMeta(seArray)
            %   NP.SE.AddRecMeta(seArray, metaTb)
            % 
            
            if nargin < 2
                metaTb = NP.Preproc.ImportMetaSheet([], "Recordings");
            end
            
            recIdList = metaTb.Subject + "_" + metaTb.Block;
            
            for i = 1 : numel(seArray)
                recId = NP.SE.GetID(seArray(i));
                m = recId == recIdList;
                if any(m)
                    seArray(i).userData.recMeta = table2struct(metaTb(m,:));
                else
                    error("Cannot find the metadata for '%s'", recId);
                end
            end
        end
        
        function tv = AddTaskValueTable(se, varargin)
            % Add taskValue table to se
            
            p = inputParser();
            p.KeepUnmatched = true;
            p.addParameter('taskNames', {}, @iscellstr);
            p.parse(varargin{:});
            taskNames = lower(cellstr(p.Results.taskNames));
            
            fprintf("\nAdding taskValue table to se\n");
            
            for i = 1 : numel(taskNames)
                switch taskNames{i}
                    case 'lmv'
                        tv = LMV.SE.AddTaskValueTable(se);
                    case 'sentgen'
                        if NP.SE.GetID(se) == "NP38_B5"
                            stimTb = readtable("sentence_generation_stim_20220830.xlsx");
                            stimTb = stimTb([5 6 8:17 3:45], :);
                            tv = NP.SenGen.AddTaskValueTable(se, stimTb);
                        else
                            tv = NP.SenGen.AddTaskValueTablePTB(se);
                        end
                    case 'timit'
                        tv = NP.TIMIT.AddTaskValueTable(se);
                    case 'mocha'
                        tv = NP.MOCHA.AddTaskValueTablePTB(se);
                    case {'semsr', 'semsr1'}
                        tv = Semsr.SE.AddTaskValueTable(se);
                    case 'cv'
                        tv = NP.CV.AddTaskValueTable(se);
                    case 'custom_cv'
                        tv = NP.CV.AddTaskValueTableCustom(se);
                    case 'arithmetic'
                        tv = NP.Arithmetic.AddTaskValueTable(se);
                    case 'timit_read'
                        tv = NP.PTBRead.AddTaskValueTable(se);
                    case 'ptb_read'
                        tv = NP.PTBRead.AddTaskValueTable(se);
                    case 'nbd_listen'
                        tv = NP.NBDListen.AddTaskValueTable(se);
                    case 'nbd'
                        tv = NP.NBD.AddTaskValueTable(se);
                    otherwise
                        warning("Cannot add info to taskValue table. '%s' is not a defined task name.", taskNames{i});
                end
            end
            
            se.Peek('taskValue');
        end
        
        function FixNonMonotonicTimestamps(se)
            % Detect and remove parts of timeseries where timestamps are not monotonically increasing
            
            for i = 1 : numel(se.tableNames)
                if ~se.isTimesSeriesTable(i)
                    continue
                end
                
                tb = se.GetTable(se.tableNames{i});
                
                for j = 1 : height(tb)
                    % Find the starting indices of non-monotonic parts
                    t = tb.time{j};
                    indNMT = find(diff(t) <= 0);
                    if isempty(indNMT)
                        continue
                    else
                        fprintf("Found %i non-monotonicity in '%s' at epoch %i.\n", numel(indNMT), se.tableNames{i}, j);
                    end
                    
                    % Set the mask of non-monotomic parts to false
                    m = true(size(t));
                    for a = indNMT(:)'
                        b = a + find(t(a+1:end) > t(a), 1);
                        m(a:b) = false;
                        fprintf("  Remove samples from %i to %i.\n", a, b);
                    end
                    
                    % Select monotonically increase parts
                    for k = 1 : width(tb)
                        tb.(k){j} = tb.(k){j}(m,:);
                    end
                end
                
                se.SetTable(se.tableNames{i}, tb);
            end
        end
        
        function UpdateConventions(se)
            % Update variable naming and other conventions in old se
            
            % Check and refactor variable names
            if isfield(se.userData, 'experimentInfo')
                warning("Refactor variable names for expInfo");
                se.userData = renameStructField(se.userData, 'experimentInfo', 'expInfo');
                se.userData.expInfo = renameStructField(se.userData.expInfo, 'experimentId', 'recId');
            end
        end
        
        % Slicing
        function TrimRecording(se, varargin)
            % Remove data at the beginning and the end
            
            p = inputParser();
            p.KeepUnmatched = true;
            p.addParameter('trimWin', 'sorting', @(x) isnumeric(x) || ismember(x, {'sorting'}));
            p.parse(varargin{:});
            win = p.Results.trimWin;
            
            if ischar(win) || isstring(win)
                switch char(lower(win))
                    case 'sorting'
                        probe = 1;
                        ksOps = se.userData.ksMeta(probe).ops;
                        tWin = ksOps.trange;
                        tWin = tWin * ksOps.fs / se.userData.apMeta(probe).imSampRate;
                    otherwise
                        error("'%s' is not a valid option for trimWin", ops.trimWin);
                end
            else
                tWin = win;
            end
            
            fprintf("\nTrim recording before %.3fs and after %.3fs\n", tWin(1), tWin(2));
            
            se.SliceSession(tWin, 'absolute');
            se.RemoveEpochs(2);
        end
        
        function SliceToTrials(se, varargin)
            % Slice recording into epochs
            
            assert(se.numEpochs==1, "The input se should only has one epoch, instead it has %i", se.numEpochs);
            
            p = inputParser();
            p.KeepUnmatched = true;
            p.addParameter('taskNames', {}, @iscellstr);
            p.parse(varargin{:});
            taskNames = lower(cellstr(p.Results.taskNames));
            
            % Different tasks use different delimiter events
            delimLookup = struct;
            delimLookup.timit = 'stim';
            delimLookup.lmv = 'cue1';
            delimLookup.sentgen = 'stim';
            delimLookup.semsr = 'stim';
            delimLookup.semsr1 = 'stim';
            delimLookup.cv = 'stim';
            delimLookup.custom_cv = 'prod';
            delimLookup.arithmetic = 'stim';
            delimLookup.timit_read = 'stim';
            delimLookup.ptb_read = 'stim';
            delimLookup.nbd_listen = 'stim';
            delimLookup.nbd = 'stim';
            delimLookup.mocha = 'stim';
            
            % Find delimiter events
            tt = se.GetTable('taskTime');
            delimEvents = cell(size(taskNames));
            for i = 1 : numel(taskNames)
                if ~isfield(delimLookup, taskNames{i})
                    warning("Not slicing the '%s' task since the delimiter event is not specified.", taskNames{i});
                    continue
                end
                delimName = delimLookup.(taskNames{i});
                delimEvts = cat(1, tt.(delimName){:});
                m = strcmp(delimEvts.GetVfield('task'), taskNames{i});
                delimEvents{i} = double(delimEvts(m));
            end
            
            % Set delimiter event times as trialOn
            trialOn = cat(1, delimEvents{:});
            trialOn = unique(trialOn); % this will also sort the timestamps
            se.SetColumn('taskTime', 'trialOn', {trialOn});
            
            % Make slicing time
            rt = se.GetReferenceTime;
            tSlice = trialOn + rt;
            tSlice(1) = rt;
            
            % Slice recording
            fprintf("\nSlice recording into trials\n");
            se.SliceSession(tSlice, 'absolute');
            
            % Realign time to trialOn
            se.AlignTime('trialOn', 'taskTime');
            
            % Make sure certain variables are in standard data types
            tt = se.GetTable('taskTime');
            if ismember('prod', tt.Properties.VariableNames)
                if ~iscell(tt.prod)
                    tt.prod = num2cell(tt.prod);
                end
            end
            se.SetTable('taskTime', tt);
        end
        
        function se = BleedTrials(se, tOffset)
            % Slice each epoch or trial to include a period before and after (i.e. bleeding into adjacant trials)
            % 
            %   se = BleedTrials(se)
            %   se = BleedTrials(se, tOffset)
            % 
            % Inputs
            %   se              The input se will not be modified.
            %   tOffset         Default is [-1 1] to include one sec before one sec after.
            % Output
            %   se              Sliced se.
            % 
            
            if nargin < 2
                tOffset = [-1 1];
            end
            
            se = se.Duplicate;
            
            % Align time at trialOn
            se.AlignTime('trialOn', 'taskTime');
            
            % Compute time windows for slicing
            trStart = zeros(se.numEpochs, 1);
            dur = diff(se.GetReferenceTime);
            trEnd = trStart + [dur; 10]; % hardcode using 10s for the last trial
            tWin = [trStart trEnd] + tOffset(:)';
            
            % Slice se
            fprintf("Slicing each trial to include a period before and after.\n")
            for i = 1 : numel(se.tableNames)
                if se.isEventValuesTable(i)
                    continue;
                end
                if se.tableNames{i} == "taskTime"
                    tb = se.GetTable(se.tableNames{i});
                elseif se.isEventTimesTable(i)
                    tb = se.SliceEventTimes(se.tableNames{i}, tWin, 'Fill', 'bleed');
                else
                    tb = se.SliceTimeSeries(se.tableNames{i}, tWin, 'Fill', 'bleed');
                end
                se.SetTable(se.tableNames{i}, tb);
            end
        end
        
        % Time morphing
        function SetMorphTimes(se, morphFor, varargin)
            % Add morphing info to the 'taskTime' and/or 'taskValue' table
            % 
            %   NP.SE.SetMorphTimes(se)
            %   NP.SE.SetMorphTimes(se, morphFor)
            %   NP.SE.SetMorphTimes(se, 'lmv')
            %   NP.SE.SetMorphTimes(se, 'lmv', targetType)
            %   NP.SE.SetMorphTimes(se, 'sentgen')
            %   NP.SE.SetMorphTimes(se, 'semsr')
            %   NP.SE.SetMorphTimes(se, 'events', eventNames)
            % 
            % Inputs
            %   se              The input MSessionExplorer object.
            %   morphFor        'lmv', 'sentgen', 'semsr', or 'events'. Will use taskName in the taskValue 
            %                   table if this argument is not provided or is empty [].
            %   targetType      When morphFor is 'lmv', targetType specifies what times to align to.
            %                   Default is 'median', or one can use 'best'. See the LMV.SE.FindMorphTimes 
            %                   for details.
            %   eventNames      When morphFor is 'events', eventNames should be an array of strings that 
            %                   specifies what events (in the taskTime table) to align to after morphing.
            
            % Use trial name in morphFor
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            rt = se.GetReferenceTime('taskTime');
            if nargin < 2 || isempty(morphFor)
                fprintf("'morphFor' is not provided. Will use task name as the value.\n");
                morphFor = unique(tv.taskName);
                assert(isscalar(morphFor), "The se contains multiple tasks. Consider spliting the task conditions first.");
            end
            
            % Match trials and find morph times
            if morphFor == "lmv"
                % LMV task - morph times are determined using data medians
                fprintf("\nMatch %s trials and find morph times\n", morphFor);
                [tt, tv] = LMV.SE.MatchProd2Stim(tt, tv, rt);
                [tt, tv] = LMV.SE.FindMorphTimes(tt, tv, rt, varargin{:});
                
            elseif morphFor == "sentgen"
                % SentGen task
                fprintf("\nMatch %s trials and find morph times\n", morphFor);
                [tt, tv] = NP.SenGen.FindVerbTimes(tt, tv);
                % You can use actionOn as the reference time instead of verbOn to align "pushed" to "pushes" instead of "is" to "pushes"
                tt = NP.TaskBaseClass.FindMorphTimes(tt, {'stimOn', 'sentOn', 'verbOn', 'stimOff'});
            
            % elseif morphFor == "mocha"
            %     fprintf("\nFind morph times to align %s events\n", morphFor);
            %     tt = NP.TaskBaseClass.FindMorphTimes(tt, {'stimOn', 'prodOn'});
                
            elseif morphFor == "semsr"
                % Semsr task
                fprintf("\nFind morph times to align %s events\n", morphFor);
                tt = NP.TaskBaseClass.FindMorphTimes(tt, {'keyOn', 'stimOff', 'tReact'});

            elseif morphFor == "arithmetic"
                % Arithmetic task
                fprintf("\nFind morph times to align %s events\n", morphFor);
                [tt, tv] = NP.Arithmetic.FindKeyTimes(tt, tv);
                tt = NP.TaskBaseClass.FindMorphTimes(tt, {'arithmetic_firstOperandOn', 'arithmetic_firstOperandOff', 'arithmetic_operationOn','arithmetic_operationOff', 'arithmetic_secondOperandOn', 'arithmetic_secondOperandOff', 'arithmetic_equalsOn','arithmetic_equalsOff','arithmetic_responseOn'});
                
            elseif morphFor == "ptb_read"
                % PTBRead task
                fprintf("\nFind morph times to align %s events\n", morphFor);
                [tt, tv] = NP.PTBRead.FindKeyTimes(tt, tv);
                tt = NP.TaskBaseClass.FindMorphTimes(tt, {'firstWordOn', 'lastWordOff'});
            
            elseif morphFor == "nbd_listen"
                % PTBRead task
                fprintf("\nFind morph times to align %s events\n", morphFor);
                [tt, tv] = NP.NBDListen.FindKeyTimes(tt, tv);
                tt = NP.TaskBaseClass.FindMorphTimes(tt, {'firstWordOn', 'lastWordOff'});

            elseif morphFor == "nbd"
                % NBD task
                fprintf("\nFind morph times to align %s events\n", morphFor);
                [tt, tv] = NP.NBD.FindKeyTimes(tt, tv);
                tt = NP.TaskBaseClass.FindMorphTimes(tt, {'firstWordOn', 'lastWordOn', 'lastWordOff'});

            elseif morphFor == "timit"
                % TIMIT task
                fprintf("\nFind morph times to align %s events\n", morphFor);
                [tt, tv] = NP.TIMIT.FindKeyTimes(tt, tv);
                tt = NP.TaskBaseClass.FindMorphTimes(tt, {'stimOn', 'stimOff'});

            elseif morphFor == "events"
                % Morph to align specific events
                fprintf("\Find morph times to align the specified events\n");
                tt = NP.TaskBaseClass.FindMorphTimes(tt, varargin{:});
                
            else
                fprintf("\nNo morphing for '%s' task\n", morphFor);
                return
            end
            se.SetTable('taskTime', tt);
            se.SetTable('taskValue', tv);
        end
        
        function se = MorphSession(se, varargin)
            % Apply time mapping to all data in a SE object. This method requires reference times
            % to be available, and the global time continuity will be maintained.
            % 
            %   seMorphed = NP.SE.MorphSession(se)
            % 
            % Input
            %   se              The original MSessionExplorer object with the following data:
            %                   1) 'morphFrom' and 'morphTo' columns in the 'taskTime' table. See 
            %                      the NP.SE.SetMorphTimes function to learn more.
            %                   2) Reference times.
            % Output
            %   seMorphed       The time-morphed MSessionExplorer object. 
            %                   A timeSeries table named 'morph' will be added or updated to track 
            %                   the change of time and speed. The 'time0' column stores the most 
            %                   original timestamps no matter how many consecutive morphings have 
            %                   been applied to the object. The 'speed' column stores the relative 
            %                   speed wrt the most original speed. For example, a speed of 2 means 
            %                   twice the original speed, and 0.5 means half of that.
            % 
            % See also NP.SE.SetMorphTimes, NP.Unit.AddSimSpikeTimeTable, NP.SE.MorphEpochs
            
            p = inputParser();
            p.KeepUnmatched = true;
            p.addParameter('isMorph', true, @islogical);
            p.parse(varargin{:});
            isMorph = p.Results.isMorph;
            
            if isMorph
                fprintf("\nApply time morphing...\n");
            else
                disp("Time morphing is not performed since 'isMorph' is set to false.");
                return
            end
            
            se = se.Duplicate();
            isVerbose = se.isVerbose;
            se.isVerbose = false;
            
            % Align time to trial onset
            % i.e. text onset in MOCHA, pic onset in SenGen, cue1 onset in LMV
            fprintf("Align trials at 'trialOn'\n");
            se.AlignTime('trialOn', 'taskTime');
            
            % Find and add up changes in trial duration
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            dt = cellfun(@(x,y) x(end)-y(end), tt.morphTo, tt.morphFrom);
            dt = cumsum([0; dt(1:end-1)]);
            tt.morphTo = cellfun(@(x,y) x+y, tt.morphTo, num2cell(dt), 'Uni', false);
            se.SetTable('taskTime', tt);
            
            % Slice se
            fprintf("Vectorize session to a single epoch\n");
            se.SliceSession(0, 'absolute');
            
            % Add a morph table to keep track of speed change
            if ismember('morph', se.tableNames)
                fprintf("Reuse the 'morph' table to keep track of speed change\n");
            else
                fprintf("Add a 'morph' table to keep track of speed change\n");
                tRange = NP.SE.GetRecTimeRange(se);
                itvl = 0.01; % sampling interval in sec
                mTb = table;
                mTb.time = {(tRange(1) : itvl : tRange(2))'};
                mTb.time0 = mTb.time;
                se.SetTable('morph', mTb, 'timeSeries', se.GetReferenceTime);
            end
            
            % Construct Morpher objects
            fprintf("Morph session\n");
            tt = se.GetTable('taskTime');
            mObj = NP.Morpher(tt.morphFrom, tt.morphTo);
            mObj.MorphSE(se);
            
            % Compute resulting speed
            fprintf("Compute speed change\n");
            mtb = se.GetTable('morph');
            mtb.speed = cellfun(@(t,t0) gradient(t)./gradient(t0), mtb.time, mtb.time0, 'Uni', false);
            se.SetTable('morph', mtb);
            
            % Slice se back
            fprintf("Slice session back to trial epochs\n");
            tt = se.GetTable('taskTime');
            se.SliceSession(tt.trialOn{1}, 'absolute');
            se.SetTable('taskValue', tv, 'eventValues');
            se.isVerbose = isVerbose;
            
            % Replace morphed spike times with simulated spike times
            if ismember('spikeRate', se.tableNames)
                NP.Unit.AddSimSpikeTimeTable(se, 'spikeTime');
            end
        end
        
        function se = MorphEpochs(se, varargin)
            % Apply time mapping epoch by epoch in a SE object. This method does not require reference 
            % times to be available, but as a reuslt the global time continuity will be broken.
            % 
            %   seMorphed = NP.SE.MorphEpochs(se)
            % 
            % Input
            %   se              The original MSessionExplorer object with the following data:
            %                   1) 'morphFrom' and 'morphTo' columns in the 'taskTime' table. See 
            %                      the NP.SE.SetMorphTimes function to learn more.
            % Output
            %   seMorphed       The time-morphed MSessionExplorer object. 
            %                   A timeSeries table named 'morph' will be added or updated to track 
            %                   the change of time and speed. The 'time0' column stores the most 
            %                   original timestamps no matter how many consecutive morphings have 
            %                   been applied to the object. The 'speed' column stores the relative 
            %                   speed wrt the most original speed. For example, a speed of 2 means 
            %                   twice the original speed, and 0.5 means half of that.
            % 
            % See also NP.SE.SetMorphTimes, NP.Unit.AddSimSpikeTimeTable, NP.SE.MorphSession
            
            p = inputParser();
            p.KeepUnmatched = true;
            p.addParameter('isMorph', true, @islogical);
            p.parse(varargin{:});
            isMorph = p.Results.isMorph;
            
            if isMorph
                fprintf("\nApply time morphing...\n");
            else
                disp("Time morphing is not performed since 'isMorph' is set to false.");
                return
            end
            
            se = se.Duplicate();
            
            % Construct Morpher objects
            fprintf("Morph epochs\n");
            tt = se.GetTable('taskTime');
            mObj = NP.Morpher(tt.morphFrom, tt.morphTo);
            
            % Morph time
            mObj.MorphSE(se);
            
            % Compute resulting speed
            if ismember('morph', se.tableNames)
                fprintf("Compute speed change\n");
                mtb = se.GetTable('morph');
                mtb.speed = cellfun(@(t,t0) gradient(t)./gradient(t0), mtb.time, mtb.time0, 'Uni', false);
                se.SetTable('morph', mtb);
            else
                fprintf("Time and speed change is not tracked because the 'morph' table is not found. " + ...
                    "It cannot be created when morphng by epochs.\n");
            end
            
            % Replace morphed spike times with simulated spike times
            if ismember('spikeRate', se.tableNames)
                NP.Unit.AddSimSpikeTimeTable(se, 'spikeTime');
            end
        end
        
        % Split and group
        function seTb = SplitConditions(se, varargin)
            % Split an SE into a table of SEs by conditions in the behavValue table
            % Trials with NaN condition will be excluded
            % 
            %   seTb = SplitConditions(se)
            %   seTb = SplitConditions(..., 'ConditionVars', {})
            %   seTb = SplitConditions(..., 'ConditionVars', {}, 'SourceTable', 'taskValue')
            %   seTb = SplitConditions(se, ops)
            % 
            
            p = inputParser();
            p.KeepUnmatched = true;
            p.addParameter('ConditionVars', {}, @(x) iscellstr(x) || isstring(x) || ischar(x));
            p.addParameter('SourceTable', 'taskValue', @(x) ismember(x, se.tableNames));
            p.parse(varargin{:});
            condVars = cellstr(p.Results.ConditionVars);
            srcTb = p.Results.SourceTable;
            
            % Find groups
            if isempty(condVars)
                % Initialize table with a dummy grouping variable
                dummyCond = ones(se.numEpochs, 1);
                T = table(dummyCond);
            else
                % Get and modify variables from behavValue table
                tv = se.GetTable(srcTb);
                T = tv(:, condVars);
            end
            [condId, seTb] = findgroups(T);
            
            % Split SE by conditions
            for i = 1 : max(condId)
                % Remove non-member
                seCopy = se.Duplicate();
                seCopy.RemoveEpochs(condId ~= i);
                
                % Add to table
                [seTb.recId{i}, seTb.subjectId{i}] = NP.SE.GetID(seCopy);
                seTb.se(i) = seCopy;
                seTb.numTrial(i) = seCopy.numEpochs;
            end
        end
        
        function condTb = CombineConditions(condTb, seTb, varargin)
            % Combine the same conditions across rows of an (usually concatenated) seTb
            
            p = inputParser;
            p.addParameter('UniformOutput', false);
            p.parse(varargin{:});
            isUni = p.Results.UniformOutput;
            
            % Take out condition columns in seTb
            isCond = ismember(seTb.Properties.VariableNames, condTb.Properties.VariableNames);
            seCond = seTb(:,isCond);
            seTb = seTb(:,~isCond);
            
            % Process each row of condTb
            nCond = width(condTb);
            for i = 1 : height(condTb)
                % Find the same condition across sessions
                isCond = true(height(seTb),1);
                for j = 1 : nCond
                    isCond = isCond & ismember(seCond.(j), condTb.(j)(i));
                end
                if ~any(isCond)
                    continue;
                end
                
                % Concatenate data
                for j = 1 : width(seTb)
                    vn = seTb.Properties.VariableNames{j};
                    col = seTb.(vn);
                    if iscell(col)
                        s2 = cellfun(@(x) size(x,2), col);
                        isCatable = numel(unique(s2)) == 1;
                        if ~iscellstr(col) && isCatable
                            condTb.(vn){i} = cat(1, col{isCond});
                            continue
                        end
                    end
                    condTb.(vn){i} = col(isCond);
                end
            end
            
            % Remove empty conditions
            condTb(cellfun(@isempty, condTb.(nCond+1)), :) = [];
            
            % Denest
            if ~isUni
                return;
            end
            for j = nCond+1 : width(condTb)
                try
                    condTb.(j) = cat(1, condTb.(j){:});
                catch
                    warning('Column %s cannot be denested', condTb.Properties.VariableNames{j});
                end
            end
        end
        
        function seCat = CatEpochs(seArray)
            % Concatenate multiple se objects into one with incremented time and added recording info
            % 
            %   seCat = NP.SE.CatEpochs(seArray)
            % 
            
            % Make a hard copy of seArray
            seArray = arrayfun(@Duplicate, seArray);
            
            rt0 = 0;
            maxTrial = 0;
            for i = 1 : numel(seArray)
                % Make sure cluster IDs are unique
                se = seArray(i);
                NP.Unit.SetUniqueClusId(se);
                
                % Offset reference times such that one recording follows another with a small gap
                tRange = NP.SE.GetRecTimeRange(se);
                rt = se.GetReferenceTime - tRange(1) + rt0;
                se.SetReferenceTime(rt);
                rt0 = rt0 + diff(tRange) + 1; % +1 for an one sec gap
                
                % Add recording IDs to trials in taskValue table
                tv = se.GetTable('taskValue');
                tv.recId(:) = string(NP.SE.GetID(se));
                
                % Renumber trials such that one recording follows another
                tv.trialNum = tv.trialNum + maxTrial;
                maxTrial = max(tv.trialNum);
                se.SetTable('taskValue', tv);
            end
            
            % Merging
            seCat = Merge(seArray);
        end
        
        function seArray = PartitionTrials(se, k)
            % Randomly split trials into partitions
            %   k is the number of fold as in KFold crossvalidation
            %   se can be an array and the partition populates the second dimension
            
            if numel(se) == 1
                % Partition
                c = cvpartition(se.numEpochs, 'Kfold', k);
                for i = k : -1 : 1
                    seArray(i) = se.Duplicate();
                    seArray(i).RemoveEpochs(~c.test(i));
                end
            else
                % Recursively process each se
                for i = numel(se) : -1 : 1
                    seArray(i,:) = SL.SE.PartitionTrials(se(i), k);
                end
            end
        end
        
        % Resampling
        function featTb = ResampleFeatures(se, ops)
            % Resample feature variables
            % 
            %   featTb = ResampleFeatures(se, ops)
            % 
            % Inputs
            %   se                  A MSessionExplorer object.
            %   ops                 A struct with the following options:
            %   ops.featVars        A struct where fieldnames denote the se tables of interest. Each field values 
            %                       is a cell array of variable (i.e. column) names of that table.
            %   ops.rsWin           A n-by-2 array of time window(s). When n == 1, the same window is used to 
            %                       resample all epochs in se; when n == se.numEpochs, each window is used for 
            %                       the corresponding epoch.
            %   ops.rsBinSize       Bin size or sample size for resampling.
            %   ops.rsShifts        A scalar or vector of time shift(s). Each value will create a time-shifted 
            %                       version of the feature timeseries. All the shifted timeseries of each feature 
            %                       will be concatenated horiontally first before concatenation across features.
            %                       Default is 0.
            %   ops.rsArgs          Interpolation parameters for resampling.
            % 
            % Output
            %   featTb              A timeSeries table of extracted and resampled features.
            % 
            
            % Check requested features
            if isempty(ops.rsVars)
                featTb = [];
                return
            end
            
            % Make bin edges
            nWin = size(ops.rsWin, 1);
            assert(nWin==1 || nWin==se.numEpochs, "The number of time windows in ops.rsWin must be one or equal the number of trials in se.");
            if nWin == 1
                ops.rsWin = repmat(ops.rsWin, [se.numEpochs 1]);
            end
            tEdges = arrayfun(@(a,b) (a : ops.rsBinSize : b)', ops.rsWin(:,1), ops.rsWin(:,2), 'Uni', false);
            
            % Initialize a table with timestamps
            featTb = table();
            featTb.time = cellfun(@MMath.BinEdges2Centers, tEdges, 'Uni', false);
            
            for i = 1 : numel(ops.rsVars)
                % Unpack struct
                s = ops.rsVars(i);
                tn = s.tableName;
                vn = s.varNames;
                rd = s.readout;
                if isempty(rd)
                    rd = "default";
                end
                fprintf("Resample %s from the '%s' table using %s readout.\n", strjoin(vn, ', '), tn, rd);
                
                % Handle special variables first
                isSpecial = true(size(vn));
                for j = 1 : numel(vn)
                    switch [tn '_' vn{j}]
                        case 'morph_speed'
                            % Derive speed from time0
                            tsTb = se.ResampleTimeSeries(tn, tEdges, [], {'time0'}, ops.rsArgs{:});
                            featTb.speed = cellfun(@(t,t0) gradient(t)./gradient(t0), tsTb.time, tsTb.time0, 'Uni', false);
                        otherwise
                            isSpecial(j) = false;
                    end
                end
                vn = vn(~isSpecial);
                if isempty(vn)
                    continue
                end
                
                % Standard resampling
                k = strcmp(tn, se.tableNames);
                if se.isEventTimesTable(k)
                    if rd == "span"
                        % Resample event spans
                        tsTb = NP.SE.ResampleEventSpans(se, tn, tEdges, [], vn);
                    else
                        % Resample event times
                        tsTb = se.ResampleEventTimes(tn, tEdges, [], vn, 'Normalization', 'countdensity');
                    end
                    featTb = mergeTables(featTb, tsTb, tn);
                    
                elseif se.isTimesSeriesTable(k)
                    if rd == "event"
                        % Resample timeseries without interpolation
                        tsTb = NP.SE.ResampleEventValues(se, tn, tEdges, [], vn);
                    else
                        % Resample timeseries with standard interpolation
                        tsTb = se.ResampleTimeSeries(tn, tEdges, [], vn, ops.rsArgs{:});
                    end
                    featTb = mergeTables(featTb, tsTb, tn);
                    
                elseif se.isEventValuesTable(k)
                    % Resample scalar attributes
                    ev = se.GetTable(tn);
                    ev = ev(:,vn);
                    nSp = cellfun(@numel, featTb.time);
                    for m = 1 : width(ev)
                        ev.(m) = arrayfun(@(x,n) repmat(x, [n 1]), ev.(m), nSp, 'Uni', false);
                    end
                    featTb = mergeTables(featTb, ev, tn);
                    
                else
                    fprintf("The requested table '%s' is invalid.\n", tn);
                end
            end
            
            function tb = mergeTables(tb1, tb2, tb2name)
                vn1 = tb1.Properties.VariableNames;
                vn2 = tb2.Properties.VariableNames;
                tb2 = tb2(:, vn2~="time");
                vn2 = vn2(vn2~="time");
                for n = 1 : numel(vn2)
                    if ismember(vn2{n}, vn1)
                        vn2{n} = [vn2{n} '_' tb2name];
                    end
                end
                tb2.Properties.VariableNames = vn2;
                tb = [tb1 tb2];
            end
        end
        
        function respTb = ResampleResponses(se, ops)
            % Resample spiking responses
            % 
            %   respTb = ResampleResponses(se, ops)
            % 
            % Inputs
            %   se                  A MSessionExplorer object.
            %   ops                 A struct with the following options:
            %   ops.rsWin           A n-by-2 array of time window(s). When n == 1, the same window is used to 
            %                       resample all epochs in se; when n == se.numEpochs, each window is used for 
            %                       the corresponding epoch.
            %   ops.rsBinSize       Bin size or sample size for resampling.
            % 
            % Output
            %   respTb              A timeSeries table of extracted and resampled neural responses.
            % 
            
            % Make bin edges
            tWin = ops.rsWin;
            tEdges = (tWin(1) : ops.rsBinSize : tWin(2))';
            
            % Resampling
            if ismember('spikeRate', se.tableNames)
                tn  ='spikeRate';
                fprintf("Resample responses from '%s' table\n", tn);
                respTb = se.ResampleTimeSeries(tn, tEdges, 'Antialiasing', true, ops.rsArgs{:});
                
            elseif ismember('spikeTime', se.tableNames)
                tn  ='spikeTime';
                fprintf("Resample responses from '%s' table\n", tn);
                respTb = se.ResampleEventTimes(tn, tEdges, 'Normalization', 'countdensity');
                
            else
                error("The se object contains neither spikeRate nor spikeTime table");
            end
        end
        
        function tsTb = ResampleEventSpans(se, tbName, tEdges, rowInd, colInd)
            % Resample event objects to timeseries of event spans
            % 
            %   tsTb = ResampleEventSpans(se, tbName, tEdges, rowInd, colInd)
            % 
            % Input convention is the same as MSessionExplorer.ResampleEventTimes
            
            etTb = se.GetTable(tbName);
            tsTb = table;
            for i = height(etTb) : -1 : 1
                % Timestamps
                tsTb.time{i} = MMath.BinEdges2Centers(tEdges{i});
                
                % Create event masks
                for j = 1 : numel(colInd)
                    vn = colInd{j};
                    if iscell(etTb.(vn))
                        evt = etTb.(vn){i};
                    else
                        evt = etTb.(vn)(i);
                    end
                    if isempty(evt)
                        m = zeros(size(tsTb.time{i}));
                    else
                        m = evt.MaskTimestamps(tsTb.time{i});
                    end
                    tsTb.(vn){i} = double(m);
                end
            end
        end
        
        function tsTb = ResampleEventValues(se, tbName, tEdges, rowInd, colInd)
            % Resample a timeseries table but treat samples as single events. 
            % If multiple events fall in one time bin, mean value will be used.
            % 
            %   tsTb = ResampleEventValues(se, tbName, tEdges, rowInd, colInd)
            % 
            % Input convention is the same as MSessionExplorer.ResampleTimeseries
            
            tsTb = se.GetTable(tbName);
            for i = 1 : height(tsTb)
                % Find which event goes to which time bin
                [N, ~, bin] = histcounts(tsTb.time{i}, tEdges{i});
                
                % Convert bin indices to group indices to be compatible with splitapply
                isOut = ~bin; % out of range bins have indices of zero
                bin(isOut) = [];
                G = findgroups(bin);
                
                % Compute means of event values falling into each time bin
                tCenters = MMath.BinEdges2Centers(tEdges{i});
                for j = 1 : numel(colInd)
                    vn = colInd{j};
                    vEvt = tsTb.(vn){i}(~isOut,:);
                    v = zeros(numel(tCenters), size(vEvt,2));
                    v(N>0,:) = splitapply(@(x) mean(x,1,'omitnan'), vEvt, G);
                    tsTb.(vn){i} = v;
                end
                
                % Replace event times with continous timestamps
                tsTb.time{i} = tCenters;
            end
        end
        
        function featTb = ResampleEventFeatures(se, ops)
            % Resample feature variables
            % 
            %   featTb = ResampleEventFeatures(se, ops)
            % 
            % Inputs
            %   se                  A MSessionExplorer object.
            %   ops                 A struct with the following options:
            %                       
            % Output
            %   featTb              A timeSeries table of extracted and resampled features.
            % 
            
            % Check inputs
            if isempty(ops.rsVars)
                featTb = [];
                return
            end
            
            % Get events and time windows
            evts = se.GetTable(ops.eventVar.tableName).(ops.eventVar.varName);
            tWins = cell(size(evts));
            if isempty(ops.eventVar.readout)
                ops.eventVar.readout = 0;
            end
            for i = 1 : numel(tWins)
                e = evts{i};
                e = e(~isnan(e));
                evts{i} = e;
                if isempty(e)
                    continue
                end
                tWins{i} = [e.GetTfield('tmin') e.GetTfield('tmax')] + ops.eventVar.readout;
            end
            evts = cat(1, evts{:});
            dur = diff(cat(1, tWins{:}),1,2);
            
            % Initialize a table with timestamps
            featTb = table();
            featTb.time = {evts};
            
            for i = 1 : numel(ops.rsVars)
                % Unpack struct
                s = ops.rsVars(i);
                tn = s.tableName;
                vn = s.varNames;
                rd = s.readout;
                if isempty(rd)
                    rd = "default";
                end
                fprintf("Resample %s from the '%s' table using %s readout.\n", strjoin(vn, ', '), tn, rd);
                
                % Standard resampling
                k = strcmp(tn, se.tableNames);
                if se.isEventTimesTable(k)
                    % Resample event times
                    tsTb = se.SliceEventTimes(tn, tWins, [], vn);
                    if rd == "default"
                        rd = "contain";
                    end
                    tsTb = computeFeatVal(tsTb, rd, dur);
                    featTb = mergeTables(featTb, tsTb, tn);
                    
                elseif se.isTimesSeriesTable(k)
                    % Slice timeseries with standard interpolation
                    tsTb = se.SliceTimeSeries(tn, tWins, [], vn);
                    if rd == "default"
                        rd = "average";
                    end
                    tsTb = computeFeatVal(tsTb(:,2:end), rd);
                    featTb = mergeTables(featTb, tsTb, tn);
                    
                elseif se.isEventValuesTable(k)
                    % Resample scalar attributes
                    error("This part is not yet implemented");
                else
                    fprintf("The requested table '%s' is invalid.\n", tn);
                end
            end
            
            function tbOut = computeFeatVal(tbIn, readout, dur)
                tbOut = table;
                for n = 1 : width(tbIn)
                    cn = tbIn.Properties.VariableNames{n};
                    col = tbIn.(n);
                    if readout == "rate"
                        col = cellfun(@(x,d) numel(x(~isnan(x)))/d, col, num2cell(dur), 'Uni', false);
                    elseif readout == "contain"
                        col = cellfun(@(x) ~isnan(x(1)), col, 'Uni', false);
                    elseif readout == "average"
                        col = cellfun(@(x) mean(x,1,'omitnan'), col, 'Uni', false);
                    else
                        error("'%s' is not a valid readout option.", readout);
                    end
                    col = cell2mat(col);
                    tbOut.(cn){1} = col;
                end
            end
            
            function tb = mergeTables(tb1, tb2, tb2name)
                vn1 = tb1.Properties.VariableNames;
                vn2 = tb2.Properties.VariableNames;
                tb2 = tb2(:, vn2~="time");
                vn2 = vn2(vn2~="time");
                for n = 1 : numel(vn2)
                    if ismember(vn2{n}, vn1)
                        vn2{n} = [vn2{n} '_' tb2name];
                    end
                end
                tb2.Properties.VariableNames = vn2;
                tb = [tb1 tb2];
            end
            
        end
        
        function mTb = MeanTimeseries(tsTb, varargin)
            % Compute mean timeseries across epochs
            % 
            %   mTb = MeanTimeseries(tsTb)
            %   mTb = MeanTimeseries(tsTb, ..., 'TrialMask', true(height(tsTb),1))
            %   mTb = MeanTimeseries(tsTb, ..., 'StatFunc', @MMath.MeanStats)
            % 
            % Inputs
            %   tsTb            A timeSeries table in which all epochs have the same number of samples.
            %   'TrialMask'     1) epoch-by-1 logical vector which masks a set of epochs for all variables.
            %                   2) epoch-by-col logical array which can mask a different set of epochs for 
            %                      each variable, including the first column 'time'.
            %                   3) epoch-by-(col-1) array, similar to 2), but masks all variables except 
            %                      for 'time'. Average time will be computed from all epochs.
            % Output
            %   mTb             A timeSeries table. Columns are the same as tsTb.
            %                   The first row stores the mean.
            %                   The second row stores the standard deviation.
            %                   The third row stores the standard error of the mean.
            % 
            % see also MMath.MeanStats
            
            p = inputParser;
            p.addParameter('TrialMask', true(height(tsTb),1), @(x) size(x,1)==height(tsTb));
            p.addParameter('StatFunc', @MMath.MeanStats, @(x) isa(x, 'function_handle'));
            p.parse(varargin{:});
            M = p.Results.TrialMask;
            func = p.Results.StatFunc;
            
            nSp = cellfun(@(x) size(x,1), tsTb.time);
            isEqualSp = isscalar(unique(nSp));
            assert(isEqualSp, "Epochs do not have equal number of sample thus cannot be averaged.");
            
            if iscolumn(M)
                % Propagate epoch mask for all variables
                M = repmat(M(:), [1 width(tsTb)]);
            elseif size(M,2) == width(tsTb)-1
                % Pad mask at the time column if not provided
                M = [true(height(tsTb),1) M];
            end
            
            mTb = table;
            for j = 1 : width(tsTb)
                vn = tsTb.Properties.VariableNames{j};
                arr = cat(3, tsTb.(vn){M(:,j)});
                [mTb.(vn){1}, mTb.(vn){2}, mTb.(vn){3}] = func(arr, 3);
            end
            
            for i = 2 : height(mTb)
                mTb.time{i} = mTb.time{1};
            end
        end
        
    end
    
end

