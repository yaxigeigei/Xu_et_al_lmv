classdef LMV < TaskBaseClass
    %TASKBASECLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Current version of stimuli
        timitListFile = "TIMIT_selected_221011.xlsx";
        dimexListFile = "DIMEX_selected_220808.xlsx";
        brListFile = "BR_selected_240823.xlsx";
        
        % Interval between beeps
        % note: interval between cue2 offset and cue3 onset = postCueDelay(2) + audio stim dur * stimDurFactor + postMentalDelay
        postCueDelays = [0.3 0.5 0];    % delay after each cue offset
        postStimDelayMu = 1;            % fixed delay after audio stim offset and "think" onset
        postStimDelayMin = 0.1;
        postStimDelayMax = 3;
        stimDurFactor = 1.2;
        postMentalDelay = 0.5;
        
        % Beep parameters
        beepFreq = [392 262 330]; % G4 C4 E4
        beepDur = 0.3; % duration of beep
        
        % Trigger
        trigDelay = 0.146;
        trigDur = 0.1;
    end
    
    methods
        function this = LMV()
            % LMV Construct an instance of this class
        end
        
        function Run(this)
            % Run LMV task
            
            %---------------
            % Input Task Options
            %---------------
            
            dlgtitle = 'Task Setup';
            dims = [1 40];
            prompt = { ...
                'Subject ID:', 'NP00'; ...
                'Block #:', '0'; ...
                'Stim set (TIMIT/DIMEX/BR/<filename>):', 'TIMIT'; ...
                'Cue type (beep or voice):', 'beep'; ...
                'Include mental repeat (No: 0; Yes: 1)', '0'; ...
                '# of rounds:', '3'; ...
                'Stereo: 0; Mono + tigger: 1', '1'; ...
                };
            
            setupCacheFile = 'last_task_setup.mat';
            if exist(setupCacheFile, 'file')
                load(setupCacheFile, 'prompt');
            end
            
            answers = inputdlg(prompt(:,1), dlgtitle, dims, prompt(:,2));
            
            if isempty(answers)
                disp('Task did not start');
                return
            else
                prompt(:,2) = answers;
                save(setupCacheFile, 'prompt');
            end
            
            k = 0;
            
            k = k + 1;
            subjId = answers{k};
            
            k = k + 1;
            blockId = ['B' answers{k}];
            
            k = k + 1;
            stimSet = answers{k};
            
            k = k + 1;
            cueType = lower(answers{k});
            if ~any(strcmp(cueType, {'beep', 'voice'}))
                error('Cue type must be ''beep'' or ''voice'', but was ''%s''', cueType);
            end
            
            k = k + 1;
            doMentalRepeat = str2double(answers{k});
            
            k = k + 1;
            nRounds = str2double(answers{k});
            if isnan(nRounds)
                nRounds = Inf;
            end
            
            k = k + 1;
            isTrigger = str2double(answers{k});
            
            
            %---------------
            % Sound Setup
            %---------------
            
            % Initialize Sounddriver
            InitializePsychSound(1);
            
            % Open Psych-Audio port, with the follow arguements
            % (1) [] = default sound device
            % (2) 1 = sound playback only
            % (3) 1 = default level of latency
            % (4) Requested frequency in samples per second
            % (5) 2 = stereo putput
            nrchannels = 2; % number of channels
            try
                freq = 48000; % sampling frequency of the sound device
                pahandle = PsychPortAudio('Open', [], 1, 1, freq, nrchannels);
            catch
                freq = 44100; % sampling frequency of the sound device
                pahandle = PsychPortAudio('Open', [], 1, 1, freq, nrchannels);
            end
            
            % Set the volume to half for this demo
            PsychPortAudio('Volume', pahandle, 0.5);
            
            % How many times to we wish to play each sound
            repetitions = 1;
            
            % Should we wait for the device to really start (1 = yes)
            % INFO: See help PsychPortAudio
            waitForDeviceStart = 1;
            
            % Construct cue objects
            switch cueType
                case 'beep'
                    % Make beep objects
                    cues = arrayfun(@(f) MakeBeep(f, this.beepDur, freq), this.beepFreq, 'Uni', false);
                case 'voice'
                    % Load prompt sounds
                    [promptListen, promptThink, promptSpeak] = this.LoadPrompts(freq);
                    cues = {promptListen, promptThink, promptSpeak};
                otherwise
                    error('Cue type must be ''beep'' or ''voice'', but was ''%s''', cueType);
            end
            
            % Load stimuli
            switch upper(stimSet)
                case 'TIMIT'
                    stimSet = this.timitListFile;
                case 'DIMEX'
                    stimSet = this.dimexListFile;
                case 'BR'
                    stimSet = this.brListFile;
            end
            stimTb = this.LoadStim(stimSet, freq);
            
            
            %---------------
            % Task
            %---------------
            
            % Open log file
            logFolder = fullfile(pwd, 'data');
            if ~exist(logFolder, 'dir')
                mkdir(logFolder);
            end
            logPath = fullfile(logFolder, [subjId '_' blockId '_' char(datetime('now', 'Format', 'yyyy-MM-dd_hh-mm-ss')) '.txt']);
            f = fopen(logPath, 'W');
            
            % Write task parameters
            fprintf(f, 'O|taskName|LMV\n');
            fprintf(f, 'O|cueType|%s\n', cueType);
            fprintf(f, 'O|beepDur|%g\n', this.beepDur);
            fprintf(f, 'O|beepFreq|%g|%g|%g\n', this.beepFreq(:));
            fprintf(f, 'O|postCueDelay|%g|%g|%g\n', this.postCueDelays(:));
            fprintf(f, 'O|postStimDelayMu|%g\n', this.postStimDelayMu);
            fprintf(f, 'O|postStimDelayMin|%g\n', this.postStimDelayMin);
            fprintf(f, 'O|postStimDelayMax|%g\n', this.postStimDelayMax);
            fprintf(f, 'O|doMentalRepeat|%g\n', doMentalRepeat);
            fprintf(f, 'O|stimDurFactor|%g\n', this.stimDurFactor);
            fprintf(f, 'O|postMentalDelay|%g\n', this.postMentalDelay);
            
            
            % The avaliable keys to press
            KbName('UnifyKeyNames');
            nextKey = KbName('RightArrow');
            escapeKey = KbName('ESCAPE');
            
            % Runtime variables
            trialNum = 0;
            exitTask = false;
            taskStartTime = GetSecs;
            
            % Wait for keypress to start the first trial
            disp("Press right arrow key to start a trial, or ESC to end the task.");
            [secs, keyCode, deltaSecs] = KbStrokeWait;
            while ~any(keyCode([nextKey escapeKey]))
                [secs, keyCode, deltaSecs] = KbStrokeWait;
            end
            if keyCode(escapeKey)
                PsychPortAudio('Close', pahandle); % close the audio device
                return
            end
            
            while exitTask == false && trialNum < nRounds*height(stimTb)
                
                %---------------
                % Preparation
                %---------------
                
                % Get stimulus data
                trialNum = trialNum + 1;
                iStim = mod(trialNum-1, height(stimTb)) + 1;
                stimId = stimTb.id{iStim};
                stimText = stimTb.text{iStim};
                stim = stimTb.audio{iStim}';
                if isTrigger
                    trig = this.MakeTriggerWave(stim, freq, this.trigDur, this.trigDelay);
                    stimAudio = PsychPortAudio('CreateBuffer', [], [stim; trig]);
                    cueAudio = cellfun(@(x) [x; zeros(size(x))], cues, 'Uni', false);
                else
                    stimAudio = PsychPortAudio('CreateBuffer', [], [stim; stim]);
                    cueAudio = cellfun(@(x) [x; x], cues, 'Uni', false);
                end
                stimDur = numel(stim) / freq;
                
                % Log trial info
                startTime = GetSecs + 0.1; % give some extra time before playing cue
                trialStartTime = startTime;
                fprintf(f, 'I|trial|%.3f|%i\n', startTime, trialNum);
                fprintf(f, 'I|stimId|%.3f|%s\n', startTime, stimId);
                fprintf(f, 'I|stim|%.3f|%s\n', startTime, stimText);
                
                % Print trial info
                fprintf('Trial %i\n', trialNum);
                fprintf('Stim #%i: %s\n', iStim, stimText);
                
                
                %---------------
                % Listening
                %---------------
                
                % Fill the audio playback buffer with the audio data, doubled for stereo presentation
                PsychPortAudio('FillBuffer', pahandle, cueAudio{1});
                
                % Start audio playback
                PsychPortAudio('Start', pahandle, repetitions, startTime, waitForDeviceStart);
                fprintf(f, 'I|cue1On|%.3f\n', startTime);
                
                % Wait for the beep to end
                [~, ~, ~, estStopTime] = PsychPortAudio('Stop', pahandle, 1, 1);
                fprintf(f, 'I|cue1Off|%.3f\n', estStopTime);
                
                % Compute new start time for follow-up beep, beepPauseTime after end of previous one
                startTime = estStopTime + this.postCueDelays(1);
                
                % Fill the audio playback buffer with the audio data, doubled for stereo presentation
                PsychPortAudio('FillBuffer', pahandle, stimAudio);
                
                % Start audio playback
                PsychPortAudio('Start', pahandle, repetitions, startTime, waitForDeviceStart);
                fprintf(f, 'I|stimOn|%.3f\n', startTime);
                
                % Wait for stop of playback
                [~, ~, ~, estStopTime] = PsychPortAudio('Stop', pahandle, 1, 1);
                fprintf(f, 'I|stimOff|%.3f\n', estStopTime);
                
                % Compute new start time for follow-up beep, beepPauseTime after end of previous one
                postStimDelay = this.SampleExp(this.postStimDelayMu, this.postStimDelayMin, this.postStimDelayMax);
                fprintf(f, 'I|postStimDelay|%.3f|%.3f\n', estStopTime, postStimDelay);
                startTime = estStopTime + postStimDelay;
                
                
                %---------------
                % Mental Repeat
                %---------------
                
                if doMentalRepeat
                    % Fill the audio playback buffer with the audio data, doubled for stereo presentation
                    PsychPortAudio('FillBuffer', pahandle, cueAudio{2});
                    
                    % Start the beep
                    PsychPortAudio('Start', pahandle, repetitions, startTime, waitForDeviceStart);
                    fprintf(f, 'I|cue2On|%.3f\n', startTime);
                    
                    % Wait for the beep to end
                    [~, ~, ~, estStopTime] = PsychPortAudio('Stop', pahandle, 1, 1);
                    fprintf(f, 'I|cue2Off|%.3f\n', estStopTime);
                    
                    % Compute new start time for follow-up beep, beepPauseTime after end of previous one
                    startTime = estStopTime + this.postCueDelays(2) + stimDur * this.stimDurFactor + this.postMentalDelay;
                end
                
                
                %---------------
                % Vocal Repeat
                %---------------
                
                % Fill the audio playback buffer with the audio data, doubled for stereo presentation
                PsychPortAudio('FillBuffer', pahandle, cueAudio{3});
                
                % Start the beep
                PsychPortAudio('Start', pahandle, repetitions, startTime, waitForDeviceStart);
                fprintf(f, 'I|cue3On|%.3f\n', startTime);
                
                % Wait for the beep to end
                [~, ~, ~, estStopTime] = PsychPortAudio('Stop', pahandle, 1, 1);
                fprintf(f, 'I|cue3Off|%.3f\n', estStopTime);
                
                
                %---------------
                % End
                %---------------
                
                % Wait for a keyboard button to start next trial
                [secs, keyCode, deltaSecs] = KbStrokeWait;
                while ~any(keyCode([nextKey escapeKey]))
                    [secs, keyCode, deltaSecs] = KbStrokeWait;
                end
                
                % Print trial info
                fprintf('Trial/total time: %.1f / %.1f s\n\n', GetSecs-trialStartTime, GetSecs-taskStartTime);
                
                % Depending on the button press, either continue or exit the task
                if keyCode(escapeKey)
                    break
                end
            end
            
            % Close the audio device
            PsychPortAudio('Close', pahandle);
            
            % Close log file
            fclose('all');
            
            % Convert log file
            tb = this.TXT2CSV(logPath);
            assignin('base', 'trialTb', tb);
        end
        
        function [L, T, S] = LoadPrompts(this, freq)
            % Load voice prompts of "listen", "think", and "speak"
            %
            %   [L, T, S] = LoadPrompts(freq)
            %
            
            promFolder = fullfile(pwd, 'prompts');
            
            [w, Fs] = audioread(fullfile(promFolder, 'listen.wav'));
            w = this.AdaptAudio(w, Fs, freq, 4);
            L = w(:,1)';
            
            [w, Fs] = audioread(fullfile(promFolder, 'think.wav'));
            w = this.AdaptAudio(w, Fs, freq, 4);
            T = w(:,1)';
            
            [w, Fs] = audioread(fullfile(promFolder, 'speak.wav'));
            w = this.AdaptAudio(w, Fs, freq, 4);
            S = w(:,1)';
        end
        
        function tb = LoadStim(this, stimListFile, toFs)
            % Load selected stimuli
            %
            %   tb = LoadStim(stimListFile, toFs)
            %
            % Inputs
            %   stimListFile    Path of the stimulus list spreadsheet.
            %   toFs            The desired sampling frequency of the audio waveform.
            %
            % Output
            %   tb              Table of stimuli with at least the following columns:
            %                   id        Stim identifier, e.g. 'fbcg1_si2242'.
            %                   text      Text of the stim, e.g. 'Twenty-two or twenty-three.'
            %                   audio     Stim wavefrom.
            %
            %                   Additional columns may be included:
            %                   name      Txt file of stim.
            %                   folder    Folder path of stim.
            % 
            
            wavSearch = struct2table(dir(fullfile(pwd, "**", "*.wav")));
            
            tb = readtable(stimListFile);
            
            for i = 1 : height(tb)
                % Find audio file
                isStim = tb.id{i}+".wav" == wavSearch.name;
                wavFile = fullfile(wavSearch.folder{isStim}, wavSearch.name{isStim});
                
                % Read audio
                [w, Fs] = audioread(wavFile);
                w = this.AdaptAudio(w, Fs, toFs, 4);
                tb.audio{i} = w;
            end
        end
        
    end
    
    methods(Static)
        function tb = TXT2CSV(txtFiles, csvDir)
            % Convert txt log of LMV task to CSV file
            %
            %   tb = TXT2CSV()
            %   tb = TXT2CSV(txtFiles)
            %   tb = TXT2CSV(txtFiles, csvDir)
            %
            % Input
            %   txtFiles        The txt file(s) to convert. A file selection window will show up if this
            %                   argument is empty []. When multiple files are selected, trials will be
            %                   concatenated and re-numbered continuously. This is useful when the task
            %                   was run multiple times during a block and generated multiple log files.
            %   csvDir          The folder path to save the CSV file in. If not provided, the CSV file
            %                   will be saved at '<preproc_root>/recId/labels/lmv' folder. The folder will
            %                   be created if it doesn't exist.
            % Output
            %   tb              The table being saved.
            %
            
            % Find txt log files
            if ~exist('txtFiles', 'var') || isempty(txtFiles)
                txtFiles = MBrowse.Files([], "Please select one or more LMV log files", {'*.txt'});
            end
            if isempty(txtFiles)
                return
            end
            
            % Read log data
            txtFiles = cellstr(txtFiles);
            tbs = cell(size(txtFiles));
            
            for k = 1 : numel(txtFiles)
                % Load mat file
                sLog = Satellites.Import(txtFiles{k}, ...
                    'TagDelimiter', '|', ...
                    'MsgDelimiter', '|', ...
                    'ForceNumericValues', false, ...
                    'DelimiterEvent', 'trial', ...
                    'TimeScaling', 1);
                
                et = sLog.timeTable;
                ev = sLog.valueTable;
                rt = sLog.episodeRefTime;
                
                % Add selected info to table
                tb = table;
                tb.trial_num = ev.trial;
                tb.stim_id = ev.stimId;
                tb.label = ev.stim;
                tb.ref_time = rt;           % reference time
                tb.cue1_onset = et.cue1On;  % cue for listening
                tb.cue3_onset = et.cue3On;  % cue for vocal repeat
                
                % Remove the pre-task episode
                tbs{k} = tb(2:end,:);
            end
            
            
            % Countinue the trial number count
            lastTrial = 0;
            for k = 1 : numel(tbs)
                tbs{k}.trial_num = tbs{k}.trial_num + lastTrial;
                lastTrial = tbs{k}.trial_num(end);
            end
            
            % Concatenate tables
            tb = cat(1, tbs{:});
            
            % Convert to global time
            tb.ref_time = tb.ref_time - tb.ref_time(1);
            tb{:,{'cue1_onset','cue3_onset'}} = tb{:,{'cue1_onset','cue3_onset'}} + tb.ref_time;
            
            
            % Save table as CSV file
            [logFolders, logNames] = cellfun(@fileparts, txtFiles, 'Uni', false);
            
            if ~exist('csvDir', 'var')
                %     recId = regexp(logNames{1}, '^NP[0-9]+_B[0-9]+', 'match', 'once');
                %     csvDir = fullfile(NP.Data.GetPreprocRoot, recId, 'labels', 'lmv');
                csvDir = logFolders{1};
            end
            if ~isempty(csvDir) && ~exist(csvDir, 'dir')
                mkdir(csvDir);
            end
            
            if numel(logNames) == 1
                csvName = logNames{1} + ".csv";
            else
                csvName = logNames{1} + "_etc.csv";
            end
            csvFile = fullfile(csvDir, csvName);
            
            writetable(tb, csvFile);
        end
    end
end

