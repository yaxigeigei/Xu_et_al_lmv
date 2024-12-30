classdef SEMaker
    % SEMaker collects and organize preprocessed data into a MSessionExplorer object
    % Object constructor
    %
    %   M = SEMaker(recId)
    %   M = SEMaker(subjectId, blockId)
    %   M = SEMaker(se)
    %
    % Inputs
    %   recId           Recording ID, e.g. 'NP38_B6'
    %   subjectId       Subject ID, e.g. 'NP38'
    %   blockId         Block ID, e.g. 'B6'
    %   se              An existing se object
    %
    
    properties
        ppRoot                  % root directory of preprocessing
        se MSessionExplorer     % the MSessionExplorer object
        isOverwrite = false     % by default, whether to overwrite existing data
        isDry = false           % whether to run the maker without actually changing data in se (i.e. dry run)
    end
    
    properties(Dependent)
        recId
        subjectId
        blockId
    end
    
    methods
        function val = get.recId(this)
            val = this.se.userData.expInfo.recId;
        end
        function val = get.subjectId(this)
            val = this.se.userData.expInfo.subjectId;
        end
        function val = get.blockId(this)
            val = this.se.userData.expInfo.blockId;
        end
    end
    
    methods
        function this = SEMaker(varargin)
            %SEMAKER Construct an instance of this class
            
            % 
            this.ppRoot = NP.Data.GetPreprocRoot();
            
            if isa(varargin{1}, 'MSessionExplorer')
                % Use the given se object
                seIn = varargin{1};
                assert(isfield(seIn.userData, 'expInfo'), 'The given se object must have valid expInfo');
                this.se = seIn;
            else
                % Create a new se and fill in
                this.se = MSessionExplorer();
                this.AddExperimentInfo(varargin{:});
            end
        end
        
        function AddExperimentInfo(this, varargin)
            % Add expInfo to se.userData
            %
            %   M.AddExperimentInfo(recId)
            %   M.AddExperimentInfo(subjectId, blockId)
            %
            % Inputs
            %   recId           Recording ID, e.g. 'NP38_B6'
            %   subjectId       Subject ID, e.g. 'NP38'
            %   blockId         Block ID, e.g. 'B6'
            %
            
            % Get IDs
            if numel(varargin) == 1
                recordId = varargin{1};
                recIdParts = strsplit(recordId, '_');
                [subjId, blkId] = recIdParts{1:2};
            else
                [subjId, blkId] = varargin{:};
                recordId = strjoin(varargin, '_');
            end
            
            % Make a struct
            s = struct;
            s.subjectId = char(subjId);
            s.blockId = char(blkId);
            s.recId = char(recordId);
            
            % Save to SE
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.userData.expInfo = s;
        end
        
        function AddNI(this, varargin)
            % Add analog data to ni table. Extract digital events and add to niTimes table. Add metadata to se.userData.niMeta.
            % Data will be organized in a single epoch with zero as the reference time.
            %
            %   M.AddNI()
            %   M.AddNI(isForce)
            %   M.AddNI(..., 'AnalogChanNames', {'sync', 'mic', 'speaker1', 'speaker2'})
            %   M.AddNI(..., 'DigitalChanNames', {'sync', 'trig1', 'trig2'})
            %
            % Inputs
            %   isForce                 Whether or not to force overwrite existing data. Default is false.
            %   'AnalogChanNames'       A list of analog channel names. The number of names determines the number of
            %                           channels to read (starting from the first channel).
            %   'DigitalChanNames'      A list of digital channel names. The number of names determines the number of
            %                           channels to read (starting from the first channel).
            %
            
            fprintf('\n');
            
            % Parse options
            p = inputParser();
            p.addOptional('isForce', this.isOverwrite, @islogical);
            p.addParameter('AnalogChanNames', {'sync', 'mic', 'speaker1', 'speaker2'}, @iscellstr);
            p.addParameter('DigitalChanNames', {'sync', 'trig1', 'trig2'}, @iscellstr);
            p.parse(varargin{:});
            isForce = p.Results.isForce;
            aiChanNames = p.Results.AnalogChanNames;
            diChanNames = p.Results.DigitalChanNames;
            
            % Check if all the data already exist
            if all(ismember({'ni', 'niTime'}, this.se.tableNames)) && ~isForce
                disp('Will not overwrite the existing NIDQ data');
                return
            end
            
            % Search for a nidq.bin file
            binPath = this.IFindSglxFiles('ni');
            if isempty(binPath)
                disp("Skip NIDQ since the required file(s) are missing.");
                return
            end
            
            % Read NIDQ data
            disp('Reading NIDQ data ...');
            [meta, ana, dig, t] = MSpikeGLX.ReadNI(binPath);
            disp('Finished reading.');
            
            % Make a table of analog timeseries
            nAna = min(size(ana,2), numel(aiChanNames));
            if nAna > 0
                ana = ana(:, 1:nAna);
                niTb = MSessionExplorer.MakeTimeSeriesTable(t, num2cell(ana,1), 'VariableNames', aiChanNames, 'Verbose', false);
            else
                niTb = [];
            end
            
            % Make a table of onset and offset times of the digital signals
            nDig = min(size(dig,2), numel(diChanNames));
            if nDig > 0
                dig = dig(:, 1:nDig);
                nitTb = table;
                for i = 1 : numel(diChanNames)
                    win = MMath.Logical2Bounds(dig(:,i));
                    win = win / str2double(meta.niSampRate);
                    nitTb.([diChanNames{i} 'On']) = {win(:,1)};
                    nitTb.([diChanNames{i} 'Off']) = {win(:,2)};
                end
            else
                nitTb = [];
            end
            
            % Add tables and metadata
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.userData.niMeta = meta;
            disp("Added metadata in se.userData.niMeta.");
            if ~isempty(niTb)
                this.se.SetTable('ni', niTb, 'timeSeries', 0);
                disp("Added 'ni' table to se.");
            end
            if ~isempty(nitTb)
                this.se.SetTable('niTime', nitTb, 'eventTimes', 0);
                disp("Added 'niTime' table to se.");
            end
        end
        
        function ReplaceMicAudio(this)
            % Replace the mic timeseries in the 'ni' table with the denoised version.
            % This function always overwrites the existing data.
            
            fprintf('\n');
            disp("Replacing 'mic' signal in 'ni' table with the denoised version.");
            if ~ismember('ni', this.se.tableNames)
                warning("Skip replacement because the 'ni' table is missing.");
            end
            
            % Read denoised mic WAV file
            expId = this.recId;
            audioFolder = this.IGetModulePath("audio_files");
            dnMicFile = fullfile(audioFolder, [expId '_mic_denoised.wav']);
            if ~exist(dnMicFile, 'file')
                disp("Skip replacement because denoised audio file does not exist.");
                return
            end
            disp("Read denoised mic audio file.");
            w = single(audioread(dnMicFile));
            
            % Check and replace data
            mic = this.se.GetColumn('ni', 'mic');
            if numel(w) ~= numel(mic{1})
                warning("Original mic audio has %i samples but the denoised has %i samples. Only replacing the first %i.", ...
                    numel(mic{1}), numel(w), numel(w));
                mic{1}(1:numel(w)) = w;
            else
                mic{1} = w;
            end
            
            % Apply change to se
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.SetColumn('ni', 'mic', mic);
            disp("Updated 'mic' data in 'ni' table.");
        end
        
        function AddLFP(this, varargin)
            % Add LFP table and metadata. Data will be organized in a single epoch with zero as the reference time.
            % 
            %   M.AddLFP()
            %   M.AddLFP(isForce)
            %   M.AddLFP(..., 'ProbeNum', [])
            %
            % Inputs
            %   isForce         Whether or not to force overwrite existing data. Default is false.
            %   'ProbeNum'     	The imec stream number, e.g. [0 1 2 ...] indicating imec0, imec1, imec2... 
            %                   Default is [] which will include all probes.
            % 
            
            fprintf('\n');
            
            % Parse options
            p = inputParser();
            p.addOptional('isForce', this.isOverwrite, @islogical);
            p.addParameter('ProbeNum', [], @isnumeric);
            p.parse(varargin{:});
            isForce = p.Results.isForce;
            probeNums = p.Results.ProbeNum;
            
            % Check if all the data already exist
            if ismember('LFP', this.se.tableNames) && ~isForce
                disp('Will not overwrite the existing LFP data');
                return
            end
            
            % Locate LFP data by searching for the lf.bin file
            [lfPaths, probeNames] = this.IFindSglxFiles('lf');
            if ~isempty(probeNums)
                m = ismember(probeNames, "imec"+probeNums);
                lfPaths = lfPaths(m);
                probeNames = probeNames(m);
            end
            if isempty(lfPaths)
                disp("Skip LFP since the required file(s) are missing.");
                return
            end
            
            disp("Adding 'LFP' table");
            lfpTb = table();
            
            for i = 1 : numel(lfPaths)
                % Read LFP data
                disp("Reading lf.meta and lf.bin ...");
                [meta, V] = MSpikeGLX.ReadLFP(lfPaths(i));
                meta = MUtil.EvalFields(meta);
                disp("Finished reading.");
                
                % Read channel table
                chanMapPath = 'NP1_NHP_HalfCol_kilosortChanMap.mat'; % -> should not be hardcoded
                disp("Read channel map from " + chanMapPath);
                [meta.chanTb, meta.chanTbKS] = MKilosort2.LoadChanMap2Table(chanMapPath);
                
                % Remove unconnected (e.g. 193th reference) channels
                disp("Remove the reference channel.");
                V = V(:, meta.chanTb.isConnected);
                
                % Sort channels by the order in chanTbKS (which is by descending distance to tip)
                disp("Sort channels by descending distance to tip.");
                V = V(:, meta.chanTbKS.sortInd);
                
                % Downsample from 2500 to 400 Hz
                fromFs = meta.imSampRate;
                toFs = fromFs * 400/2500; % this handles non-integer fs
                fprintf('Downsampling LFP from %g to %gHz ...\n', fromFs, toFs);
                V = resample(double(V), round(toFs), round(fromFs));
                t = ((1:size(V,1))-0.5) / toFs;
                
                % Make LFP timeseries table
                if i == 1
                    lfpTb.time{1} = t';
                else
                    V = interp1(t', V, lfpTb.time{1}); % unify sampling to the first probe
                end
                lfpTb.(probeNames(i)){1} = single(V);
                
                % Add metadata struct
                lfMeta(i) = meta;
            end
            
            % Add table and metadata
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.userData.lfMeta = lfMeta;
            this.se.SetTable('LFP', lfpTb, 'timeSeries', 0);
            disp("Added 'LFP' table to se, and metadata to se.userData.lfMeta");
        end
        
        function AddAPMeta(this, varargin)
            % Add AP metadata
            % 
            %   M.AddAPMeta()
            %   M.AddAPMeta(isForce)
            %   M.AddAPMeta(..., 'ProbeNum', [])
            %
            % Inputs
            %   isForce         Whether or not to force overwrite existing data. Default is false.
            %   'ProbeNum'     	The imec stream number, e.g. [0 1 2 ...] indicating imec0, imec1, imec2... 
            %                   Default is [] which will include all probes.
            % 
            
            fprintf('\n');
            
            % Parse options
            p = inputParser();
            p.addOptional('isForce', this.isOverwrite, @islogical);
            p.addParameter('ProbeNum', [], @isnumeric);
            p.parse(varargin{:});
            isForce = p.Results.isForce;
            probeNums = p.Results.ProbeNum;
            
            % Check if the data already exist
            if isfield(this.se.userData, 'apMeta') && ~isForce
                disp('Will not overwrite the existing AP metadata');
                return
            end
            
            % Locate LFP data by searching for the lf.bin file
            [apPaths, probeNames] = this.IFindSglxFiles('ap');
            if ~isempty(probeNums)
                m = ismember(probeNames, "imec"+probeNums);
                apPaths = apPaths(m);
                probeNames = probeNames(m);
            end
            if isempty(apPaths)
                disp("Skip AP metadata since the required file(s) are missing.");
                return
            end
            
            % Read AP metadata
            disp("Reading ap.meta file(s)");
            for i = numel(apPaths) : -1 : 1
                apMeta(i) = MSpikeGLX.ReadMetaEZ(apPaths(i));
            end
            apMeta = MUtil.EvalFields(apMeta);
            disp("Finished reading");
            
            % Add table and metadata
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.userData.apMeta = apMeta;
            disp("Added metadata to se.userData.apMeta");
        end
        
        function AddSortingResult(this, varargin)
            % Add spikeTime table. Data will be organized in a single epoch with zero as the reference time.
            %
            %   M.AddSortingResult()
            %   M.AddSortingResult(isForce)
            %   M.AddSortingResult(..., 'SelectionPrompt', true)
            %
            % Inputs
            %   isForce             Whether or not to force overwrite existing data. Default is false.
            %   'SelectionPrompt'   Whether or not to prompt a listbox for users to select all or a subset of 
            %                       Kilosort outputs. Default is true. If set to false, all outputs will be 
            %                       imported (treated as different probes).
            % 
            
            fprintf('\n');
            
            % Parse options
            p = inputParser();
            p.addOptional('isForce', this.isOverwrite, @islogical);
            p.addParameter('SelectionPrompt', true, @islogical);
            p.addParameter('ProbeNum', [], @isnumeric);
            p.parse(varargin{:});
            isForce = p.Results.isForce;
            isPrompt = p.Results.SelectionPrompt;
            probeNums = p.Results.ProbeNum;
            
            % Check if all the data already exist
            if ismember('spikeTime', this.se.tableNames) && ~isForce
                disp('Will not overwrite the existing spikeTime data');
                return
            end
            
            % Searching for the kilosort output folder
            disp('Search for the kilosort output folders');
            ksRoot = this.IGetModulePath("kilosort");
            ksSearch = MBrowse.Dir2Table(fullfile(ksRoot, "*", "*", "cluster_group.tsv"));
            if height(ksSearch) == 0
                disp("Skip sorting results since no Kilosort output can be found.");
                return
            end
            ksFolders = ksSearch.folder;
            disp(ksFolders);
            
            % Allow user to select the result to import when more than one
            if isPrompt
                [ind, isSelect] = listdlg('ListString', ksFolders, ...
                    'PromptString', 'Select kilosort output folders', ...
                    'SelectionMode', 'multiple', ...
                    'ListSize', [800 100]);
                if ~isSelect
                    return
                end
                ksFolders = ksFolders(ind);
            end
            
            stTb = table();
            
            for k = 1 : numel(ksFolders)
                % Import Kilosort/Phy output
                disp("Reading Kilosort/Phy output ...");
                sr = MTracer.KilosortResult();
                sr.ImportData(ksFolders{k}, fullfile(ksFolders{k}, 'chanMap.mat'));
                disp("Finished reading.");
                sr.ComputeAll();
                
                % Remove bulky data
                spkTb = sr.spkTb;
                col2rm = {'pcWeights', 'pcFeatures'};
                for i = 1 : numel(col2rm)
                    if ismember(col2rm{i}, spkTb.Properties.VariableNames)
                        spkTb.pcWeights = [];
                    end
                end
                tempTb = sr.tempTb;
                col2rm = {'pcChanInd'};
                for i = 1 : numel(col2rm)
                    if ismember(col2rm{i}, tempTb.Properties.VariableNames)
                        tempTb.pcWeights = [];
                    end
                end
                
                % Remove noise clusters
                clusTb = sr.clusTb;
                clusTb = clusTb(~strcmp(clusTb.group, 'noise'), :);
                spkTb = spkTb(ismember(spkTb.clusId, clusTb.clusId), :);
                
                % Increment cluster IDs by 1e4 after each probe
                inc = (k-1)*1e4;
                clusTb.clusId = clusTb.clusId + inc;
                spkTb.clusId = spkTb.clusId + inc;
                
                % Make spike time table
                if isnan(spkTb.timeSecAdj)
                    st = spkTb.timeSec;
                else
                    st = spkTb.timeSecAdj;
                end
                
                clusTb = sortrows(clusTb, 'depth', 'descend');
                uid = clusTb.clusId;
                ust = cell(1, numel(uid));
                for i = 1 : numel(uid)
                    ust{i} = st(spkTb.clusId == uid(i));
                end
                
                % Append columns to the table
                stTb = [stTb MSessionExplorer.MakeEventTimesTable(ust, 'VariableNames', "u" + uid, 'Uni', false)];
                
                % Put sr properties in a struct
                meta.ksDir = ksFolders{k};
                meta.ops = sr.rezLite.ops;
                meta.chanMapFile = sr.chanMapFile;
                meta.chanMapName = sr.chanMapName;
                meta.chanMapTb = sr.chanTb;
                meta.clusTb = clusTb;
                meta.spkTb = spkTb;
                meta.tempTb = tempTb;
                ksMeta(k) = meta;
            end
            
            % Add table and metadata
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.userData.ksMeta = ksMeta;
            this.se.SetTable('spikeTime', stTb, 'eventTimes', 0);
            disp("Added 'spikeTime' table to se, and metadata in se.userData.ksMeta");
        end

        function AddSortingResult_SHJ(this, varargin)
            % Add spikeTime table. Data will be organized in a single epoch with zero as the reference time.
            %
            %   M.AddSortingResult()
            %   M.AddSortingResult(isForce)
            %   M.AddSortingResult(..., 'SelectionPrompt', true)
            %
            % Inputs
            %   isForce             Whether or not to force overwrite existing data. Default is false.
            %   'SelectionPrompt'   Whether or not to prompt a listbox for users to select all or a subset of 
            %                       Kilosort outputs. Default is true. If set to false, all outputs will be 
            %                       imported (treated as different probes).
            % 
            
            fprintf('\n');
            
            % Parse options
            p = inputParser();
            p.addOptional('isForce', this.isOverwrite, @islogical);
            p.addParameter('SelectionPrompt', true, @islogical);
            p.addParameter('ProbeNum', [], @isnumeric);
            p.parse(varargin{:});
            isForce = p.Results.isForce;
            isPrompt = p.Results.SelectionPrompt;
            probeNums = p.Results.ProbeNum;
            
            % Check if all the data already exist
            if ismember('spikeTime', this.se.tableNames) && ~isForce
                disp('Will not overwrite the existing spikeTime data');
                return
            end
            
            % Searching for the kilosort output folder
            disp('Search for the kilosort output folders');

            % shj's preproc root is different
            ksRoot = fullfile('/data_store2/data_share/hervey-jumper_neuropixels/preproc', this.recId, 'kilosort');
            
            ksSearch = MBrowse.Dir2Table(fullfile(ksRoot, "*", "*", "cluster_group.tsv"));
            if height(ksSearch) == 0
                disp("Skip sorting results since no Kilosort output can be found.");
                return
            end
            ksFolders = ksSearch.folder;
            disp(ksFolders);
            
            % Allow user to select the result to import when more than one
            if isPrompt
                [ind, isSelect] = listdlg('ListString', ksFolders, ...
                    'PromptString', 'Select kilosort output folders', ...
                    'SelectionMode', 'multiple', ...
                    'ListSize', [800 100]);
                if ~isSelect
                    return
                end
                ksFolders = ksFolders(ind);
            end
            
            stTb = table();
            
            for k = 1 : numel(ksFolders)
                % Import Kilosort/Phy output
                disp("Reading Kilosort/Phy output ...");
                sr = MTracer.KilosortResult();
                sr.ImportData(ksFolders{k}, fullfile(ksFolders{k}, 'chanMap.mat'));
                disp("Finished reading.");
                sr.ComputeAll();
                
                % Remove bulky data
                spkTb = sr.spkTb;
                col2rm = {'pcWeights', 'pcFeatures'};
                for i = 1 : numel(col2rm)
                    if ismember(col2rm{i}, spkTb.Properties.VariableNames)
                        spkTb.pcWeights = [];
                    end
                end
                tempTb = sr.tempTb;
                col2rm = {'pcChanInd'};
                for i = 1 : numel(col2rm)
                    if ismember(col2rm{i}, tempTb.Properties.VariableNames)
                        tempTb.pcWeights = [];
                    end
                end
                
                % Remove noise clusters
                clusTb = sr.clusTb;
                clusTb = clusTb(~strcmp(clusTb.group, 'noise'), :);
                spkTb = spkTb(ismember(spkTb.clusId, clusTb.clusId), :);
                
                % Increment cluster IDs by 1e4 after each probe
                inc = (k-1)*1e4;
                clusTb.clusId = clusTb.clusId + inc;
                spkTb.clusId = spkTb.clusId + inc;
                
                % Make spike time table
                if isnan(spkTb.timeSecAdj)
                    st = spkTb.timeSec;
                else
                    st = spkTb.timeSecAdj;
                end
                
                clusTb = sortrows(clusTb, 'depth', 'descend');
                uid = clusTb.clusId;
                ust = cell(1, numel(uid));
                for i = 1 : numel(uid)
                    ust{i} = st(spkTb.clusId == uid(i));
                end
                
                % Append columns to the table
                stTb = [stTb MSessionExplorer.MakeEventTimesTable(ust, 'VariableNames', "u" + uid, 'Uni', false)];
                
                % Put sr properties in a struct
                meta.ksDir = ksFolders{k};
                meta.ops = sr.rezLite.ops;
                meta.chanMapFile = sr.chanMapFile;
                meta.chanMapName = sr.chanMapName;
                meta.chanMapTb = sr.chanTb;
                meta.clusTb = clusTb;
                meta.spkTb = spkTb;
                meta.tempTb = tempTb;
                ksMeta(k) = meta;
            end
            
            % Add table and metadata
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.userData.ksMeta = ksMeta;
            this.se.SetTable('spikeTime', stTb, 'eventTimes', 0);
            disp("Added 'spikeTime' table to se, and metadata in se.userData.ksMeta");
        end
        
        function AddTaskEvents(this)
            % Create task events as MSessionExplorer.Event objects and add them to the 'taskTime' table.
            % Numeric onset and offset times will also be added for compatibility and convenience.
            % Data will be organized in a single epoch with zero as the reference time.
            % 
            %   AddTaskEvents()
            % 
            
            fprintf('\n');
            
            % Read labels
            labelFolder = this.IGetModulePath("labels");
            labelFile = fullfile(labelFolder, 'combined_task_labels.csv');
            if exist(labelFile, 'file')
                lbTb = readtable(labelFile, 'Delimiter', ',');
            else
                disp("Skip task events since combined_task_labels.csv is missing.");
                return
            end
            
            % Round times to millisecond
            lbTb.OnsetTime = round(lbTb.OnsetTime, 3);
            lbTb.OffsetTime = round(lbTb.OffsetTime, 3);
            
            % Make NP.TGEvent objects
            evt = MSessionExplorer.Event(lbTb.OnsetTime);
            evt = evt.SetTfield('tOn', lbTb.OnsetTime);
            evt = evt.SetTfield('tOff', lbTb.OffsetTime);
            evt = evt.SetVfield('task', string(lbTb.Task));
            evt = evt.SetVfield('source', string(lbTb.Source));
            evt = evt.SetVfield('type', string(lbTb.Type));
            evt = evt.SetVfield('labelId', lbTb.(1));
            
            % taskTime table
            if ismember('taskTime', this.se.tableNames)
                ttTb = this.se.GetTable('taskTime');
                disp('Use the existing taskTime table.');
            else
                ttTb = table;
                disp('Create a new taskTime table.');
            end
            
            % Add events by Type (change to Source in the future)
            evtNames = string(unique(lbTb.Type));
            for n = evtNames(:)'
                isEvt = strcmp(lbTb.Type, n);
                ttTb.(n) = {evt(isEvt)};
                ttTb.(n + "On") = {lbTb.OnsetTime(isEvt)};
                ttTb.(n + "Off") = {lbTb.OffsetTime(isEvt)};
            end
            
            % Save taskTime table to SE
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.SetTable('taskTime', ttTb, 'eventTimes', 0);
            disp("Added task events to the 'taskTime' table.");
        end
        
        function AddReadingEvents(this)
            % Create reading events as NP.TGEvent objects and add them to the 'taskTime' table.
            % Numeric onset and offset times will also be added for compatibility and convenience.
            % Data will be organized in a single epoch with zero as the reference time.
            % Same as AddTaskEvents but with a different label file to keep the two separate.
            % 
            %   AddReadingEvents()
            % 
            
            fprintf('\n');
            
            % Read labels
            labelFolder = this.IGetModulePath("labels");
            labelFile = fullfile(labelFolder, 'combined_reading_labels.csv');
            if exist(labelFile, 'file')
                lbTb = readtable(labelFile, 'Delimiter', ',');
            else
                disp("Skip reading events since combined_reading_labels.csv is missing.");
                return
            end
            
            % Read reading data
            tgFolder = fullfile(this.IGetModulePath("reading"), 'textgrids');
            hasTG = exist(tgFolder, 'dir');
            if hasTG
                % Import TextGrid structs
                [tg, meta] = NP.Preproc.ImportTextgrid(tgFolder);
                
                % Exclude missing items from lbTb
                D = abs(lbTb.OnsetTime - meta.tOn');
                d = min(D, [], 2);
                m = d < 0.003; % use a tolerance of 3ms
                if ~all(m)
                    arrayfun(@(x) warning("Skip the label which onsets at %g for its TextGrid file is missing.", x), lbTb.OnsetTime(~m));
                end
                lbTb = lbTb(m,:);
                
                % Make NP.TGEvent objects
                evt = NP.TGEvent(cat(1, tg{:}));
                evt = evt + meta.tOn; % convert to global time
                evt = Trim(evt); % trim silent ends
                evt = round(evt, 3); % round to millisecond
                evt = evt.SetVfield('task', string(lbTb.Task));
                evt = evt.SetVfield('source', string(lbTb.Source));
                evt = evt.SetVfield('type', string(lbTb.Type));
                evt = evt.SetVfield('labelId', lbTb.(1));
            else
                % No TextGrid data
                tg = [];
                meta = [];
                
                % Make NP.TGEvent objects
                evt = MSessionExplorer.Event(lbTb.OnsetTime);
                evt = evt.SetTfield('tmin', lbTb.OnsetTime);
                evt = evt.SetTfield('tmax', lbTb.OffsetTime);
                evt = evt.SetVfield('task', string(lbTb.Task));
                evt = evt.SetVfield('source', string(lbTb.Source));
                evt = evt.SetVfield('type', string(lbTb.Type));
                evt = evt.SetVfield('labelId', lbTb.(1));
            end
            
            % Get taskTime table
            if ismember('taskTime', this.se.tableNames)
                ttTb = this.se.GetTable('taskTime');
                disp('Update the existing taskTime table.');
            else
                ttTb = table;
                disp('Create a new taskTime table.');
            end
            
            % Add events by Source
            evtNames = string(unique(lbTb.Source));
            for n = evtNames(:)'
                isEvt = strcmp(lbTb.Source, n);
                % Check if the source is already in the table
                if ismember(n, ttTb.Properties.VariableNames)
                    ttTb.(n) = {sort(vertcat(ttTb.(n){1}, evt(isEvt)))};
                else
                    ttTb.(n) = {evt(isEvt)};
                end
                ttTb.(n + "On") = {ttTb.(n){1}.GetTfield('tmin')};
                ttTb.(n + "Off") = {ttTb.(n){1}.GetTfield('tmax')};
            end
            
            % Save taskTime table to SE
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.SetTable('taskTime', ttTb, 'eventTimes', 0);
            this.se.userData.phoneMeta = meta;
            this.se.userData.phoneMeta.tg = tg;
            disp("Added reading events to the 'taskTime' table.");
            
        end

        function AddSpeechEvents(this)
            % Create speech events as NP.TGEvent objects and add them to the 'taskTime' table.
            % Numeric onset and offset times will also be added for compatibility and convenience.
            % Data will be organized in a single epoch with zero as the reference time.
            % 
            %   AddSpeechEvents()
            % 
            
            fprintf('\n');
            
            % Read labels
            labelFolder = this.IGetModulePath("labels");
            labelFile = fullfile(labelFolder, 'combined_speech_labels.csv');
            if exist(labelFile, 'file')
                lbTb = readtable(labelFile, 'Delimiter', ',');
            else
                disp("Skip speech events since combined_speech_labels.csv is missing.");
                return
            end
            
            % Read speech data
            phoneFolder = fullfile(this.IGetModulePath("phones"), 'results', 'all');
            hasTG = exist(phoneFolder, 'dir');
            if hasTG
                % Import TextGrid structs
                [tg, meta] = NP.Preproc.ImportTextgrid(phoneFolder);
                
                % Exclude missing items from lbTb
                D = abs(lbTb.OnsetTime - meta.tOn');
                d = min(D, [], 2);
                m = d < 0.003; % use a tolerance of 3ms
                if ~all(m)
                    arrayfun(@(x) warning("Skip the label which onsets at %g for its TextGrid file is missing.", x), lbTb.OnsetTime(~m));
                end
                lbTb = lbTb(m,:);
                
                % Make NP.TGEvent objects
                evt = NP.TGEvent(cat(1, tg{:}));
                evt = evt + meta.tOn; % convert to global time
                evt = Trim(evt); % trim silent ends
                evt = round(evt, 3); % round to millisecond
                evt = evt.SetVfield('task', string(lbTb.Task));
                evt = evt.SetVfield('source', string(lbTb.Source));
                evt = evt.SetVfield('type', string(lbTb.Type));
                evt = evt.SetVfield('labelId', lbTb.(1));
            else
                % No TextGrid data
                tg = [];
                meta = [];
                
                % Make NP.TGEvent objects
                evt = MSessionExplorer.Event(lbTb.OnsetTime);
                evt = evt.SetTfield('tmin', lbTb.OnsetTime);
                evt = evt.SetTfield('tmax', lbTb.OffsetTime);
                evt = evt.SetVfield('task', string(lbTb.Task));
                evt = evt.SetVfield('source', string(lbTb.Source));
                evt = evt.SetVfield('type', string(lbTb.Type));
                evt = evt.SetVfield('labelId', lbTb.(1));
            end
            
            % Get taskTime table
            if ismember('taskTime', this.se.tableNames)
                ttTb = this.se.GetTable('taskTime');
                disp('Update the existing taskTime table.');
            else
                ttTb = table;
                disp('Create a new taskTime table.');
            end
            
            % Add events by Source
            evtNames = string(unique(lbTb.Source));
            for n = evtNames(:)'
                isEvt = strcmp(lbTb.Source, n);
                % Check if the source is already in the table
                if ismember(n, ttTb.Properties.VariableNames)
                    ttTb.(n) = {sort(vertcat(ttTb.(n){1}, evt(isEvt)))};
                else
                    ttTb.(n) = {evt(isEvt)};
                end
                ttTb.(n + "On") = {ttTb.(n){1}.GetTfield('tmin')};
                ttTb.(n + "Off") = {ttTb.(n){1}.GetTfield('tmax')};
            end
            
            % Save taskTime table to SE
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.SetTable('taskTime', ttTb, 'eventTimes', 0);
            this.se.userData.phoneMeta = meta;
            this.se.userData.phoneMeta.tg = tg;
            disp("Added speech events to the 'taskTime' table.");
        end
        
        function AddIntensity(this, varargin)
            % Add 'inten' table. Data will be organized in a single epoch with zero as the reference time.
            % 
            %   M.AddIntensity()
            % 
            
            fprintf('\n');
            
            % Parse options
            p = inputParser();
            p.addOptional('isForce', this.isOverwrite, @islogical);
            p.parse(varargin{:});
            isForce = p.Results.isForce;
            
            % Check if the data already exist
            tbName = 'inten';
            if ismember(tbName, this.se.tableNames) && ~isForce
                fprintf("Will not overwrite the existing %s data\n", tbName);
                return
            end
            
            % Read acous files
            intensFolder = this.IGetModulePath("intensity");
            if isempty(intensFolder)
                disp("Skip intensity features since the data folder is missing.");
                return
            end
            isCombineEpochs = true;
            [tb, meta] = NP.Preproc.ImportIntensity(intensFolder, isCombineEpochs);
            
            % Save acous table to SE
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.SetTable(tbName, tb, 'timeSeries', 0);
            this.se.userData.intenMeta = meta;
            fprintf("Added '%s' table to se.\n", tbName);
        end
        
        function AddPitch(this, varargin)
            % Add 'pitch' table and columns in taskTime. 
            % Data will be organized in a single epoch with zero as the reference time.
            % 
            %   M.AddPitch()
            % 
            
            fprintf('\n');
            
            % Parse options
            p = inputParser();
            p.addOptional('isForce', this.isOverwrite, @islogical);
            p.parse(varargin{:});
            isForce = p.Results.isForce;
            
            % Check if the data already exist
            tbName = 'pitch';
            if ismember(tbName, this.se.tableNames) && ~isForce
                fprintf("Will not overwrite the existing %s data\n", tbName);
                return
            end
            
            % Read files
            pitchFolder = this.IGetModulePath(tbName);
            if isempty(pitchFolder)
                disp("Skip pitch features since the data folder is missing.");
                return
            end
            isCombineEpochs = true;
            [tsTb, meta] = NP.Preproc.ImportPitch(pitchFolder, isCombineEpochs);
            
            % Save table to SE
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.SetTable(tbName, tsTb, 'timeSeries', 0);
            this.se.userData.pitchMeta = meta;
            fprintf("Added '%s' table to se, and %s events to 'taskTime' table.\n", tbName, tbName);
        end
        
        function AddArtics(this, varargin)
            % Add 'artic' table. Data will be organized in a single epoch with zero as the reference time.
            % 
            %   M.AddArtics()
            % 
            
            fprintf('\n');
            
            % Parse options
            p = inputParser();
            p.addOptional('isForce', this.isOverwrite, @islogical);
            p.parse(varargin{:});
            isForce = p.Results.isForce;
            
            % Check if the data already exist
            tbName = 'artic';
            if ismember(tbName, this.se.tableNames) && ~isForce
                fprintf("Will not overwrite the existing %s data\n", tbName);
                return
            end
            
            % Read artic files
            articFolder = this.IGetModulePath("artics");
            if isempty(articFolder)
                disp("Skip articulatory features since the data folder is missing.");
                return
            end
            subFolders = fullfile(articFolder, ...
                {'F01_indep_Haskins_loss_90_filter_fix_bn_False_0_setting2', ...
                'hprc_train_no_m1f2_h2tv_gru_tanh_nogan'});
            
            k = 0;
            for i = 1 : numel(subFolders)
                if ~exist(subFolders{i}, 'dir')
                    continue
                end
                k = k + 1;
                isCombineEpochs = true;
                [tb, meta] = NP.Preproc.ImportArtics(subFolders{i}, isCombineEpochs);
                
                % Save artic table to SE
                if this.isDry
                    this.INoteIfDryRun();
                    return
                end
                if k == 1
                    tbName = 'artic';
                else
                    tbName = "artic"+k;
                end
                this.se.SetTable(char(tbName), tb, 'timeSeries', 0);
                this.se.userData.articMeta(k) = meta;
                fprintf("Added '%s' table to se.\n", tbName);
            end
        end
        
        function AddLandmark(this, varargin)
            % Add landmark table. Data will be organized in a single epoch with zero as the reference time.
            % 
            %   M.AddLandmark()
            %   M.AddLandmark(isForce)
            %
            % Input
            %   isForce         Whether or not to force overwrite existing data. Default is false.
            % 
            
            fprintf('\n');
            
            % Parse options
            p = inputParser();
            p.addOptional('isForce', this.isOverwrite, @islogical);
            p.parse(varargin{:});
            isForce = p.Results.isForce;
            
            % Check if data already exist
            if ismember('landmark', this.se.tableNames) && ~isForce
                disp('Will not overwrite the existing landmark data');
                return
            end
            
            % Find landmark folder
            lmDir = this.IGetModulePath("landmarks");
            if isempty(lmDir)
                disp("Skip landmark features since the data folder is missing.");
                return
            end
            
            % Get ni table
            if ~ismember('ni', this.se.tableNames)
                warning("Skip landmark features since the 'ni' table is missing.");
                return
            end
            [niTb, tt] = this.se.GetTable('ni', 'taskTime');
            
            % Get the onset and offset times of speech production
            if ~ismember('taskTime', this.se.tableNames)
                warning("Skip landmark features since the 'taskTime' table is missing.");
                return
            end
            if ~all(ismember({'prodOn', 'prodOff'}, tt.Properties.VariableNames))
                warning("Skip landmark features since 'prodOn' and/or 'prodOff' are missing from 'taskTime' table.");
                return
            end
            prodWin = cell2mat([tt.prodOn, tt.prodOff]);
            
            % Import landmark data
            disp("Processing landmarks");
            [tb, meta] = NP.Preproc.ImportLandmarks(lmDir, this.recId, niTb, prodWin);
            
            % Save landmark table to SE
            if this.isDry
                this.INoteIfDryRun();
                return
            end
            this.se.SetTable('landmark', tb, 'timeSeries', 0);
            this.se.userData.landmarkMeta = meta;
            disp("Added 'landmark' table to se.");
        end
        
        function varargout = PosthocFix(this)
            % Recording-specific posthoc fix
            
            fprintf('\n');
            
            switch this.recId
                case 'NP30_B12'
                    % Add texts to reading (text display) events
                    
                    textTb = readtable("C:\chang_lab\project_np\preproc\NP30_B12\labels\mocha\text_timing.txt");
                    tt = this.se.GetTable('taskTime');
                    T = tt.reading{1};
                    T = T.SetVfield('source', repmat("text", size(T)));
                    T = T.SetVfield('type', string(textTb.(3)));
                    tt.reading{1} = T;
                    tt.Properties.VariableNames = strrep(tt.Properties.VariableNames, 'reading', 'text');
                    
                    this.se.SetTable('taskTime', tt);
                    disp("Added texts to reading (text display) events.");
                    
                    varargout{1} = tt;
            end
        end
        
        function [filePaths, imNames] = IFindSglxFiles(this, streamType)
            % Find the most processed SGLX files of ap, lf, or ni stream
            
            streamType = string(lower(streamType));
            
            sglxFolder = this.IGetModulePath("sglx");
            
            if ismember(streamType, ["lf", "ap"])
                % Find file paths from all probes
                fprintf("Searching for the most processed %s file\n", streamType);
                filePattern = this.recId + "_*." + streamType + ".meta"; % e.g. NP00_B0_*.lf.meta
                imSearch = MBrowse.Dir2Table(fullfile(sglxFolder, '*', '*_imec*'));
                imDirNames = unique(string(imSearch.name));
                for i  = numel(imDirNames) : -1 : 1
                    imNames(i) = regexp(imDirNames(i), 'imec\d+$', 'match');
                    fileSearch = MBrowse.Dir2Table(fullfile(sglxFolder, 'mc_*', imDirNames(i), filePattern));
                    if height(fileSearch) == 0
                        fprintf("Motion corrected %s file is not found. Search for CatGT %s file instead.\n", streamType, streamType);
                        fileSearch = MBrowse.Dir2Table(fullfile(sglxFolder, 'catgt_*', imDirNames(i), filePattern));
                    end
                    if height(fileSearch) == 0
                        fprintf("CatGT %s file is not found. Search for original %s file instead.\n", streamType, streamType);
                        fileSearch = MBrowse.Dir2Table(fullfile(sglxFolder, '*', imDirNames(i), filePattern));
                    end
                    if height(fileSearch) == 0
                        filePaths(i) = string(NaN);
                        warning("Cannot find any %s file for %s.", streamType, imNames);
                    else
                        filePaths(i) = string(fullfile(fileSearch.folder{1}, fileSearch.name{1}));
                        fprintf("Found %s\n", filePaths(i));
                    end
                end
                
            elseif streamType == "ni"
                % Find file path from any copy
                disp("Searching for a nidq.bin file");
                filePattern = this.recId + "*.nidq.bin"; % e.g. NP00_B0_*.nidq.bin
                fileSearch = MBrowse.Dir2Table(fullfile(sglxFolder, ['*' this.recId '*'], filePattern));
                if height(fileSearch) == 0
                    filePaths = string(NaN);
                    warning('Cannot find any nidq.bin file.');
                else
                    filePaths = string(fullfile(fileSearch.folder{1}, fileSearch.name{1}));
                    if height(fileSearch) > 1
                        disp('Found more than one nidq.bin files (they should be identical).');
                    end
                    fprintf("Use %s\n", filePaths);
                end
                imNames = string(NaN);
            else
                error("'%s' is not a valid stream type. Must be 'ap', 'lf', or 'ni'", streamType);
            end
            
            % Remove missing paths from outputs
            m = ismissing(filePaths);
            filePaths(m) = [];
            imNames(m) = [];
        end
        
        function p = IGetModulePath(this, moduleName)
            % Build the folder path for a given module based on the convention
            p = fullfile(this.ppRoot, this.recId, moduleName);
            if ~exist(p, 'dir')
                p = '';
            end
        end

        function b = ICheckExpInfo(this, isErr)
            % Check if expInfo is present or not
            b = isfield(this.se.userData, 'expInfo');
            if exist('isErr', 'var') && isErr
                assert(b, 'expInfo must be added before running the extraction');
            end
        end
        
        function INoteIfDryRun(this)
            % Display a message noting this is a dry run
            if this.isDry
                disp('This is a dry run. No data will be added.');
            end
        end
        
        % Not in use
        function [filePaths, probeIds] = FindKilosortResults(this, streamType)
            % Find the most processed SGLX files of ap, lf, or ni stream
            
            fprintf("Searching for most processed KS outputs\n");
            ksFolder = this.IGetModulePath("kilosort");
            ksSearch = MBrowse.Dir2Table(fullfile(ksFolder, "*", "*", "cluster_group.tsv"));
            [~, dirNames] = cellfun(@fileparts, ksSearch.folder, 'Uni', false);
            imNames = regexp(imDirNames(i), 'imec\d+$', 'match');
            
                
        end
    end
end

