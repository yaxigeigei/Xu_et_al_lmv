classdef Preproc
    % Import preprocessing outputs to data tables and metadata structs
    
    methods(Static)
        % Importers
        function [lfTb, meta] = ImportLFP(lfPath, chanMapPath)
            % Read, downsample, and condition SpikeGLX LFP data
            % 
            %   [lfTb, meta] = ImportLFP(lfPath, chanMapPath)
            % 
            
            % Check files
            [folder, baseName, ext] = fileparts(lfPath);
            if strcmp(ext, 'lf')
                baseName = lfPath;
            end
            metaPath = fullfile(folder, [baseName '.meta']);
            binPath = fullfile(folder, [baseName '.bin']);
            assert(exist(metaPath, 'file'), 'The lf.meta file cannot be found.');
            assert(exist(binPath, 'file'), 'The lf.bin file cannot be found.');
            
            % Load LFP data
            disp('Reading lf.meta and lf.bin ...');
            [meta, v, t] = MSpikeGLX.ReadLFP(metaPath);
            
            % Remove the reference channel
            v(:,193) = [];
            
            % Sort channels by the order in chanTb (i.e. by depth)
            [~, chanTbKS] = MKilosort2.LoadChanMap2Table(chanMapPath);
            v = v(:, chanTbKS.sortInd);
            
            % Downsample from 2500 to 400 Hz
            fromFs = str2double(meta.imSampRate);
            toFs = fromFs * 400/2500; % this handles non-integer fs
            fprintf('Downsampling LFP from %g to %gHz ...\n', fromFs, toFs);
            v = resample(double(v), round(toFs), round(fromFs));
            t = ((1:size(v,1))-0.5) / toFs;
            
            lfTb = table();
            lfTb.time{1} = t';
            lfTb.v{1} = single(v);
        end
        
        function ss = ImportKilosortResults(ksFolders)
            % Import Kilosort results in a struct array
            %
            %   ss = NP.Preproc.ImportKilosortResults(ksFolders)
            % 
            
            ksFolders = cellstr(ksFolders);
            ss = cell(size(ksFolders));
            
            for i = 1 : numel(ksFolders)
                % Find or make SortingResults object
                srSearch = MBrowse.Dir2Table(fullfile(ksFolders{i}, '**', 'sr_*.mat'));
                if ~isempty(srSearch)
                    % Check available sr files
                    srPath = fullfile(srSearch.folder{end}, srSearch.name{end});
                    if height(srSearch) > 1
                        warning('%i sr files are found. Use the last file %s', height(srSearch), srPath);
                    end

                    % Load existing object
                    fprintf('Reading sr file\n  %s\n', srPath);
                    load(srPath, 'sr');
                    
                    % Recompute all metrics
                    sr.ComputeAll();
                else
                    % Import Kilosort/Phy output
                    fprintf('Reading Kilosort/Phy output\n');
                    sr = MTracer.KilosortResult();
                    sr.ImportData(ksFolders{i});
                    sr.ComputeAll();
                    
                    % Cache sr object
                    srPath = fullfile(ksFolders{i}, 'sr.mat');
                    save(srPath, "sr");
                end

                % Store varibales in a struct
                s.srFile = srPath;
                s.chanMapFile = sr.chanMapFile;
                s.chanMapName = sr.chanMapName;
                s.chanMapTb = sr.chanTb;
                s.spkTb = sr.spkTb;
                s.clusTb = sr.clusTb;
                s.tempTb = sr.tempTb;
                ss{i} = s;
            end
            
            ss = [ss{:}];
        end
        
        function [stTb, ksData] = ImportSpikes(ksFolders)
            % Import Kilosort results to a spikeTime table and a metadata struct
            % 
            %   [stTb, ksData] = NP.Preproc.ImportSpikes(ksFolders)
            % 
            
            ksData = NP.Preproc.ImportKilosortResults(ksFolders);
            
            % Remove some bulky data
            col2rm = {'pcWeights', 'pcFeatures'};
            for i = 1 : numel(col2rm)
                if ismember(col2rm{i}, ksData.spkTb.Properties.VariableNames)
                    ksData.spkTb.pcWeights = [];
                end
            end
            col2rm = {'pcChanInd'};
            for i = 1 : numel(col2rm)
                if ismember(col2rm{i}, ksData.tempTb.Properties.VariableNames)
                    ksData.tempTb.pcWeights = [];
                end
            end
            
            % Keep only SU and MU data
            clusTb = ksData.clusTb;
            clusTb = sortrows(clusTb, 'depth', 'descend');
            clusTb = clusTb(~strcmp(clusTb.group, 'noise'), :);
            ksData.clusTb = clusTb;
            
            spkTb = ksData.spkTb;
            spkTb = spkTb(ismember(spkTb.clusId, clusTb.clusId), :);
            ksData.spkTb = spkTb;
            
            % Make spike time table
            if isnan(spkTb.timeSecAdj)
                st = spkTb.timeSec;
            else
                st = spkTb.timeSecAdj;
            end
            
            uid = clusTb.clusId;
            for i = numel(uid) : -1 : 1
                ust{i} = st(spkTb.clusId == uid(i));
            end
            
            stTb = MSessionExplorer.MakeEventTimesTable(ust, 'VariableNames', "u" + uid, 'Uni', false);
        end
        
        function [tg, meta] = ImportTextgrid(phoneDir)
            % Read TextGrid files to mPraat TextGrid structs
            % 
            %   [tg, meta] = NP.Preproc.ImportTextgrid(phoneDir)
            % 
            
            % Find data files with the correct name pattern
            fileSearch = NP.Preproc.GetFeatFiles(fullfile(phoneDir, '*.TextGrid'));
            
            % Store file info
            meta.filename = fileSearch.name;
            meta.folder = phoneDir;
            meta.tOn = fileSearch.tOn;
            meta.tOff = fileSearch.tOff;
            
            % Read text files
            for i = height(fileSearch) : -1 : 1
                tgPath = fullfile(fileSearch.folder{i}, fileSearch.name{i});
                try
                    tg{i,1} = tgRead(tgPath);
                catch
                    warning('Failed to read %s', fileSearch.name{i});
                end
            end
        end
        
        function [tb, meta] = ImportIntensity(folderPath, isCombineEpochs)
            % Read intensity features from NPY files to a table
            % 
            %   [tb, meta] = NP.Preproc.ImportIntensity(folderPath, isCombineEpochs)
            % 
            
            if ~exist('isCombineEpochs', 'var')
                isCombineEpochs = false;
            end
            
            % Find data files with the correct name pattern
            fileTb = NP.Preproc.GetFeatFiles(fullfile(folderPath, '*_*_*.npy'));
            meta.file_name = fileTb.name;
            meta.folder = folderPath;
            
            % Read npy files
            v = cellfun(@readNPY, fullfile(fileTb.folder, fileTb.name), 'Uni', false);
            
            % Make timeSeries table
            varNames = {'env', 'peakEnv', 'peakRate'};
            tb = NP.Preproc.MakeTimeseriesTable(fileTb.tOn, fileTb.tOff, v, isCombineEpochs, varNames);
        end
        
        function [tb, meta] = ImportPitch(folderPath, isCombineEpochs)
            % Read pitch features from NPY files to a table
            % 
            %   [tb, meta] = NP.Preproc.ImportPitch(folderPath, isCombineEpochs)
            % 
            
            if ~exist('isCombineEpochs', 'var')
                isCombineEpochs = false;
            end
            
            % Find data files with the correct name pattern
            fileTb = NP.Preproc.GetFeatFiles(fullfile(folderPath, '*_*_*.npy'));
            meta.file_name = fileTb.name;
            meta.folder = folderPath;
            
            % Read npy files
            v = cellfun(@readNPY, fullfile(fileTb.folder, fileTb.name), 'Uni', false);
            
            % Make timeSeries table
            varNames = {'F0raw', 'F0'};
            tb = NP.Preproc.MakeTimeseriesTable(fileTb.tOn, fileTb.tOff, v, isCombineEpochs, varNames);
        end
        
        function [tb, meta] = ImportArtics(folderPath, isCombineEpochs)
            % Read inversed articulatory features from NPY files to a table
            % 
            %   [tb, meta] = NP.Preproc.ImportArtics(folderPath, isCombineEpochs)
            % 
            
            if ~exist('isCombineEpochs', 'var')
                isCombineEpochs = false;
            end
            
            % Find data files with the correct name pattern
            fileTb = NP.Preproc.GetFeatFiles(fullfile(folderPath, '*_*_*.npy'));
            meta.file_name = fileTb.name;
            meta.folder = folderPath;
            
            % Read npy files
            v = cellfun(@readNPY, fullfile(fileTb.folder, fileTb.name), 'Uni', false);
            
            % Make timeSeries table
            switch size(v{1},2)
                case 18
                    % For F01_indep_Haskins_loss_90_filter_fix_bn_False_0_setting2
                    varNames = {'tt_x' 'tt_y' 'td_x' 'td_y' 'tb_x' 'tb_y' 'li_x' 'li_y' 'ul_x' 'ul_y' 'll_x' 'll_y' 'la' 'pro' 'ttcl' 'tbcl' 'v_x' 'v_y'};
                case 9
                    % For hprc_train_no_m1f2_h2tv_gru_tanh_nogan
                    varNames = {'la', 'lp', 'ja', 'ttcl', 'ttcd', 'tmcl', 'tmcd', 'trcl', 'trcd'};
                otherwise
                    varNames = "akt" + (1:size(v{1},2));
            end
            tb = NP.Preproc.MakeTimeseriesTable(fileTb.tOn, fileTb.tOff, v, isCombineEpochs, varNames);
        end
        
        function [tb, meta] = ImportLandmarks(lmDir, recId, niTb, tRanges)
            % TBW
            %
            %   [tb, meta] = NP.Preproc.ImportLandmarks(lmDir, recId, niTb, tRanges)
            % 
            
            % Parameters
            dsFs = 4e3; % downsampling Fs for faster and better xcorr result
                        % results are stable at least from 2e3 to 8e3 Hz
            
            
            % Load camera audio
            [aCam, fsCam] = audioread(fullfile(lmDir, [recId '_audio.wav']));
            aCam = aCam(:,1); % taking one channel of the stereo
            
            % Downsample video audio
            rDs = round(fsCam / dsFs);
            aCam = decimate(aCam, rDs);
            fsCam = fsCam / rDs;
            tCam = (0:numel(aCam)-1)' / fsCam;
            
            
            % Get NI mic audio
            aNI = niTb.mic{1};
            tNI = niTb.time{1};
            
            % Downsample NI audio
            fsNI = 1 / diff(tNI(1:2));
            rDs = round(fsNI / dsFs);
            aNI = decimate(double(aNI), rDs);
            tNI = downsample(tNI, rDs);
            fsNI = fsNI / rDs;
            
            
            % Find time offset epoch by epoch
            itvl = diff(tRanges, 1, 2);
            tRanges(itvl < 1, :) = []; % ignore very short epochs
            t0NI = mean(tRanges, 2);
            t0Cam = zeros(size(tRanges, 1), 1);
            for i = 1 : size(tRanges, 1)
                % Find samples in this epoch
                tInd = tNI >= tRanges(i,1) & tNI < tRanges(i,2);
                
                % Resample camera audio to match NI audio
                aCamRs = interp1(tCam, aCam, tNI(tInd), 'linear', 'extrap');
                
                % Compute crosscorrelation
                maxLag = round(fsNI * 2); % allows up to a 2-second offset
                [r, lags] = xcorr(aCamRs, aNI(tInd), maxLag);
                
                % Find time offset
                [~, I] = max(r);
                t0Cam(i) = t0NI(i) + lags(I) / fsNI;
                
%                 MPlot.Figure(i); clf
%                 subplot(3,1,1)
%                 plot(tNI(tInd), zscore(aNI(tInd)))
%                 subplot(3,1,2)
%                 plot(tNI(tInd), zscore(aCamRs))
%                 subplot(3,1,3)
%                 plot(lags, r)
            end
            
            % Fit a curve between epoch time offsets and video epoch start times
            dt = t0NI - t0Cam;
            
            MPlot.Figure(889); clf
            plot(t0Cam, dt, 'o'); hold on
            
            if recId == "NP30_B12"
                dt = medfilt1(dt, 3);
                plot(t0Cam, dt, 'o');
                ft = fittype('a + b*normcdf(x,mu,sig) + c*x', 'independent', 'x');
                mdl = fit(t0Cam, dt, ft, 'start', [-0.007, -0.9, 0, 465, 1.3]);
            else
                mdl = fit(t0Cam, dt, 'poly1', 'Robust', 'Bisquare');
            end
            plot(mdl);
            
            
            % Load tracking results
            s = load(fullfile(lmDir, [recId '.mat']));
            coords = s.coords;
            tFr0 = (0.5:size(coords,3))' / s.frame_rate; % use mid-frame times
            
            % Remove empty frames
            isEmp = isnan(squeeze(coords(1,1,:)));
            coords(:,:,isEmp) = [];
            tFr0(isEmp) = [];
            
            % Apply time correction
            tFrX = tFr0 + mdl(tFr0);
            
            % Extract coordinates by the groups of interest
            sXY = NP.LM.FindXY(coords);
            
            % Make landmark table
            colData = [tFr0 struct2cell(sXY)'];
            colNames = [{'camTime'} fieldnames(sXY)'];
            tb = MSessionExplorer.MakeTimeSeriesTable(tFrX, colData, 'VariableNames', colNames);
            
            % Add derived quantities
            [tb.mouthHeight, tb.mouthWidth] = arrayfun(@NP.LM.GetMouthSize, sXY, 'Uni', false);
            tb.chin2nose = arrayfun(@NP.LM.GetChin2NoseDist, sXY, 'Uni', false);
            [tb.leftPupilX, tb.leftPupilY, tb.rightPupilX, tb.rightPupilY] = arrayfun(@NP.LM.GetPupilPositions, sXY, 'Uni', false);
            
            % Edit metadata
            s.audio_fs = fsCam;
            meta = rmfield(s, 'coords');
            meta.dtMdl = mdl;
        end
        
        function tb = ImportMetaSheet(id, tabName)
            % Read the master metadata spreadsheet.
            % It's downloaded from https://docs.google.com/spreadsheets/d/1cpAUkyLmi17d3mwPXLCd5q8Tk0Q-RYYPbwXQah3A8OE
            % and saved using the name in NP.Data.metaSheetName
            % 
            %   tb = NP.Preproc.ImportMetaSheet()
            %   tb = NP.Preproc.ImportMetaSheet(id, tabName)
            %   tb = NP.Preproc.ImportMetaSheet([], tabName)
            % 
            
            if nargin < 2
                tabName = "Recordings";
            end
            warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
            tb = readtable(NP.Data.metaSheetName, 'Sheet', tabName);
            warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');
            
            if ismember(tabName, ["Recordings", "Archived Recordings"])
                allIds = cellfun(@(x,y) x + "_" + y, tb.Subject, tb.Block);
            elseif tabName == "Subjects"
                allIds = tb.NP;
            else
                error("Importing metadata from '%s' tab is not supported.", tabName);
            end
            
            if nargin < 1 || isempty(id)
                id = allIds;
            end
            id = string(id);
            
            indHit = [];
            for i = 1 : numel(id)
                isHit = strcmp(id(i), allIds);
                if ~any(isHit)
                    warning("No entry matches with '%s' in the '%s' tab.", id(i), tabName);
                    continue
                elseif sum(isHit) > 1
                    error("More than one entries match with '%s' in the '%s' tab.", id(i), tabName);
                else
%                     fprintf("Found %s metadata for %s\n", tabName, id(i));
                end
                isHit = find(isHit, 1);
                indHit(end+1) = isHit;
            end
            tb = tb(indHit,:);
            
            tb = table2struct(tb);
            tb = MUtil.EvalFields(tb);
            tb = struct2table(tb, 'AsArray', true);
        end
        
        % Utilties
        function fileSearch = GetFeatFiles(filePattern)
            % Search for feature files, parse filenames for onset and offset times, and return
            % file paths in the order of onset times
            %
            %   fileSearch = GetFeatFiles(filePattern)
            % 
            % Input
            %   filePattern     A path pattern with wildcards to match the paths of feature files.
            % Outputs
            %   fileSearch      A table that contains basic file info and the following columns.
            %                   tOn         A vector of onset times.
            %                   tOff        A vector of offset times.
            %                   Rows are sorted by increasing tOn.
            %
            
            % Find data files with the correct name pattern
            fileSearch = MBrowse.Dir2Table(filePattern);
            
            % Extract onset and offset time from file name
            [~, bareNames] = cellfun(@fileparts, fileSearch.name, 'Uni', false); % remove extension
            nameParts = cellfun(@(x) strsplit(x, '_'), bareNames, 'Uni', false);
            tOn = cellfun(@(x) str2double(x(2)), nameParts);
            tOff = cellfun(@(x) str2double(x(3)), nameParts);
            
            % Sort with increasing onset time
            fileSearch.bareName = bareNames;
            fileSearch.tOn = tOn;
            fileSearch.tOff = tOff;
            fileSearch = sortrows(fileSearch, 'tOn');
        end
        
        function tb = MakeTimeseriesTable(tOn, tOff, v, isCombineEpochs, varNames)
            % Make timeSeries table from data
            % 
            %   tb = MakeTimeseriesTable(tOn, tOff, v, varNames, isCombineEpochs)
            % 
            
            % Make timestamps
            nEpSp = cellfun(@(x) size(x,1), v);
            t = arrayfun(@(a,b,n) linspace(a,b,n)', tOn, tOff, nEpSp, 'Uni', false);
            
            % Re-partition cell array
            v = cell2mat(v);
            [nSp, nVar] = size(v);
            if isCombineEpochs
                t = {cell2mat(t)};
                v = mat2cell(v, nSp, ones(nVar,1));
            else
                v = mat2cell(v, nEpSp, ones(nVar,1));
            end
            
            % Construct a timeseries table
            tb = table;
            tb.time = t;
            for i = 1 : numel(varNames)
                tb.(varNames{i}) = v(:,i);
            end
        end
        
        function tg = TIMIT2TextGrid(timitDir, id)
            % Construct textgrid compatible struct from wrd, phn
            
            id = cellstr(id);
            tg = cell(size(id));
            
            for i = 1 : numel(id)
                % Read wrd and phn files as tables
                wrdFile = fullfile(timitDir, id{i} + ".wrd");
                phnFile = fullfile(timitDir, id{i} + ".phn");
                if ~exist(wrdFile, 'file') || ~exist(wrdFile, 'file')
                    warning('Cannot find wrd and/or phn files with stim ID: ''%s''', id{i});
                    continue
                end
                wrdTb = readtable(wrdFile, 'FileType', 'text');
                phnTb = readtable(phnFile, 'FileType', 'text');
                
                % Remove silent phones
                phnTb(strcmp(phnTb.(3), 'h#'), :) = [];
                phnTb.(3) = upper(phnTb.(3));
                
                % Hardcode the audio sampling rate of TIMIT
                fs = 16000;
                
                % Construct TextGrid compatible struct
                w = struct;
                w.name = 'words';
                w.type = 'interval';
                w.T1 = wrdTb.(1) / fs;
                w.T2 = wrdTb.(2) / fs;
                w.Label = wrdTb.(3)';
                
                p = struct;
                p.name = 'phones';
                p.type = 'interval';
                p.T1 = phnTb.(1) / fs;
                p.T2 = phnTb.(2) / fs;
                p.Label = phnTb.(3)';
                
                s = struct;
                s.tier = {w, p};
                s.tmin = w.T1(1);
                s.tmax = w.T2(end);
                
                tg{i} = s;
            end
            
            % Remove empty elements
            isEmp = cellfun(@isempty, tg);
            tg(isEmp) = [];
        end
        
        function result = ComputeDelimiter(sig, fs, varargin)
            % Process signal used for delimiting
            
            % Handle user inputs
            p = inputParser();
            p.addParameter('ValueFunc', @(x) x, @(x) isa(x, 'function_handle'));
            p.parse(varargin{:});
            valFunc = p.Results.ValueFunc;
            
            % Generate timestamps based on sampling rate
            t = (0:numel(sig)-1)' / fs;
            
            % Find rising edges in delimiter signal
            indEdge = MMath.Logical2Bounds(sig); % find pairs of rising and falling edges
            indRise = indEdge(:,1);
            indFall = indEdge(:,2) + 1;
            indFall(end) = min(indFall(end), numel(sig));
            tRise = t(indRise);
            tFall = t(indFall);
            
            % Calculate delimiter values
            dur = tFall - tRise;
            val = valFunc(dur);
            
            % Output
            result.delimiterRiseTime = tRise;
            result.delimiterDur = dur;
            result.delimiterValue = val;
        end
        
    end
    
end
