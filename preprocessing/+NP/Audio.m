classdef Audio
    
    methods(Static)
        % Waveform
        function FilterSpeech(se)
            % Highpass filter audio signal at 60Hz for speech analysis
            % 
            %   NP.Audio.FilterSpeech(se)
            % 
            
            % Vectorize data
            seLite = se.Duplicate({'ni'}, false);
            seLite.SliceSession(0, 'absolute');
            
            % Prepare parameters for high-pass filtering
            %   The lower cutoff is based on the typical lower bound of male pitch
            %   The following filter has an order of 7
            fs = se.userData.niMeta.niSampRate;
            if ischar(fs)
                fs = str2double(fs);
            end
            D = designfilt('highpassiir', ...
                'PassbandFrequency', 60, ...
                'StopbandFrequency', 50, ...
                'StopbandAttenuation', 70, ...
                'PassbandRipple', 0.1, ...
                'SampleRate', fs, ...
                'DesignMethod', 'ellip');
            
            % Filter audio
            ni = seLite.GetTable('ni');
            if nargin < 2
                colNames = {'mic', 'speaker1'};
            end
            for i = 1 : numel(colNames)
                vn = colNames{i};
                w = ni.(vn){1};
                dataType = class(w);
                w = filtfilt(D, double(w));
                ni.(vn){1} = cast(w, dataType);
            end
            seLite.SetTable('ni', ni);
            
            % Reslice
            tRef = se.GetReferenceTime('ni');
            tSlice = tRef;
            tSlice(1) = 0;
            seLite.SliceSession(tSlice, 'absolute');
            seLite.AlignTime(tRef-tSlice);
            ni = seLite.GetTable('ni');
            se.SetTable('ni', ni);
            se.userData.niMeta.hpFilter = D;
        end
        
        function WriteTrial(se, varargin)
            % Write waveform from trial(s) to individual audio file(s)
            %
            %   NP.Audio.WriteTrial(se)
            %   NP.Audio.WriteTrial(se, trialInd)
            %   NP.Audio.WriteTrial(se, trialInd, tWin)
            %   NP.Audio.WriteTrial(..., 'Channel', "mic")
            %   NP.Audio.WriteTrial(..., 'FolderPath', '')
            %   NP.Audio.WriteTrial(..., 'FileNames', '')
            % 
            
            p = inputParser();
            p.addOptional('trialInd', 1:se.numEpochs, @isnumeric);
            p.addOptional('tWin', [-Inf Inf], @isnumeric);
            p.addParameter('Channel', 'mic', @(x) isstring(x) || ischar(x));
            p.addParameter('FolderPath', '', @(x) isstring(x) || ischar(x));
            p.addParameter('FileNames', '', @(x) isstring(x) || ischar(x) || iscellstr(x));
            p.parse(varargin{:});
            trialInd = p.Results.trialInd;
            tWin = p.Results.tWin;
            chan = string(p.Results.Channel);
            folderPath = p.Results.FolderPath;
            fileNames = cellstr(p.Results.FileNames);
            
            if ~isempty(folderPath) && ~exist(folderPath, 'dir')
                mkdir(folderPath);
            end
            
            ni = se.SliceTimeSeries('ni', tWin, trialInd, chan, 'Fill', 'none');
            tv = se.GetTable('taskValue');
            tv = tv(trialInd,:);
            
            % Check signal range
            w = cat(1, ni.(chan){:});
            wMax = max(abs(w));
            if wMax > 1
                k = 1/wMax/2;
            else
                k = 1;
            end
            
            for i = 1 : height(ni)
                % Get data
                t = ni.time{i};
                w = ni.(chan){i} * k;
                
                % Fill NaN window limits
                tWinVal = tWin;
                if isinf(tWinVal(1))
                    tWinVal(1) = t(1);
                end
                if isinf(tWinVal(2))
                    tWinVal(2) = t(end);
                end
                tWinVal = round(tWinVal, 3);
                
                % Resample waveform to span the window
                t0 = t;
                t = (tWinVal(1) : diff(t(1:2)) : tWinVal(2))';
                w  = interp1(t0, w, t, 'linear', 0);
                
                % Generate audio file
                if i > numel(fileNames) || isempty(fileNames{i})
                    recId = NP.SE.GetID(se);
                    fileName = sprintf("%s_trial%i_%.3f-%.3fs.wav", recId, tv.trialNum(i), tWinVal(1), tWinVal(2));
                else
                    fileName = fileNames{i};
                end
                filePath = fullfile(folderPath, fileName);
                NP.Audio.Write(filePath, w, t);
            end
        end
        
        function Write(filePath, w, timeORfs, toFs)
            % Save audio waveform as a WAV file (default at 44100 Hz)
            %
            %   NP.Audio.Write(filePath, w, timeORfs)
            %   NP.Audio.Write(filePath, w, timeORfs, toFs)
            % 
            if isscalar(timeORfs)
                fromFs = round(timeORfs);
            else
                fromFs = round(1 / diff(timeORfs(1:2)));
            end
            if ~exist('toFs', 'var')
                toFs = 44100;
            end
            w = resample(double(w), toFs, fromFs); % upsample to the target frequency
            audiowrite(filePath, w, toFs);
        end
        
        % Spectrograms
        function AddMelTable(se)
            % Add Mel spectrogram table to se
            % 
            %   NP.Audio.AddMelTable(se)
            %   NP.Audio.AddMelTable(se, colNames)
            % 
            
            if nargin < 2
                colNames = {'mic', 'speaker1'};
            end
            % Vectorize data
            seLite = se.Duplicate({'ni'}, false);
            seLite.SliceSession(0, 'absolute');
            
            % Compute Mel spectrogram
            ni = seLite.GetTable('ni');
            mel = table;
            for i = 1 : numel(colNames)
                vn = colNames{i};
                [S, F, T] = NP.Audio.ComputeMelSpectrogram(ni.(vn){1}, ni.time{1});
                mel.time = {T};
                mel.(vn) = {S'};
            end
            seLite.SetTable('mel', mel, 'timeSeries', seLite.GetReferenceTime);
            seLite.RemoveTable('ni');
            
            % Reslice
            tRef = se.GetReferenceTime('ni');
            tSlice = tRef;
            tSlice(1) = 0;
            seLite.SliceSession(tSlice, 'absolute');
            seLite.AlignTime(tRef-tSlice);
            mel = seLite.GetTable('mel');
            se.SetTable('mel', mel, 'timeSeries', tRef);
            se.userData.melMeta.F = F;
        end
        
        function [S, F, T] = ComputeMelSpectrogram(w, tORfs)
            % Compute Mel spectrogram with preset parameters.
            % If waveform timestamps are provided, the spectrogram timestamps will be kept consistent.
            % 
            %   [S, F, T] = NP.Audio.ComputeMelSpectrogram(w, tORfs)
            % 
            % Inputs
            %   w           A vector of audio waveform.
            %   tORfs       1) A vector of sample timestamps that matches the size of w. In this case,
            %                  the spectrogram time T will be kept consistent with this waveform time.
            %               2) A scalar specifying the sampling frequency.
            % Outputs
            %   S           A m-by-n matrix of magnitude values. m is the number of frequency bins, 
            %               n is the number of time bins.
            %   F           A m-element vector of frequencies in Herz.
            %   T           A n-element vector of timestamps in second.
            % 
            % See also melSpectrogram, hamming
            
            % Compute Mel spectrogram
            if isscalar(tORfs)
                fs = tORfs;
            else
                fs = 1 / diff(tORfs(1:2));
            end
            win = hamming(round(0.025 * fs), 'periodic');
            overlap = round(length(win) / 2);
            [S, F, T] = melSpectrogram(w, fs, 'NumBands', 80, 'Window', win, 'OverlapLength', overlap, 'FrequencyRange', [0 8e3], 'FFTLength', 2048);
            
            % Scale power to dB
            S =  10*log10(S + cast(eps,'like',S));
            
            % Convert relative spectrogram time to absolute time
            if ~isscalar(tORfs)
                T = double(T); % avoid out of range
                T = T + tORfs(round(length(win)/2));
            end
        end
        
        function PlotMelSpectrogram(varargin)
            % Plot Mel strectrogram with 
            % 
            %   NP.Audio.PlotMelSpectrogram(ax, S, F, T)
            %   NP.Audio.PlotMelSpectrogram(ax, S, F, T, dt)
            %   NP.Audio.PlotMelSpectrogram(ax, S, F, T, dt, yOffset)
            % 
            % Inputs
            %   ax          Axes to plot the spectrogram in.
            %   S           A m-by-n matrix of magnitude values. m is the number of frequency bins, 
            %               n is the number of time bins.
            %   F           A m-element vector of frequencies in Herz.
            %   T           A n-element vector of timestamps in second.
            %   dt          A scalar in sec to interpolate S along the time axis. Default is 0.005s.
            %   yOffset     Offset of the spectgrogram in the y-axis (in Mels).
            % 
            % See also NP.Audio.ComputeMelSpectrogram
            
            NP.Audio.PlotSpectrogram(varargin{:}, 'Scale', 'mel');
        end
        
        function PlotSpectrogram(ax, S, F, T, varargin)
            % Plot strectrogram
            % 
            %   NP.Audio.PlotSpectrogram(ax, S, F, T)
            %   NP.Audio.PlotSpectrogram(ax, S, F, T, dt)
            %   NP.Audio.PlotSpectrogram(ax, S, F, T, dt, yOffset)
            %   NP.Audio.PlotSpectrogram(..., 'Scale', 'linear')
            % 
            % Inputs
            %   ax          Axes to plot the spectrogram in.
            %   S           A m-by-n matrix of magnitude values. m is the number of frequency bins, 
            %               n is the number of time bins.
            %   F           A m-element vector of frequencies in Herz.
            %   T           A n-element vector of timestamps in second.
            %   dt          A scalar in sec to interpolate S along the time axis. Default is 0.005s.
            %   yOffset     Offset of the spectgrogram in the y-axis (in Mels).
            %   'Scale'     The scale the spectrogram is plotted in, 'mel' or 'linear'. Note that F 
            %               is always in Hz regardless of the chosen scale.
            % 
            % See also NP.Audio.ComputeMelSpectrogram
            
            p = inputParser;
            p.addOptional('dt', [], @(x) isnumeric(x));
            p.addOptional('yOffset', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Scale', 'linear', @(x) ischar(x) || isstring(x));
            p.addParameter('FTick', [], @isnumeric);
            p.parse(varargin{:});
            dt = p.Results.dt;
            yOffset = p.Results.yOffset;
            scaleName = lower(p.Results.Scale);
            tickVal = p.Results.FTick;
            
            if isempty(dt)
                dt = 0.005; % sec
            end
            
            % Interpolate time
            T = single(T);
            Tq = (T(1) : dt : T(end))';
            if scaleName == "mel"
                F = hz2mel(F);
            end
            Fq = linspace(F(1), F(end), numel(F)*4);
            S = interp2(F, T, S, Fq, Tq, 'linear');
            
            % Get tick frequencies
            Fq = Fq + yOffset;
            if isempty(tickVal)
                nFq = numel(Fq);
                nTick = 5;
                tickInd = 1 : ceil(nFq/nTick) : nFq;
                tickVal = Fq(tickInd);
            end
            
            % Plot
            imagesc(ax, Tq, Fq, -S');
            
            ax.XLim = [Tq(1)-dt/2, Tq(end)+dt/2];
            ax.YDir = 'normal';
            ax.YTick = tickVal;
            if scaleName == "mel"
                tickVal = mel2hz(tickVal - yOffset);
            else
                tickVal = tickVal - yOffset;
            end
            if max(tickVal) > 1e3
                ax.YTickLabel = round(tickVal/1e3, 1); % convert to kHz
                ax.YLabel.String = 'Frequency (kHz)';
            else
                ax.YTickLabel = round(tickVal); % convert to kHz
                ax.YLabel.String = 'Frequency (Hz)';
            end
            colormap(ax, "gray");
            cLims = prctile(-S(:), [1 99]);
            ax.CLim = cLims + [-10 -15];
            % ax.CLim(2) = max(ax.CLim(1)+40, cLims(2)-10);
            % ax.CLim(2) = 70 * sign(ax.CLim(2)); % sign deals with inverted color
            MPlot.Axes(ax);
        end
        
        % Features
        function AddEnvelopRate(se)
            % Add the time derivative of envelop, dEnv, to 'inten' table
            % 
            %   AddEnvelopRate(se)
            % 
            inten = se.GetTable('inten');
            inten.dEnv = cellfun(@(y,t) max(gradient(y)./gradient(t), 0), inten.env, inten.time, 'Uni', false);
            se.SetTable('inten', inten);
        end
        
    end
    
end
