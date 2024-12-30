classdef Pitch
    
    methods(Static)
        function [f0, T, F, S] = ComputePitch(w, tORfs, method)
            % Estimate pitch or F0 from audio waveform
            % 
            %   [f0, T, F, S] = NP.Pitch.ComputePitch(w, t)
            %   [f0, T, F, S] = NP.Pitch.ComputePitch(w, fs)
            %   [f0, T, F, S] = NP.Pitch.ComputePitch(..., method)
            % 
            % Inputs
            %   w           A vector of audio waveform.
            %   t           A vector of sample timestamps that matches the size of w. In this case,
            %               the spectrogram time T will be kept consistent with this waveform time.
            %   fs          A numeric scalar specifying the sampling frequency in Hz.
            %   method      1) Any method supported by MTALAB's built-in pitch function. 
            %               2) 'pitchnn'. This calls MATLAB's pitchnn function, which requires prior 
            %                             installation of the Crepe network. 
            %               3) 'STRAIGHT'. This calls the func_StraightPitch function from VoiceSauce 
            %                              package. The package folder path should be added. 
            %               Default is 'NCF' from the built-in pitch function. Method names are not
            %               case-sensitive.
            % Outputs
            %   f0          A n-element vector of pitch contour in Hz. n is the number of time bins.
            %   T           A n-element vector of timestamps in second.
            %   F           A m-element vector of frequencies in Hz. m is the number of frequency bins. 
            %   S           A m-by-n matrix of salience.
            % 
            % See also pitch, pitchnn, func_StraightPitch
            
            if isscalar(tORfs)
                fs = tORfs;
            else
                fs = 1 / diff(tORfs(1:2));
            end
            
            if nargin < 3
                method = 'NCF';
            end
            
            F = [];
            S = [];
            
            if strcmpi(method, 'pitchnn')
                % Extract pitch contour and salience map
                [f0, T, S] = pitchnn(w, round(fs), 'ModelCapacity', 'full', 'ConfidenceThreshold', 0.5);
                
                % Get frequency scale using a diagnal salience map
                F = crepePostprocess(diag(ones(360,1)));
                
            elseif strcmpi(method, 'STRAIGHT')
                % Extract pitch contour
                ops.maxstrdur = 100;    % batch size in sec
                ops.minstrF0 = 60;      % lower F0 bound in Hz
                ops.maxstrF0 = 300;     % higher F0 bound in Hz
                [f0, ~] = func_StraightPitch(double(w), fs, ops);
                T = (1 : numel(f0))' * 1e-3; % f0 is by default estimated at 1kHz
                
            else
                % Extract pitch contour
                [f0, loc] = pitch(w, round(fs), 'Method', method, 'Range', [50 250]);
                T = (loc-1) / fs;
            end
            
            % Convert relative spectrogram time to absolute time
            if ~isscalar(tORfs)
                T = double(T); % avoid out of range
                T = T + tORfs(1);
            end
        end
        
        function f0 = RefinePitch(f0, t, sentObj)
            % Setting voiceless periods to NaN followed by piecewise Svitzky-Golay smoothing
            % 
            %   f0 = NP.Pitch.RefinePitch(f0, t, sentObj)
            % 
            
            % Set voiceless periods to NaN
            phn = Cut(Cut(sentObj));
            m = phn(phn.IsVoiced).MaskTimestamps(t);
            f0(~m) = NaN;
            
            % Remove outliers piecewise
            bb = MMath.Logical2Bounds(m);
            len = diff(bb, 1, 2) + 1;
            dt = diff(t(1:2));
            for i = 1 : size(bb, 1)
                % Detect and fill outliers
                ind = bb(i,1) : bb(i,2);
                f0(ind) = filloutliers(f0(ind), "nearest");
                
                % Smoothing
                order = 2;
                framelen = round(0.1/dt/2)*2+1; % 0.1 sec; sgolayfilt requires even frame length
                if len(i) >= framelen
                    f0(ind) = sgolayfilt(f0(ind), order, framelen);
                end
                
                % Set boundary samples to NaN for easier interpolation and plotting
                f0(ind([1 end])) = NaN;
            end
        end
        
        function EnrichPitchTable(se, speechEvents)
            % Post-process and derive features in the pitch of an se
            % 
            %   NP.Pitch.EnrichPitchTable(se)
            %   NP.Pitch.EnrichPitchTable(se, speechEvents)
            % 
            
            if nargin < 2
                speechEvents = {'stim', 'prod'};
            end
            
            % Add trial events for sentence masking
            if ~ismember('trial', se.tableNames)
                LMV.SE.AddTrialEvents(se);
            end
            
            % Vectorize se
            tbNames = {'taskTime', 'trial', 'pitch'};
            se1 = se.Duplicate(tbNames, false);
            se1.SliceSession(0, 'absolute');
            se1.Column2Cell('trial');
            [tt, tr, pitch] = se1.GetTable(tbNames{:});
            
            % Compute relative pitch
            t = pitch.time{1};
            F0 = pitch.F0{1};
            lnF0 = log(F0);
            rLnF0 = NaN(size(lnF0));
            for i = 1 : numel(speechEvents)
                tg = tt.(speechEvents{i}){1};
                tgMask = tg.MaskTimestamps(t);
                if speechEvents{i} == "stim"
                    % Normalize ln pitch by each sentence during stim
                    for j = 1 : width(tr)
                        te = tr.(j){1};
                        teMask = te.MaskTimestamps(t);
                        mm = tgMask & teMask;
                        lnF0(mm) = filloutliers(lnF0(mm), NaN);
                        rLnF0(mm) = MMath.Normalize(lnF0(mm), 'zscore');
                    end
                elseif speechEvents{i} == "prod"
                    % Normalize ln pitch during prod
                    lnF0(tgMask) = filloutliers(lnF0(tgMask), NaN);
                    rLnF0(tgMask) = MMath.Normalize(lnF0(tgMask), 'zscore');
                else
                    error("'%s' is not a supported speechEvent", speechEvents{i});
                end
            end
            rF0 = exp(rLnF0);
            pitch.rF0{1} = rF0;
            
            % Binning
            pitch.brF0{1} = NP.Pitch.BinPitch(pitch.rF0{1});
            
            % Compute the rate of change
            dt = gradient(t);
            pitch.dF0{1} = gradient(pitch.F0{1}) ./ dt;
            pitch.drF0{1} = gradient(pitch.rF0{1}) ./ dt;
            
            % Fujisaki decomposition
            phr = NaN(size(F0));
            acc = NaN(size(F0));
            for i = 1 : numel(speechEvents)
                tg = tt.(speechEvents{i}){1};
                for j = 1 : numel(tg)
                    m = tg(j).MaskTimestamps(t);
                    [phr(m), acc(m)] = NP.Pitch.FujisakiDecomp(rF0(m), t);
                end
            end
            pitch.voicing{1} = ~isnan(F0);
            pitch.phrase{1} = phr;
            pitch.accent{1} = acc;
            
            % Slice back
            se1.SetTable('pitch', pitch);
            se1.RemoveTable('taskTime');
            tRef = se.GetReferenceTime('pitch');
            tSlice = tRef;
            tSlice(1) = 0;
            se1.SliceSession(tSlice, 'absolute');
            se1.AlignTime(tRef-tSlice);
            pitch = se1.GetTable('pitch');
            se.SetTable('pitch', pitch);
        end
        
        function [brF0, binEdges] = BinPitch(rF0)
            % Bin relative pitch using the method from Tang et al. 2017 with two modifications to overcome outliers.
            % 1. Cap the upper limit to 6, since most of the real rF0 variations are within 4 based on NP41_B1 (male) and NP44_B2 (female).
            % 2. Exclude values outside the upper limit rather than grouping them into the last bin.
            % 
            %   [brF0, binEdges] = NP.Pitch.BinPitch(rF0)
            % 
            
            % Determine binEdges
            lims = prctile(rF0, [2.5 97.5]);
            lims(2) = min(lims(2), 6);
            nBins = 10;
            binEdges = linspace(lims(1), lims(2), nBins);
            binEdges(1) = -Inf; 
            
            % Find bin memberships
            [~, ~, iBin] = histcounts(rF0, binEdges);
            brF0 = NaN(numel(rF0), nBins);
            for i = 1 : nBins
                brF0(iBin==i, i) = 1;
            end
        end
        
        function [phr, acc, fb] = FujisakiDecomp(f0, t)
            % Decompose pitch contour (f0) into phrase (p), accent (ac), and baseline (fb)
            % 
            %   [phr, acc, fb] = FujisakiDecomp(f0, t)
            % 
            % Inputs
            %   f0          A vector of ln(pitch). Voiceless periods should be in NaN.
            %   t           A vector of timestamps with the same size as f0.
            % Outputs
            %   phr         Phrase.
            %   acc         Accent.
            %   fb          Baseline - this is not implemented and will return a vector of zeros.
            % 
            
            % Bridge the gaps in f0
            f0 = fillmissing(f0, 'linear', 'EndValues', 'none');
            
            % Get a proxy for phrase by highpass filtering f0
            acc = NaN(size(f0));
            m = ~isnan(f0);
            fs = 1/diff(t(1:2));
            acc(m) = highpass(f0(m), 1.5, fs);
            
            % Take the remaining compoent in f0 as a proxy for accent
            phr = f0 - acc;
            
            % Rectify accent
            acc(acc < 0) = 0;
            
            % Not implemented yet
            fb = zeros(size(f0));
        end
        
        function PlotCrepeActivation(ax, S, F, T, f0)
            % Plot Crepe network output as a heatmap and overlay with thresholded pitch contours
            % 
            %   PlotCrepeActivation(ax, S, F, T)
            %   PlotCrepeActivation(ax, S, F, T, f0)
            % 
            % Inputs
            %   ax          Axes to plot the spectrogram in.
            %   S           A m-by-n matrix of magnitude values. m is the number of frequency bins, 
            %               n is the number of time bins.
            %   F           A m-element vector of frequencies in Hz.
            %   T           A n-element vector of timestamps in second.
            %   f0          A n-element vector of pitch contour in Hz.
            %
            % See also NP.Audio.ComputeMelSpectrogram
            
            % Resample activation S with even Mel frequency bins from F(1) to 500Hz
            F = hz2mel(F);
            Fq = flip(hz2mel(500) : -5 : F(1));
            S = interp1(F, S', Fq', 'linear');
            
            % Plot heatmap
            imagesc(ax, T, Fq, -S);
            if exist('f0', 'var') && ~isempty(f0)
                plot(ax, T, hz2mel(f0), 'Color', 'b', 'LineWidth', 2);
            end
            
            % Get tick frequency
            nFq = numel(Fq);
            nTick = 5;
            tickInd = 1 : ceil(nFq/nTick) : nFq;
            tickVal = Fq(tickInd);
            
            ax.XLim = T([1 end]);
            ax.YDir = 'normal';
            ax.YTick = tickVal;
            ax.YTickLabel = round(mel2hz(tickVal));
            ax.YLabel.String = 'Frequency (Hz)';
            colormap(ax, "gray");
%             ax.CLim = [-1 0];
            MPlot.Axes(ax);
        end
        
        % One-time analyses
        function ComparePitchMethods(ax, se, trialIdx, varargin)
            % 
            
            p = inputParser();
            p.addOptional('tWin', [], @isnumeric);
            p.addParameter('Background', 'mel', @(x) ischar(x) || isstring(x));
            p.parse(varargin{:});
            tWin = p.Results.tWin;
            bg = lower(p.Results.Background);
            
            tbNames = {'ni', 'pitch', 'pitch2'};
            tbs = cell(size(tbNames));
            for i = 1 : numel(tbNames)
                if ~ismember(tbNames{i}, se.tableNames)
                    continue
                end
                if ~isempty(tWin)
                    % Reslice data if time window is provided
                    tbs{i} = se.SliceTimeSeries(tbNames{i}, tWin, trialIdx, 'Fill', 'bleed');
                else
                    % Take the trial as is
                    tb = se.GetTable(tbNames{i});
                    tbs{i} = tb(trialIdx,:);
                end
            end
            niTb = tbs{1};
            
            % Plot background
            switch bg
                case 'mel'
                    % Mel spectrogram
                    melTb = se.SliceTimeSeries('mel', tWin, trialIdx, 'Fill', 'bleed');
                    S = melTb.mic{1};
                    F = se.userData.melMeta.F;
                    isLow = F < 500; % zoom in to the lower freq range
                    F = F(isLow);
                    S = S(:, isLow);
                    T = melTb.time{1};
                    NP.Audio.PlotMelSpectrogram(ax, S, F, T);
                case 'pitch'
                    % Salience map from pitchnn
                    [~, T, F, S] = NP.Pitch.ComputePitch(niTb.mic{1}, niTb.time{1}, 'pitchnn');
                    NP.Pitch.PlotPitchMap(ax, S, F, T);
            end
            hold(ax, 'on');
            
            % Plot NCF pitch contours
            [f0, t] = NP.Pitch.ComputePitch(niTb.mic{1}, niTb.time{1}, 'NCF');
            m = ~isnan(f0);
            f0(m) = hz2mel(f0(m));
            plot(ax, t, f0, 'Color', 'g', 'LineWidth', 1);
            
            % Plot STRAIGHT pitch contours
            f0Tb = tbs{3};
            t = f0Tb.time{1};
            f0 = f0Tb.pitch{1};
            m = ~isnan(f0);
            f0(m) = hz2mel(f0(m));
            plot(ax, t, f0, 'Color', 'b', 'LineWidth', 1);
            
            % Plot old pitch contours
            f0Tb = tbs{2};
            t = f0Tb.time{1};
            f0 = f0Tb.pitch{1};
            m = ~isnan(f0);
            f0(m) = hz2mel(f0(m));
            plot(ax, t, f0, 'Color', 'r', 'LineWidth', 1);
            
            if ~isempty(tWin)
                ax.XLim = tWin;
            end
            ax.YLim(1) = 0;
            ax.XLabel.String = 'Time (s)';
            MPlot.Axes(ax);
        end
        
        function ComparePitchRefinement(ax, se, trialIdx, varargin)
            % 
            
            p = inputParser();
            p.addOptional('tWin', [], @isnumeric);
            p.addParameter('Background', 'mel', @(x) ischar(x) || isstring(x));
            p.parse(varargin{:});
            tWin = p.Results.tWin;
            bg = lower(p.Results.Background);
            
            tbNames = {'pitch2'};
            tbs = cell(size(tbNames));
            for i = 1 : numel(tbNames)
                if ~ismember(tbNames{i}, se.tableNames)
                    continue
                end
                if ~isempty(tWin)
                    % Reslice data if time window is provided
                    tbs{i} = se.SliceTimeSeries(tbNames{i}, tWin, trialIdx, 'Fill', 'bleed');
                else
                    % Take the trial as is
                    tb = se.GetTable(tbNames{i});
                    tbs{i} = tb(trialIdx,:);
                end
            end
            
            % Plot background
            switch bg
                case 'mel'
                    % Mel spectrogram
                    melTb = se.SliceTimeSeries('mel', tWin, trialIdx, 'Fill', 'bleed');
                    S = melTb.mic{1};
                    F = se.userData.melMeta.F;
                    isLow = F < 500; % zoom in to the lower freq range
                    F = F(isLow);
                    S = S(:, isLow);
                    T = melTb.time{1};
                    NP.Audio.PlotMelSpectrogram(ax, S, F, T);
                case 'pitch'
                    % Salience map from pitchnn
                    niTb = se.SliceTimeSeries('ni', tWin, trialIdx, 'Fill', 'bleed');
                    [~, T, F, S] = NP.Pitch.ComputePitch(niTb.mic{1}, niTb.time{1}, 'pitchnn');
                    NP.Pitch.PlotPitchMap(ax, S, F, T);
            end
            hold(ax, 'on');
            
            % Plot STRAIGHT pitch contours
            f0Tb = tbs{1};
            t = f0Tb.time{1};
            f0 = f0Tb.pitch{1};
            m = ~isnan(f0);
            f0(m) = hz2mel(f0(m));
            plot(ax, t, f0, 'Color', 'b', 'LineWidth', 1);
            
            % Refined pitch
            f0Tb = tbs{1};
            t = f0Tb.time{1};
            f0 = f0Tb.pitch{1};
            tt = se.GetTable('taskTime');
            tt = tt(trialIdx,:);
            tg = cat(1, tt.stim, tt.prod{1});
            f0 = NP.Pitch.RefinePitch(f0, t, tg);
            
            m = ~isnan(f0);
            f0(m) = hz2mel(f0(m));
            plot(ax, t, f0, 'Color', 'm', 'LineWidth', 1);
            
            if ~isempty(tWin)
                ax.XLim = tWin;
            end
            ax.YLim(1) = 0;
            ax.XLabel.String = 'Time (s)';
            MPlot.Axes(ax);
        end
        
        % Not in use
        function AddPitchTable(se, chanNames, speechEvents)
            % Add pitch table to se
            % 
            %   NP.Audio.AddPitchTable(se)
            % 
            
            if nargin < 2
                chanNames = {'mic', 'speaker1'};
            end
            if nargin < 3
                speechEvents = {'prod', 'stim'};
            end
            
            % Vectorize data
            tbNames = {'taskTime', 'ni'};
            seLite = se.Duplicate(tbNames, false);
            seLite.SliceSession(0, 'absolute');
            
            % Find periods of voicing
            [tt, ni] = seLite.GetTable(tbNames{:});
            
            % Make pitch table
            pitch = table;
            for i = 1 : numel(chanNames)
                % Compute pitch
                cn = chanNames{i};
                [f0, T] = NP.Pitch.ComputePitch(ni.(cn){1}, ni.time{1});
                
                % Find periods without voicing and set pitch to NaN
                sen = tt.(speechEvents{i}){1};
                phn = Cut(Cut(sen)); % cut down to phones
                phn(~phn.IsVoiced) = [];
                m = phn.MaskTimestamps(T);
                f0(m) = NaN;
                
                pitch.time = {T};
                pitch.(cn) = {f0};
            end
            ff = cell2mat(pitch{1,2:end});
            pitch.merged = {mean(ff, 2, 'omitnan')};
            
            seLite.SetTable('pitch', pitch, 'timeSeries', seLite.GetReferenceTime);
            seLite.RemoveTable(tbNames{:});
            
            % Reslice
            tRef = se.GetReferenceTime('ni');
            tSlice = tRef;
            tSlice(1) = 0;
            seLite.SliceSession(tSlice, 'absolute');
            seLite.AlignTime(tRef-tSlice);
            pitch = seLite.GetTable('pitch');
            se.SetTable('pitch', pitch, 'timeSeries', tRef);
        end
        
    end
    
end
