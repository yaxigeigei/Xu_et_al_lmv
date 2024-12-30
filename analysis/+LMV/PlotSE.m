classdef PlotSE
    
    methods(Static)
        % Unit
        function Scalogram(ax, se, trialInd, tWin, uIdx)
            % Plot scalogram with optional panels for source timeseries and frequency scores
            % 
            %   PlotScalogram(S, T, F, y, scores)
            % 
            
            % Slice out spike times data
            st = se.SliceEventTimes('spikeTime', tWin, [], uIdx, 'Fill', 'none');
            
            % Compute PETH
            fs = 100;
            tEdges = tWin(1) : 1/fs : tWin(2);
            y = MNeuro.MeanEventRate(st, tEdges);
            t = MMath.BinEdges2Centers(tEdges);
            
            % Smoothing
            for i = 1 : size(y,2)
                y(:,i) = MNeuro.Filter1(y(:,i), fs, 'gaussian', 0.025);
            end
            
            % Soft normalization to peak response
            y = MMath.Normalize(y, 'minmaxsoft', 5);
            
            % Compute scalogram
            [S, T, F] = LMV.Sparsity.WaveletTransform(y, t);
            
%             % Compute cleaned 
%             [~, S] = LMV.Sparsity.MeasureFrequencyActivation(S);
            
            % Plot scalogram
            pcolor(t, F, S);
            shading(ax, 'flat');
            ax.XLim = tWin;
%             ax.YLim = F([end 1]);
            ax.YLim = [F(end) 8];
%             ax.XLabel.String = "Time (s)";
            ax.YLabel.String = "Frequency (Hz)";
%             colorbar
            ax.YScale = "log";
            ax.YTick = [1 2 4 8 16];
            ax.YTickLabel = ax.YTick;
            MPlot.Axes(ax);
        end
        
        function PETH(ax, se, trialInd, tWin, uIdx)
            % 
            
            % Slice out spike times data
            st = se.SliceEventTimes('spikeTime', tWin, [], uIdx, 'Fill', 'none');
            
            % Compute PETH
            fs = 100;
            tEdges = tWin(1) : 1/fs : tWin(2);
            Y = MNeuro.MeanEventRate(st, tEdges);
            t = MMath.BinEdges2Centers(tEdges);
            
            % Smoothing
            for i = 1 : size(Y,2)
                Y(:,i) = MNeuro.Filter1(Y(:,i), fs, 'gaussian', 0.025);
            end
            
            % Plot 
            plot(ax, t, Y, 'k');
            ax.YLabel.String = "Spk/s";
            ax.XLim = tWin;
            MPlot.Axes(ax);
        end
        
        % Task events
        function Events(ax, se, trialIdx, tWin, colNames)
            % Plot task and speech events of a trial
            % 
            %   PlotEvents(ax, se, trialIdx, tWin, colNames)
            % 
            
            if ~exist('colNames', 'var') || isempty(colNames)
                colNames = {'task', 'speech'};
            end
            colNames = cellstr(colNames);
            
            if ~exist('tWin', 'var') || isempty(tWin)
                tWin = [-Inf Inf];
            end
            tt = se.SliceEventTimes('taskTime', tWin, trialIdx, colNames); % not using 'bleed' otherwise TGEs will get converted to double
            
            % Plot task events
            n = 0;
            cc = lines(100); % just make a bunch
            evtNames = cell(100,1);
            hold(ax, 'on');
            
            for k = 1 : numel(colNames)
                evt = tt.(colNames{k});
                if iscell(evt)
                    evt = cat(1, evt{:});
                end
                
                if isa(evt, 'NP.TGEvent')
                    % Speech event objects
                    L = evt.GetVfield('source');
                    [G, L] = findgroups(L);
                    evt = splitapply(@(x) {x}, evt, G);
                    
                    for i = 1 : numel(L)
                        n = n + 1;
                        evtNames{n} = L{i};
                        evt{i}.Plot(n+[-.5 .5], 'Color', cc(n,:), 'FontSize', 8, 'Parent', ax);
                    end
                    
                elseif isa(evt, 'MSessionExplorer.Event')
                    % Task event objects
                    tOnOff = [evt.GetTfield('tOn') evt.GetTfield('tOff')];
                    L = evt.GetVfield('type');
                    [G, L] = findgroups(L);
                    tOnOff = splitapply(@(x) {x}, tOnOff, G);
                    
                    for i = 1 : numel(L)
                        n = n + 1;
                        evtNames{n} = L{i};
                        plot(ax, tOnOff{i}', [n n]', 'o', 'Color', cc(n,:));
                    end
                    
                elseif isnumeric(evt)
                    % Numeric events
                    n = n + 1;
                    evtNames{n} = colNames{k};
                    plot(ax, evt, n, 'o', 'Color', cc(n,:));
                    
                else
                    warning("'%s' is not a supported data type", class(evt));
                end
            end
            
            m = ~isinf(tWin);
            ax.XLim(m) = tWin(m);
            ax.XLabel.String = 'Time (s)';
            ax.YDir = 'reverse';
            ax.YLim = [0 n+1];
            ax.YTick = 1 : n;
            ax.YTickLabel = evtNames;
            MPlot.Axes(ax);
        end
        
        function EventBoundary(ax, se, trialInd, tWin, colNames, cmapFunc)
            % Plot lines of the same events across trials, typically on top of raster
            % 
            %   PlotEventBoundary(ax, se, trialInd, tWin, colNames)
            % 
            
            if ~exist('cmapFunc', 'var')
                cmapFunc = @lines;
            end
            
            colNames = cellstr(colNames);
            tt = se.SliceEventTimes('taskTime', tWin, trialInd, colNames); % not using 'bleed' otherwise TGEs will get converted to double
            
            y = (1 : se.numEpochs)';
            cc = cmapFunc(numel(colNames));
            
            hold(ax, 'on');
            for k = 1 : numel(colNames)
                evt = tt.(colNames{k});
                if iscell(evt)
%                     evt = cat(1, evt{:});
                    evt = cellfun(@(x) x(1), evt);
                end
                t = double(evt);
                plot(ax, t, y, '-', 'Color', cc(k,:), 'LineWidth', 1);
            end
            
            m = ~isinf(tWin);
            ax.XLim(m) = tWin(m);
            ax.XLabel.String = 'Time (s)';
            ax.YDir = 'reverse';
            ax.YLim = [0 se.numEpochs+1];
            MPlot.Axes(ax);
        end
        
        function BlockBoundary(ax, se, trialInd, tWin, colNames, varargin)
            % Plot seperator between different trial blocks
            % 
            %   PlotBlockBoundary(ax, se, trialInd, tWin, colNames)
            %   PlotBlockBoundary(ax, se, trialInd, tWin, colNames, ..., 'Label', true)
            %   PlotBlockBoundary(ax, se, trialInd, tWin, colNames, ..., LineArgs)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('Color', lines(numel(colNames)), @(x) isnumeric(x) || ischar(x));
            p.addParameter('Label', true, @islogical);
            p.parse(varargin{:});
            isLabel = p.Results.Label;
            cc = p.Results.Color;
            
            tv = se.GetTable('taskValue');
            if ~isempty(trialInd)
                tv = tv(trialInd,:);
            end
            
            if size(cc,1) == 1
                cc = repmat(cc, [numel(colNames) 1]);
            end
            
            hold(ax, 'on');
            for k = 1 : numel(colNames)
                cn = colNames{k};
                val = tv.(cn);
                [C, ~, ic] = unique(val, 'stable');
                d = [1; diff(ic); 1];
                y = find(d) - 0.5;
                
                plot(ax, tWin', [y y]', '-', 'Color', cc(k,:), p.Unmatched);
                
                ym = mean([y(1:end-1) y(2:end)], 2);
                tm = repmat(mean(tWin), size(ym));
                
                if ~isLabel
                    continue
                end
                text(ax, tm, ym, cn+": "+string(C), 'HorizontalAlignment', 'center', 'Interpreter', 'none');
            end
            
            m = ~isinf(tWin);
            ax.XLim(m) = tWin(m);
            ax.XLabel.String = 'Time (s)';
            ax.YTick = 1 : 2 : se.numEpochs;
            ax.YDir = 'reverse';
            ax.YLim = [0 se.numEpochs+1];
            ax.YLabel.String = 'Trial';
            MPlot.Axes(ax);
        end
        
        function EventWindows(ax, se, trialIdx, tWin, colNames, varargin)
            % Plot task and speech event objects for a single trial
            % 
            %   PlotEventWindows(ax, se, trialIdx, tWin, colNames)
            %   PlotEventWindows(..., 'Style', 'block')
            %   PlotEventWindows(..., 'YRange', trialIdx+[-0.5 0.5])
            %   PlotEventWindows(..., 'Colors', lines(numel(colNames)))
            %   PlotEventWindows(..., 'Alpha', 0.2)
            %   PlotEventWindows(..., 'Text', false)
            % 
            
            if ~exist('tWin', 'var') || isempty(tWin)
                tWin = [-Inf Inf];
            end
            
            colNames = cellstr(colNames);
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('Style', 'block', @(x) ischar(x) || isstring(x));
            p.addParameter('YRange', trialIdx + [-.5 .5], @(x) isnumeric(x) && numel(x)==2);
            p.addParameter('Colors', lines(numel(colNames)));
            p.addParameter('Alpha', .2, @isnumeric);
            p.addParameter('Text', false, @islogical);
            p.parse(varargin{:});
            style = p.Results.Style;
            cc = p.Results.Colors;
            yRange = p.Results.YRange;
            alpha = p.Results.Alpha;
            isText = p.Results.Text;
            
            if numel(alpha) == 1 && numel(colNames) > 1
                alpha = repmat(alpha, size(colNames));
            end
            
            % Get data
            tt = se.SliceEventTimes('taskTime', tWin, trialIdx, colNames); % not using 'bleed' otherwise TGEs will get converted to double
            
            % Plot task events
            hold(ax, 'on');
            
            for k = 1 : numel(colNames)
                evt = tt.(colNames{k});
                if iscell(evt)
                    evt = cat(1, evt{:});
                end
                
                if isa(evt, 'NP.TGEvent')
                    % Speech event objects
                    tOnOff = [evt.GetTfield('tmin') evt.GetTfield('tmax')];
                elseif isa(evt, 'MSessionExplorer.Event')
                    % Task event objects
                    tOnOff = [evt.GetTfield('tOn') evt.GetTfield('tOff')];
                else
                    % Numeric values
                    tOnOff = evt;
                end
                
                switch style
                    case 'block'
                        MPlot.Blocks(tOnOff, yRange, cc(k,:), 'FaceAlpha', alpha(k), 'Parent', ax);
                    case 'line'
                        xx = repmat(tOnOff(:)', [2 1]);
                        yy = repmat(yRange(:), [1 size(xx,2)]);
                        plot(ax, xx, yy, 'Color', [cc(k,:) alpha(k)], p.Unmatched);
                    otherwise
                        error("'%s' is not a valid style.", style);
                end
                
                if isText
                    ss = repelem(string(colNames{k}), size(tOnOff,1))';
                    x = mean(tOnOff, 2);
                    y = repelem(mean(yRange), numel(x))';
                    text(ax, x, y, ss, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                end
            end
            
            m = ~isinf(tWin);
            ax.XLim(m) = tWin(m);
            ax.XLabel.String = 'Time (s)';
%             ax.YLim = yRange;
            MPlot.Axes(ax);
        end
        
        % Phonetics
        function TGE(ax, se, trialInd, tWin, varargin)
            % Plot phonetic marks across trials that are aligned in time
            % 
            %   PlotTGE(ax, se, trialInd, tWin, tierName)
            %   PlotTGE(ax, se, trialInd, tWin, tierName, colNames)
            % 
            % Inputs
            %   ax          Axes to plot in.
            %   se          MSessionExplorer object
            %   trialInd    Indices of trials to plot for.
            %   tWin        Time window. 
            %   tierName    At what tier to plot the labels.
            %   colNames    Column names of the taskTime table to read NP.TGEvent objects from.
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('tierName', [], @(x) ischar(x) || isstring(x));
            p.addOptional('colNames', 'prod', @(x) ischar(x) || iscellstr(x) || isstring(x));
            p.parse(varargin{:});
            tierName = string(p.Results.tierName);
            colNames = cellstr(p.Results.colNames);
            
            tt = se.SliceEventTimes('taskTime', tWin, trialInd, colNames); % not using 'bleed' otherwise TGEs will get converted to double
            nTrial = height(tt);
            cc = lines(nTrial);
            
            for k = 1 : nTrial
                for i = 1 : width(tt)
                    % Get TGEvent object(s) of the current trial
                    if iscell(tt.(i))
                        tge = tt.(i){k};
                    else
                        tge = tt.(i)(k);
                    end
                    if isnan(tge)
                        continue
                    end
                    
                    % Iterate through tiers
                    nTier = tgGetNumberOfTiers(tge(1));
                    for T = 0 : nTier
                        % Plot marks at the phones tier, or last tier if phones tier is not avalable
                        if isempty(tierName) || tge(1).tier1 == tierName || T == nTier
                            tge.Plot(k+[-.5 .5], 'Color', cc(k,:), 'Parent', ax, p.Unmatched);
                            % tge.Plot(k+[-.5 .5], 'Style', 'patch', 'FontSize', 0, 'Parent', ax);
                            break
                        end
                        if T < nTier
                            tge = Cut(tge);
                        end
                    end
                end
            end
            
            ax.XLim = [min(tWin(:,1)) max(tWin(:,2))];
            ax.XLabel.String = 'Time (s)';
            ax.YDir = 'reverse';
            ax.YTick = 1:2:nTrial;
            ax.YLabel.String = 'Trials';
            ax.YLim = [0 nTrial+1];
            MPlot.Axes(ax);
        end
        
        function TGEHier(ax, se, trialIdx, tWin, colNames, varargin)
            % Plot a hierarchy of speech labels for a given trial
            % 
            %   PlotTGEHier(ax, se, trialInd, tWin, colNames)
            % 
            % Inputs
            %   ax          Axes to plot in.
            %   se          MSessionExplorer object
            %   trialIdx    Index of trial to plot for.
            %   tWin        Time window. 
            %   colNames    Column names of the taskTime table to read NP.TGEvent objects from.
            % 
            
            if ~exist('colNames', 'var') || isempty(colNames)
                colNames = {'stim', 'prod'};
            end
            colNames = cellstr(colNames);
            
            % Get data
            if isempty(tWin)
                tt = se.GetTable('taskTime');
                tt = tt(trialIdx,:);
            else
                tt = se.SliceEventTimes('taskTime', tWin, trialIdx, colNames); % not using 'bleed' otherwise TGEs will get converted to double
            end
            
            % Collect objects from all columns
            tge = cell(size(colNames));
            for i = 1 : numel(colNames)
                n = colNames{i};
                if iscell(tt.(n))
                    e = tt.(n){1}(:);
                else
                    e = tt.(n)(1);
                end
                tge{i} = e;
            end
            tge = cat(1, tge{~cellfun(@(x) all(isnan(x)), tge)});
            if isempty(tge) || all(isnan(tge))
                return
            end
            
            % Plot tiers
            tge.PlotTiers('Parent', ax, varargin{:});
            
            if ~isempty(tWin)
                ax.XLim = tWin;
            end
            ax.XLabel.String = 'Time (s)';
            ax.YDir = 'reverse';
            ax.YTick = 1 : tgGetNumberOfTiers(tge(1))+1;
            ax.YTickLabel = {'Sent', 'Word', 'Phone'};
            MPlot.Axes(ax);
        end
        
        % Acoustics
        function MelSpectrogram(ax, se, trialIdx, varargin)
            % Plot Mel spectralgram for a given trial with optional waveform overlay
            % 
            %   PlotMelSpectrogram(ax, se, trialIdx)
            %   PlotMelSpectrogram(ax, se, trialIdx, tWin)
            %   PlotMelSpectrogram(..., 'WaveformHeight', 0.15)
            %   PlotMelSpectrogram(..., 'ShowIntensity', false)
            % 
            % Inputs
            %   tWin                1-by-2 numeric array for the time window to plot.
            %   'WaveformHeight'    The height of waveform as a fraction of the spectrogram. 
            %                       Default is 0.15. Waveform will not be plotted if set to zero.
            %   'ShowIntensity'     Whether or not to show intensity features including envelope, 
            %                       peak evelope, and peak rate. Default is false.
            % 
            
            p = inputParser();
            p.addOptional('tWin', [], @isnumeric);
            p.addParameter('WaveformHeight', 0.15, @isnumeric);
            p.addParameter('ShowIntensity', false, @islogical);
            p.parse(varargin{:});
            tWin = p.Results.tWin;
            hWf = p.Results.WaveformHeight;
            isInten = p.Results.ShowIntensity;
            
            % Get data tables
            tbNames = {'mel', '', ''};
            if hWf
                tbNames{2} = 'ni';
                if isInten
                    tbNames{3} = 'inten';
                end
            end
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
            
            % Always plot spectrogram
            mel = tbs{1};
            if isempty(mel)
                fprintf("Cannot plot Mel spectrogram because the se does not have 'mel' table.\n");
                return
            end
            S = mel.(2){1};
            F = se.userData.melMeta.F;
            T = mel.time{1};
            NP.Audio.PlotMelSpectrogram(ax, S, F, T);
            hold(ax, 'on');
            melMax = max(hz2mel(F));
            
            ni = tbs{2};
            if ~isempty(ni)
                % Plot microphone waveform
                tNI = ni.time{1};
                w = ni.mic{1};
                tf = @(y) y/max(abs(w)) * hWf/2*melMax + (1-hWf/2)*melMax;
%                 w = w/max(abs(w)) * hWf/2*melMax; % normalize to the desired height
%                 w = w + (1-hWf/2)*melMax; % lift it to the upper boarder of the spectrogram
                w = tf(w);
                plot(ax, tNI, w, 'k');
                
                it = tbs{3};
                if ~isempty(it)
                    % Plot envelope, its peaks, and the peaks of d(env)/dt
                    tIt = it.time{1};
                    plot(ax, tIt, tf(it.env{1}), 'r', 'LineWidth', 2);
                    plot(ax, tIt, tf(it.peakEnv{1}/4), 'g', 'LineWidth', 2);
                    plot(ax, tIt, tf(it.peakRate{1}), 'y', 'LineWidth', 2);
                end
            end
            
            if ~isempty(tWin)
                ax.XLim = tWin;
            end
%             ax.YLim = [0 max(yNI)];
            ax.XLabel.String = 'Time (s)';
            MPlot.Axes(ax);
        end
        
        function Intensity(ax, se, trialIdx, tWin)
            % Plot timeseries of intensity features for a given trial
            % 
            %   PlotIntensity(ax, se, trialIdx, tWin)
            % 
            
            if nargin > 3
                % Reslice data if time window is provided
                ni = se.SliceTimeSeries('ni', tWin, trialIdx, 'Fill', 'bleed');
                it = se.SliceTimeSeries('inten', tWin, trialIdx, 'Fill', 'bleed');
            else
                % Take data as is
                [ni, it] = se.GetTable('ni', 'inten');
                ni = ni(trialIdx,:);
                it = it(trialIdx,:);
            end
            
            % Plot microphone waveform
            tNI = ni.time{1};
            wav = ni.mic{1};
            maxVal = max(abs(wav));
            plot(ax, tNI, wav, 'k');
            hold(ax, 'on');
            
            % Plot envelope, its peaks, and the peaks of d(env)/dt
            tIt = it.time{1};
            plot(ax, tIt, it.env{1}, 'r', 'LineWidth', 2);
            plot(ax, tIt, it.peakEnv{1} / 4, 'g', 'LineWidth', 2);
            plot(ax, tIt, it.peakRate{1}, 'y', 'LineWidth', 2);
            
            ax.XLim = tWin;
            ax.XLabel.String = 'Time (s)';
            ax.YLim = [-maxVal maxVal] * 1.2;
            ax.YLabel.String = 'AU';
            MPlot.Axes(ax);
        end
        
        function Pitch(ax, se, trialIdx, varargin)
            % Plot acoustic timeseries for a given trial
            % 
            %   PlotPitch(ax, se, trialIdx)
            %   PlotPitch(ax, se, trialIdx, tWin)
            %   PlotPitch(..., 'Background', 'mel');
            % 
            
            p = inputParser();
            p.addOptional('tWin', [], @isnumeric);
            p.addParameter('Background', 'mel', @(x) ischar(x) || isstring(x));
            p.parse(varargin{:});
            tWin = p.Results.tWin;
            bg = lower(p.Results.Background);
            
            tbNames = {'pitch'};
            switch bg
                case 'mel'
                    tbNames{2} = 'mel';
                case 'ni'
                    tbNames{2} = 'ni';
            end
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
                    melTb = tbs{2};
                    T = melTb.time{1};
                    S = melTb.mic{1};
                    F = se.userData.melMeta.F;
                    isLow = F < 400; % zoom in to the lower freq range
                    NP.Audio.PlotSpectrogram(ax, S(:,isLow), F(isLow), T, 'FTick', 0:50:400);
                    ax.YLim = [50 350];
                case 'pitch'
                    % Salience map from pitchnn
                    niTb = tbs{2};
                    [~, T, F, S] = NP.Pitch.ComputePitch(niTb.mic{1}, niTb.time{1}, 'pitchnn');
                    NP.Pitch.PlotPitchMap(ax, S, F, T);
                    ax.YLim = hz2mel([50 350]);
            end
            hold(ax, 'on');
            
            % Plot pitch contours
            f0Tb = tbs{1};
            t = f0Tb.time{1};
            f0 = f0Tb.F0{1};
            plot(ax, t, f0, 'LineWidth', 2);
            
            if ~isempty(tWin)
                ax.XLim = tWin;
            end
        end
        
        function RelativePitch(ax, se, trialIdx, varargin)
            % Plot acoustic timeseries for a given trial
            % 
            %   PlotRelativePitch(ax, se, trialIdx)
            %   PlotRelativePitch(ax, se, trialIdx, tWin)
            % 
            
            p = inputParser();
            p.addOptional('tWin', [], @isnumeric);
            p.parse(varargin{:});
            tWin = p.Results.tWin;
            
            tbNames = {'pitch'};
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
            
            % Plot pitch contours
            f0Tb = tbs{1};
            t = f0Tb.time{1};
            rF0 = f0Tb.rF0{1};
            drF0 = f0Tb.drF0{1};
            cc = lines;
            
            yyaxis(ax, 'left');
            plot(ax, t, rF0, 'Color', cc(1,:), 'LineWidth', 1); hold(ax, 'on');
            ax.YLabel.String = "Relative F0";
            
            yyaxis(ax, 'right');
            plot(ax, t, drF0, 'Color', cc(2,:), 'LineWidth', 1);
            ax.YLabel.String = "d(relative F0)/dt";
            
            if ~isempty(tWin)
                ax.XLim = tWin;
            else
                ax.XLim = t([1 end]);
            end
            MPlot.Axes(ax);
        end
        
        function BinnedPitch(ax, se, trialIdx, varargin)
            % Plot binary map of binned pitch timeseries for a given trial
            % 
            %   PlotBinnedPitch(ax, se, trialIdx)
            %   PlotBinnedPitch(ax, se, trialIdx, tWin)
            % 
            
            p = inputParser();
            p.addOptional('tWin', [], @isnumeric);
            p.parse(varargin{:});
            tWin = p.Results.tWin;
            
            tbNames = {'pitch'};
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
            
            % Resample 
            f0Tb = tbs{1};
            if isempty(tWin)
                tWin = f0Tb.time{1}([1 end]);
            end
            tEdges = tWin(1) : diff(f0Tb.time{1}(1:2)) : tWin(2);
            f0Tb = se.ResampleTimeSeries(f0Tb, tEdges);
            t = f0Tb.time{1};
            brF0 = f0Tb.brF0{1};
            
            % Plot pitch contours
            imagesc(ax, t, [], brF0');
            colormap(ax, 'gray');
            ax.YDir = 'normal';
            ax.XLim = tWin;
            ax.YLabel.String = "Relative F0 bins";
            MPlot.Axes(ax);
        end
        
        function Fujisaki(ax, se, trialIdx, varargin)
            % Plot voicing, phrase, and accent for a given trial
            % 
            %   PlotFujisaki(ax, se, trialIdx)
            %   PlotFujisaki(ax, se, trialIdx, tWin)
            % 
            
            p = inputParser();
            p.addOptional('tWin', [], @isnumeric);
            p.parse(varargin{:});
            tWin = p.Results.tWin;
            
            tbNames = {'pitch'};
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
            
            % Prepare variables
            f0Tb = tbs{1};
            t = f0Tb.time{1};
            voi = f0Tb.voicing{1};
            rF0 = f0Tb.rF0{1};
            phr = f0Tb.phrase{1};
            acc = f0Tb.accent{1};
            
            voiRanges = MMath.Logical2Bounds(voi);
            voiRanges = t(voiRanges);
            yRange = [0 max(rF0, [], 'omitnan')*1.1];
            
            % Plot timeseries
            cc = lines;
            MPlot.Blocks(voiRanges, yRange, [0 0 0], 'FaceAlpha', 0.1, 'Parent', ax);
            hold(ax, 'on');
            plot(ax, t, rF0, 'Color', cc(1,:), 'LineWidth', 2);
            plot(ax, t, phr, 'Color', cc(3,:), 'LineWidth', 1);
            plot(ax, t, acc, 'Color', cc(4,:), 'LineWidth', 1);
            if ~isempty(tWin)
                ax.XLim = tWin;
            else
                ax.XLim = t([1 end]);
            end
%             ax.YLim = yRange;
            ax.YLim = [0 10];
            ax.YLabel.String = "Relative F0";
            MPlot.Axes(ax);
        end
        
        function MelSpectrograms(ax, se, trialInd, varargin)
            % Plot a stack of Mel spectralgrams with optional waveform overlay
            % 
            %   PlotMelSpectrograms(ax, se, trialInd)
            %   PlotMelSpectrograms(ax, se, trialInd, tWin)
            %   PlotMelSpectrograms(..., 'Waveform', true)
            % 
            
            if ~ismember('mel', se.tableNames)
                return
            end
            
            p = inputParser();
            p.addOptional('tWin', [], @isnumeric);
            p.addParameter('VariableName', 'mic', @(x) ischar(x) || isstring(x));
            p.addParameter('Waveform', true, @islogical);
            p.parse(varargin{:});
            vn = p.Results.VariableName;
            tWin = p.Results.tWin;
            isWf = p.Results.Waveform;
            
            if ~isempty(tWin)
                % Reslice data if time window is provided
                mel = se.SliceTimeSeries('mel', tWin, trialInd, 'Fill', 'none');
                ni = se.SliceTimeSeries('ni', tWin, trialInd, 'Fill', 'none');
            else
                % Take data as is
                [mel, ni] = se.GetTable('mel', 'ni');
                if exist('trialInd', 'var') && ~isempty(trialInd)
                    mel = mel(trialInd,:);
                    ni = ni(trialInd,:);
                end
            end
            
            yBase = 0;
            
            for i = height(mel) : -1 : 1
                % Plot spectrogram
                S = mel.(vn){i};
                F = se.userData.melMeta.F;
                tMel = mel.time{i};
                NP.Audio.PlotMelSpectrogram(ax, S, F, tMel, [], yBase); hold(ax, 'on')
                hold(ax, 'on');
                melMax = max(hz2mel(F));
                
                if isWf
                    % Plot microphone waveform
                    tNI = ni.time{i};
                    wf = ni.(vn){i};
                    r = 0.15;
                    wf = wf/max(abs(wf)) * r/2*melMax; % make the height of waveform a fraction of the spectrogram
                    wf = wf + (1-r/2)*melMax; % lift it to the upper boarder of the spectrogram
                    plot(ax, tNI, wf+yBase, 'k');
                end
                
                yBase = yBase + melMax*1.05;
            end
            
            if ~isempty(tWin)
                ax.XLim = tWin;
            end
            ax.YLim = [0 yBase-melMax*0.05];
            ax.XLabel.String = 'Time (s)';
            MPlot.Axes(ax);
        end
        
        % Facial tracking
        function FacialFeatures(ax, se, trialIdx, tWin)
            % Plot Mel spectralgram for a given trial
            
%             varNames = {'mouthHeight', 'mouthWidth', 'chin2nose', 'leftPupilX', 'rightPupilX'};
            varNames = {'mouthHeight', 'mouthWidth', 'chin2nose'};
            
            if nargin > 3
                % Reslice data if time window is provided
                lm = se.SliceTimeSeries('landmark', tWin, trialIdx, varNames, 'Fill', 'bleed');
            else
                % Take data as is
                lm = se.GetTable('landmark');
                lm = lm(trialIdx, [{'time'} varNames]);
            end
            
            edgeLen = double(se.userData.landmarkMeta.img_size(1));
            t = lm.time{1};
            X = cell2mat(lm{1,2:end}) * edgeLen;
            plot(ax, t, X);
            
%             legend(ax, {'Mouth height', 'Mouth width'}, 'Location', 'west');
            ax.XLim = tWin;
            ax.XLabel.String = 'Time (s)';
            ax.YLabel.String = 'Distance (px)';
            MPlot.Axes(ax);
        end
        
        function Saccade(ax, se, trialInd, tWin)
            % Plot timeseries of saccade velocity for given trials
            
            varNames = {'pupilV'};
            
            if nargin > 3
                % Reslice data if time window is provided
                lm = se.SliceTimeSeries('landmark', tWin, trialInd, varNames, 'Fill', 'none');
            else
                % Take data as is
                lm = se.GetTable('landmark');
                lm = lm(trialInd, [{'time'} varNames]);
            end
            lm.color = lines(height(lm));
            
            vMax = cellfun(@max, lm.pupilV);
            vMaxMed = median(vMax, 'omitnan');
            yOffset = (0:height(lm)-1) * vMaxMed;
            
            lm = flip(lm);
            MPlot.PlotTraceLadder(lm.time, lm.pupilV, yOffset, 'ColorArray', lm.color, 'Parent', ax);
            
            ax.XLim = tWin;
            ax.XLabel.String = 'Time (s)';
            ax.YLabel.String = 'Velocity (px/s)';
            MPlot.Axes(ax);
        end
        
        % Time morphing
        function Alignment(seTb, varargin)
            % Plot alignment scores, speed map, words, and spectrograms
            % 
            %   PlotAlignment(seTb)
            % 
            % Inputs
            %   seTb            A table of se objects and other information for task conditions.
            %   
            
            % Find indices of the template trials and their time windows
            for k = height(seTb) : -1 : 1
                [tt, tv] = seTb.se(k).GetTable('taskTime', 'taskValue');
                tr = find(tv.trialNum == tv.tempTrialNum(1), 1);
                trialInd(k,:) = tr;
                tWin(k,:) = [tt.matchOn(tr) tt.matchOff(tr)] + [-0.5 0.5];
            end
            
            % Plot parameters
            numRows = 8;
            numCols = height(seTb);
            
            % Features
            for i = 1 : height(seTb)
                se = seTb.se(i);
                tr = trialInd(i);
                w = tWin(i,:);
                k = 0;
                
                % Alignment scores
                k = k + 1;
                loc = zeros(numRows, numCols);
                loc(k,i) = 1;
                loc = find(loc');
                ax = MPlot.Axes(numRows, numCols, loc); cla(ax)
                NP.TaskBaseClass.PlotAlignmentScores(ax, se);
                ax.Title.String = seTb.stimText(i);
                
                % Speech map
                k = k + 1;
                loc = zeros(numRows, numCols);
                loc(k,i) = 1;
                loc = find(loc');
                ax = MPlot.Axes(numRows, numCols, loc); cla(ax)
                NP.TaskBaseClass.PlotSpeedMap(ax, se, [], w);
                ax.XLabel.String = [];
                
                % Speech labels
                k = k + 1;
                loc = zeros(numRows, numCols);
                loc(k,i) = 1;
                loc = find(loc');
                ax = MPlot.Axes(numRows, numCols, loc); cla(ax)
                NP.TaskBaseClass.PlotAlignedWords(ax, se, [], w, 'combPhn');
                ax.XLabel.String = [];
                
                % Mel spectrogram
                k = k + 1;
                loc = zeros(numRows, numCols);
                loc(k:end,i) = 1;
                loc = find(loc');
                ax = MPlot.Axes(numRows, numCols, loc); cla(ax)
                NP.TaskBaseClass.PlotMelSpectrograms(ax, se, [], w);
            end
        end
        
        function AlignmentScores(ax, se, trialInd)
            % Plot alignment scores of stim-production and template-repeat
            
            tv = se.GetTable('taskValue');
            if exist('trialInd', 'var') && ~isempty(trialInd)
                tv = tv(trialInd,:);
            end
            
            nTrial = height(tv);
            y = 1 : nTrial;
            x1 = min(tv.alignScore, 1);
            x2 = cellfun(@(x) x.nscore, tv.alignInfo);
            
            cc = lines(2);
            plot(x1, y, 'o', 'Color', cc(1,:)); hold(ax, 'on');
            plot(x2, y, '*', 'Color', cc(2,:));
            legend(ax, {'stim', 'template'}, 'Location', 'west');
            ax.XLim = [0 1];
            ax.XLabel.String = 'Frac. of max score';
            ax.YDir = 'reverse';
            ax.YTick = 1:2:nTrial;
            ax.YLabel.String = 'Trials';
            MPlot.Axes(ax);
        end
        
        function AlignedWords(ax, se, trialInd, tWin, colName)
            % Plot phonetic marks across trials that are aligned in time
            
            if ~exist('colName', 'var') || isempty(colName)
                colName = 'prod';
            end
            
            tt = se.SliceEventTimes('taskTime', tWin, trialInd, {colName}); % not using 'bleed' otherwise TGEs will get converted to double
            nTrial = height(tt);
            cc = lines(nTrial);
            
            for k = 1 : nTrial
                % Get TGEvent object(s) of the current trial
                if iscell(tt.(colName))
                    tge = tt.(colName){k};
                else
                    tge = tt.(colName)(k);
                end
                
                % Check each tier
                nTier = tgGetNumberOfTiers(tge(1));
                for t = 0 : nTier
                    % Plot marks at the phones tier, or last tier if phones tier is not avalable
                    if tge(1).tier1 == "phone" || t == nTier
                        tge.Plot(k+[-.5 .5], 'Color', cc(k,:), 'Parent', ax);
                        break
                    end
                    if t < nTier
                        tge = Cut(tge);
                    end
                end
            end
            
            ax.XLim = tWin;
            ax.XLabel.String = 'Time (s)';
            ax.YDir = 'reverse';
            ax.YTick = 1:2:nTrial;
            MPlot.Axes(ax);
        end
        
        function SpeedMap(ax, se, trialInd, tWin)
            % Plot phonetic marks across trials that are aligned in time
            
            if ~ismember('morph', se.tableNames)
                return
            end

            t = tWin(1):0.01:tWin(2);
            tEdges = MMath.BinCenters2Edges(t);
            tt = se.ResampleTimeSeries('morph', tEdges, trialInd);
            nTrial = height(tt);
            
            y = 1 : nTrial;
            S = cat(2, tt.speed{:})';
            imagesc(ax, t, y, S);
            
            ax.XLim = tWin;
            ax.XLabel.String = 'Time (s)';
            ax.YDir = 'reverse';
            ax.YTick = 1:2:nTrial;
            ax.YLabel.String = 'Trials';
            colormap(ax, MPlot.PolarMap);
            caxis(ax, [0 2]);
            MPlot.Axes(ax);
        end
        
    end
    
end
