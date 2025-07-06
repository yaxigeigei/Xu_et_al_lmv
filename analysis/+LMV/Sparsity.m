classdef Sparsity
    
    methods(Static)
        % 
        
        function PlotKurtosis(ce, S)
            % Plot the kurtosis of population activity
            
            % Reshape vector
            respTb = ce.GetTable("resp");
            t = respTb.time{1};
            S = reshape(S, [], height(respTb));
            
            plot(t, S); hold on
            ax = gca;
            
            colNames = {'cue1', 'stim', 'prod', 'cue3'};
            cc = LMV.Param.GetTaskPhaseColors(colNames);
            NP.TaskBaseClass.PlotEventWindows(ax, ce, 1, t([1 end])', colNames, 'YRange', ax.YLim, 'Colors', cc, 'Alpha', 0.1, 'Text', false);
            
            ax.Title.String = "Population sparseness";
            % ax.XTick = [];
            % ax.YTick = [];
            % ax.XLabel.String = [];
            % ax.YLabel.String = [];
            % axis off
        end
        
        
        % not in use
        function [scores, F] = ComputeScores(Y, fs)
            % Compute sparsity scores over a range of frequencies
            % 
            %   [scores, F] = ComputeScores(Y, fs)
            % 
            % Inputs
            %   Y           A t-by-u matrix of spike rates for t time points and u units.
            %   fs          Sampling rate in Hz.
            % Outputs
            %   scores      A u-by-f matrix of sparsity scores measuring the average power 
            %               (above threshold) as a function of f transient frequencies.
            %   F           A 1-by-f vector of frequency values.
            % 
            
            % Smoothing
            if isvector(Y)
                Y = Y(:);
            end
            for i = 1 : size(Y,2)
                Y(:,i) = MNeuro.Filter1(Y(:,i), fs, 'gaussian', 0.025);
            end
            
            % Soft normalization to peak response
            Y = MMath.Normalize(Y, 'minmaxsoft', 5);
            
            % Compute scores with wavelet decomposition
            [S, T, F] = LMV.Sparsity.WaveletTransform(Y, fs);
            scores = LMV.Sparsity.MeasureFrequencyActivation(S);
            scores = scores';
            F = F';
        end
        
        function [S, T, F] = WaveletTransform(Y, tOrFs)
            % Perform wavelet transformation on signals
            % 
            %   [S, T, F] = WaveletTransform(Y, t)
            %   [S, T, F] = WaveletTransform(Y, fs)
            % 
            % Inputs
            %   Y           A t-by-u matrix of spike rates for t time points u units.
            %   t           Timestamps of 
            %   fs          Sampling rate in Hz.
            % Outputs
            %   S           A f-by-t-by-u array of scalograms for f frequencies, t time points, and u units.
            %   T           A t-by-1 vector of timestamps for S.
            %   F           A f-by-1 vector of frequency values for S.
            % 
            
            % Parameters
            frequencyLimits = [0.2 16]; % Hz
            voicesPerOctave = 8;
            timeBandWidth = 10;
            if isscalar(tOrFs)
                fs = tOrFs;
            else
                fs = 1 / diff(tOrFs(1:2));
            end
            
            if isvector(Y)
                Y = Y(:);
            end
            
            % Limit the cwt frequency limits
            frequencyLimits(1) = max(frequencyLimits(1), ...
                cwtfreqbounds(numel(Y), fs, 'TimeBandWidth', timeBandWidth));
            
            nUnit = size(Y, 2);
            S = cell(nUnit, 1);
            for i = 1 : nUnit
                [WT, F] = cwt(Y(:,i), fs, ...
                    'VoicesPerOctave', voicesPerOctave, ...
                    'TimeBandWidth', timeBandWidth, ...
                    'FrequencyLimits', frequencyLimits);
                S{i} = abs(WT);
            end
            
            S = cat(3, S{:});
            T = (1:size(S,2))' / fs;
            F = F(:);
        end
        
        function [act, S] = MeasureFrequencyActivation(S)
            % Detect and measure activations as a function of frequency from scalogram(s)
            % 
            %   [act, S] = MeasureScalogramActivation(S)
            % 
            % Input
            %   S           A f-by-t-by-u array of scalograms for f frequencies, t time points, and u units.
            % Outputs
            %   act         A f-by-u matrix of sparsity scores measuring the average power (above threshold) as 
            %               a function of f transient frequencies.
            %   S           Similar to input S except that 1) values at each frequency are subtracted by median,
            %               and 2) all subthreshold or negative values are set to zeros.
            % 
            
            % Find subthrehold values
            [~, ~, th, med] = isoutlier(S, 2); % default is 3*MAD
            isSub = S < th;
            
            % Set subthrehold values to zero after background subtraction
            S = S - med;
            S(isSub) = 0;
            
            % Compute mean of squares
            isAllSub = all(isSub, 1); % size 1-by-t-by-u
            nSupra = sum(~isAllSub, 2); % size 1-by-1-by-u
            act = sum(S.^2, 2) ./ nSupra; % size f-by-1-by-u
            act = permute(act, [1 3 2]); % size f-by-u
        end
        
        function PlotScalogram(S, T, F, y, scores)
            % Plot scalogram with optional panels for source timeseries and frequency scores
            % 
            %   PlotScalogram(S, T, F, y, scores)
            % 
            
            ax = nexttile;
            pcolor(T, F, S)
            shading flat
            xlabel("Time (s)")
            ylabel("Frequency (Hz)")
            colorbar
            ax.YScale = "log";
            ax.YTick = [1 2 4 8 16];
            ax.YTickLabel = ax.YTick;
            MPlot.Axes(ax);
            
            if exist('y', 'var') && ~isempty(y)
                ax = nexttile("north");
                for i = 1 : size(y,2)
                    plot(T, y(:,i)); hold on
                end
                ax.XLim = T([1 end]);
                ax.XTickLabel = [];
                ax.YLabel.String = "Spike rate";
                MPlot.Axes(ax);
            end
            
            if exist('scores', 'var') && ~isempty(scores)
                ax = nexttile("east");
                for i = 1 : size(scores,2)
                    plot(scores(:,i), F(:)); hold on
                end
%                 ax.XLim = [0 1];
                ax.YScale = "log";
                ax.YLim = F([end 1]);
                ax.YTick = [1 2 4 8 16];
                ax.YTickLabel = ax.YTick;
                ax.Title.String = "Scores";
                MPlot.Axes(ax);
                ax.Box = 'on';
            end
        end
        
        function ClickHeatmap(hImg, eventData)
            % 
            
            % Clears the selection box
            ax = hImg.Parent;
            if isfield(ax.UserData, 'hRect')
                delete(ax.UserData.hRect);
            end
            if eventData.Button ~= 1
                return
            end
            clusTb = ax.UserData.clusTb;
            
            % Find the range of units
            k = round(eventData.IntersectionPoint(2));
            ind = (0 : 9) + k;
            ind = MMath.Bound(ind, [1 height(clusTb)]);
            ind = unique(ind);
            
            % Plot a selection box on the heatmap
            ax.UserData.hRect = rectangle(ax, 'Position', [ax.XLim(1) ind(1)-0.5 diff(ax.XLim) ind(end)-ind(1)+1], 'EdgeColor', 'r');
            
            % Plot raster
            MPlot.Figure(123); clf
            LMV.Overview.SessionFromCache(clusTb.clusId(ind), 'DataSource', 'm2');
        end
        
    end
    
end