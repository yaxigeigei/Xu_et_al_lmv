classdef Chunk
    
    methods(Static)
        function gapData = AnalyzeWordBoundaries(seSen, morphTb, sen)
            % Analyze word boundaries for chunking patterns
            % 
            %   gapData = AnalyzeWordBoundaries(seSen, morphTb, sen)
            % 
            % This method extracts the word boundary analysis from the original script.
            
            words = Cut(sen);
            phones = arrayfun(@(x) Cut(Trim(x)), words, 'Uni', false);
            gapWins = cellfun(@(x,y) [x(end), y(1)], phones(1:end-1), phones(2:end), 'Uni', false);
            gapWins = double(cat(1, gapWins{:}));
            
            nGaps = size(gapWins, 1);
            gapData = cell(nGaps, 1);
            
            for i = 1 : nGaps
                % Slice original timestamps
                gapTb = seSen.SliceTimeSeries(morphTb, gapWins(i,:));
                if any(cellfun(@isempty, gapTb.time))
                    continue
                end
                
                % Get the speed data for this window
                x = (1 : height(gapTb))';
                y = cellfun(@(x) x(end) - x(1), gapTb.time0);
                
                % Fit linear trend
                [ft, gof] = fit(x, y, 'poly1');
                
                % Store the data
                s = struct();
                s.win = gapWins(i,:);
                s.words = words(i:i+1);
                s.repInd = x;
                s.dur = y;
                s.durHat = ft(x);
                s.r2 = gof.adjrsquare;
                s.slope = ft.p1;
                s.intercept = ft.p2;
                gapData{i} = s;
            end
            
            gapData = gapData(~cellfun(@isempty, gapData));
        end
        
        function summary = CalculateSummaryStats(gapData)
            % Calculate summary statistics for gap data
            % 
            %   summary = CalculateSummaryStats(gapData)
            
            summary = struct();
            
            if isempty(gapData)
                return;
            end
            
            % Extract slopes and R² values
            slopes = cellfun(@(x) x.slope, gapData);
            r2Values = cellfun(@(x) x.r2, gapData);
            meanDurations = cellfun(@(x) mean(x.dur), gapData);
            
            summary.meanSlope = mean(slopes);
            summary.stdSlope = std(slopes);
            summary.meanR2 = mean(r2Values);
            summary.stdR2 = std(r2Values);
            summary.meanDuration = mean(meanDurations);
            summary.stdDuration = std(meanDurations);
            summary.nGaps = numel(gapData);
            summary.slopes = slopes;
            summary.r2Values = r2Values;
            summary.meanDurations = meanDurations;
        end
        
        function PlotSpeedMaps(recResults)
            % Generate speed map plots for all sentences in one figure
            % 
            %   PlotSpeedMaps(recResults, anaDir)
            
            % Filter out empty results
            recResults = recResults(~cellfun(@isempty, recResults));
            if isempty(recResults)
                return
            end
            
            % Create figure with subplots for each sentence
            nSen = numel(recResults);
            tl = tiledlayout(nSen*2, 1);
            tl.Padding = 'compact';
            tl.TileSpacing = "none";
            
            for i = 1 : nSen
                result = recResults{i};
                
                % Plot labels
                ax = nexttile;
                result.sen.PlotTiers();
                axis(ax, 'off');
                ax.XLim = result.tWin;
                MPlot.Axes(ax);
                
                % Plot speed map
                ax = nexttile;
                t = result.time;
                S = result.speed;
                imagesc(t-t(1), 1:size(S,2), S');
                colormap('jet');
                ax.XTick = [];
                % xlabel('Time (s)');
                ylabel('Repeats');
                MPlot.Axes(ax);

                c = colorbar(ax, 'Location', 'eastoutside');
                ylabel(c, 'Speed (X)');
            end
        end
        
        function PlotGapChanges(recResults)
            % Generate gap changes plots for all sentences in one figure
            % 
            %   PlotGapChanges(recResults)
            
            % Filter out empty results
            recResults = recResults(~cellfun(@isempty, recResults));
            if isempty(recResults)
                return
            end
            
            % Collect all gap data
            allGapData = {};
            for i = 1 : numel(recResults)
                result = recResults{i};
                if isempty(result.gapData)
                    continue
                end
                
                % Add sentence info to each gap
                for j = 1:numel(result.gapData)
                    gap = result.gapData{j};
                    gap.recId = result.recId;
                    gap.stimId = result.stimId;
                    allGapData{end+1} = gap;
                end
            end
            
            % Create figure with subplots for each gap
            tl = tiledlayout("flow");
            tl.Padding = 'compact';
            tl.TileSpacing = "tight";
            
            for i = 1 : numel(allGapData)
                ax = nexttile;
                gData = allGapData{i};
                
                plot(gData.repInd, gData.dur, 'k.', 'MarkerSize', 8);
                hold on;
                plot(gData.repInd, gData.durHat, 'r-', 'LineWidth', 1.5);
                
                ax.XTick = gData.repInd;
                ax.XLim = [0 gData.repInd(end)+1];
                ax.YLim(1) = 0;
                ax.YLim(2) = max([ax.YLim(2) 0.15]);
                title(sprintf('%s %s (R²=%.3f)', gData.words(1).GetParentLabel, gData.words(2).GetParentLabel, gData.r2));
                xlabel('Repeats');
                ylabel('Duration (s)');
                MPlot.Axes(ax);
            end
        end
        
        function PlotSummaryResults(allResults, recIds, stimIds)
            % Generate summary plots across all recordings and sentences
            % 
            %   PlotSummaryResults(allResults, recIds, stimIds)
            
            anaDir = LMV.Data.GetAnalysisDir("chunking");
            
            % Extract summary statistics
            nRecs = numel(recIds);
            nStims = numel(stimIds);
            
            meanSlopes = nan(nRecs, nStims);
            meanR2 = nan(nRecs, nStims);
            meanDurations = nan(nRecs, nStims);
            nGaps = nan(nRecs, nStims);
            
            for i = 1:nRecs
                for j = 1:nStims
                    if ~isempty(allResults{i, j}) && allResults{i, j}.success
                        meanSlopes(i, j) = allResults{i, j}.summary.meanSlope;
                        meanR2(i, j) = allResults{i, j}.summary.meanR2;
                        meanDurations(i, j) = allResults{i, j}.summary.meanDuration;
                        nGaps(i, j) = allResults{i, j}.summary.nGaps;
                    end
                end
            end
            
            % Plot summary heatmaps
            f = MPlot.Figure(2293); clf
            tl = tiledlayout(2, 2);
            tl.Padding = 'compact';
            
            % Mean slopes
            ax = nexttile;
            imagesc(meanSlopes);
            colorbar;
            title('Mean Slopes');
            xlabel('Sentence');
            ylabel('Recording');
            ax.XTick = 1:nStims;
            ax.XTickLabel = stimIds;
            ax.YTick = 1:nRecs;
            ax.YTickLabel = recIds;
            MPlot.Axes(ax);
            
            % Mean R²
            ax = nexttile;
            imagesc(meanR2);
            colorbar;
            title('Mean R²');
            xlabel('Sentence');
            ylabel('Recording');
            ax.XTick = 1:nStims;
            ax.XTickLabel = stimIds;
            ax.YTick = 1:nRecs;
            ax.YTickLabel = recIds;
            MPlot.Axes(ax);
            
            % Mean durations
            ax = nexttile;
            imagesc(meanDurations);
            colorbar;
            title('Mean Gap Durations');
            xlabel('Sentence');
            ylabel('Recording');
            ax.XTick = 1:nStims;
            ax.XTickLabel = stimIds;
            ax.YTick = 1:nRecs;
            ax.YTickLabel = recIds;
            MPlot.Axes(ax);
            
            % Number of gaps
            ax = nexttile;
            imagesc(nGaps);
            colorbar;
            title('Number of Gaps');
            xlabel('Sentence');
            ylabel('Recording');
            ax.XTick = 1:nStims;
            ax.XTickLabel = stimIds;
            ax.YTick = 1:nRecs;
            ax.YTickLabel = recIds;
            MPlot.Axes(ax);
            
            MPlot.Paperize(f, 1.5, 1.2);
            exportgraphics(f, fullfile(anaDir, 'summary_heatmaps.png'));
            close(f);
        end
        
    end
    
end 