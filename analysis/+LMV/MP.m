classdef MP
    
    methods(Static)
        % Labels
        function TrialEvents(ax, se)
            % A wrapper for NP.TaskBaseClass.PlotEventWindows
            
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > se.numEpochs || trialIdx < 1
                return;
            end
            
            if ~ismember('taskTime', se.tableNames)
                return;
            end
            
            colNames = {'atten', 'stim', 'delay', 'init', 'prod'};
            cc = LMV.Param.GetTaskPhaseColors(colNames);
            
            NP.TaskBaseClass.PlotEventWindows(ax, se, trialIdx, tWin, colNames, ...
                'YRange', [3 4], 'Style', 'block', 'Color', cc, 'Alpha', .5, 'Text', true, 'TextArgs', {'FontSize', 12});
            
            ax.XLabel.String = [];
            ax.XTickLabel = [];
            ax.YLabel.String = [];
            ax.YTick = [];
            axis(ax, 'off');
        end
        
        function TrialPhone(ax, se)
            % A wrapper for NP.TaskBaseClass.PlotTGEHier
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > se.numEpochs || trialIdx < 1
                return;
            end
            if ~ismember('taskTime', se.tableNames)
                return;
            end
            % cla(ax);
            NP.TaskBaseClass.PlotTGEHier(ax, se, trialIdx, tWin, [], 0.5, 1, 1, 'TierInd', 2:3, 'FontSize', 12);
            ax.Title.String = [];
            ax.XTickLabel = [];
            ax.XLabel.String = [];
            ax.XAxis.Visible = 'off';
            axis(ax, 'off');
        end
        
        function SentenceTitle(ax, se)
            % Add the text of sentence as axes title
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > se.numEpochs || trialIdx < 1
                return
            end
            
            % if ~ismember('taskTime', se.tableNames)
            %     return
            % end
            % tt = se.SliceEventTimes('taskTime', tWin, trialIdx, {'prod'}); % not using 'bleed' otherwise TGEs will get converted to double
            % if iscell(tt.prod)
            %     tge = tt.prod{1};
            % else
            %     tge = tt.prod;
            % end
            % if isnan(tge)
            %     return
            % end
            % s = char(tge.GetAllParentLabel);
            
            if ~ismember('taskValue', se.tableNames)
                return
            end
            s = se.GetTable("taskValue").stimText{trialIdx};
            
            s(1) = upper(s(1));
            ax.Title.String = ['"' s '"'];
            ax.Title.FontAngle = "italic";
        end
        
        % Features
        function PlotIntensity(ax, se)
            % A wrapper for NP.TaskBaseClass.PlotIntensity
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > se.numEpochs || trialIdx < 1
                return;
            end
            if ~ismember('taskTime', se.tableNames)
                return;
            end
            cla(ax);
            NP.TaskBaseClass.PlotIntensity(ax, se, trialIdx, tWin);
        end
        
        function TrialMelSpectrogram(ax, se)
            % A wrapper for NP.TaskBaseClass.PlotTrialMelSpectrogram
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > se.numEpochs || trialIdx < 1
                return;
            end
            if ~ismember('taskTime', se.tableNames)
                return;
            end
            cla(ax);
            NP.TaskBaseClass.PlotMelSpectrogram(ax, se, trialIdx, tWin);
            ax.XTickLabel = [];
            ax.XLabel.String = [];
        end
        
        function TrialFacialFeatures(ax, se)
            % A wrapper for NP.TaskBaseClass.PlotTrialFacialFeatures
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > se.numEpochs || trialIdx < 1
                return;
            end
            if ~ismember('taskTime', se.tableNames)
                return;
            end
            cla(ax);
            NP.TaskBaseClass.PlotFacialFeatures(ax, se, trialIdx, tWin);
            ax.XTickLabel = [];
            ax.YLim = [0 700];
        end
        
        function TrialSaccade(ax, se)
            % A wrapper for NP.TaskBaseClass.PlotTrialSaccade
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > se.numEpochs || trialIdx < 1
                return;
            end
            if ~ismember('taskTime', se.tableNames)
                return;
            end
            cla(ax);
            NP.MOCHA.PlotSaccade(ax, se, trialIdx, tWin);
            ax.XTickLabel = [];
            ax.XLabel.String = [];
%             ax.YLim = [0 700];
        end
        
        function Video(ax, se)
            
            trialIdx = ax.UserData.epoch;
            t = ax.UserData.time;
            if trialIdx > se.numEpochs || trialIdx < 1
                return
            end
            if ~isfield(se.userData, 'landmarkMeta')
                return
            end
            
            % Find video file
            expId = se.userData.experimentInfo.experimentId;
            vidFile = se.userData.landmarkMeta.file_path;
            if ~exist(vidFile, 'file')
                vidFile = fullfile(NP.Data.GetPreprocRoot, expId, 'landmarks', expId+".mp4");
            end
            if ~exist(vidFile, 'file')
                return
            end
            se.userData.landmarkMeta.file_path = vidFile;
            
            % Initialize image object (must be done before initializing overlaying objects)
            imgSize = [1000 1000]; % [height, width]
            if ~isfield(ax.UserData, 'img') || isempty(ax.UserData.img) || ~ishandle(ax.UserData.img)
                ax.UserData.img = imshow(zeros(imgSize), 'Border', 'tight', 'Parent', ax);
                ax.CLim = [0 255];
                colormap(ax, 'gray');
                axis(ax, 'image');
                ax.XAxis.Visible = 'off';
                ax.YAxis.Visible = 'off';
                ax.TickLength(1) = 0;
                hold(ax, 'on');
            end
            
            % Load the video file
            if ~isfield(ax.UserData, 'vidObj')
                ax.UserData.vidObj = VideoReader(vidFile);
                title(ax, ax.UserData.vidObj.Name, 'Interpreter', 'none');
                ax.TitleFontSizeMultiplier = 1.2;
            end
            vidObj = ax.UserData.vidObj;
            
            % Find video time
            tRef = se.GetReferenceTime();
            tGlb = t + tRef(trialIdx);
            tGlb = tGlb + 0.07; % why this big offset?
            
            % Read video frame
            if tGlb < 0 || tGlb > vidObj.Duration
                return
            end
            vidObj.CurrentTime = tGlb;
            fr = rgb2gray(readFrame(vidObj));
            
            % Transform frame
            if any(~isfield(ax.UserData, {'frTfObj', 'faceTfObj'}))
                switch char(expId)
                    case 'NP11_B4'
                        [fr, frTfObj] = Img23.Transform(fr, 'Translate', [-100 0], 'Rotate', -70, 'Crop', [1280 1280]);
                        [fr, faceTfObj] = Img23.Transform(fr, 'Translate', [25 -100], 'Crop', imgSize);
                    case 'NP30_B12'
                        [fr, frTfObj] = Img23.Transform(fr, 'Translate', [200 50], 'Rotate', 90, 'Crop', [1280 1280]);
                        [fr, faceTfObj] = Img23.Transform(fr, 'Translate', [-150 -150], 'Crop', imgSize);
                    case 'NP35_B2'
                        [fr, frTfObj] = Img23.Transform(fr, 'Translate', [-400 100], 'Rotate', -90, 'Crop', [1280 1280]);
                        [fr, faceTfObj] = Img23.Transform(fr, 'Translate', [25 -100], 'Crop', imgSize);
                end
                ax.UserData.frTfObj = frTfObj;
                ax.UserData.faceTfObj = faceTfObj;
            else
                fr = Img23.Transform(fr, ax.UserData.frTfObj, 'Crop', [1280 1280]);
                fr = Img23.Transform(fr, ax.UserData.faceTfObj, 'Crop', imgSize);
            end
            
            % Update plot
            ax.UserData.img.CData = fr;
        end
        
        function Landmark(ax, se)
            % 
            
            trialIdx = ax.UserData.epoch;
            t = ax.UserData.time;
            if trialIdx > se.numEpochs || trialIdx < 1
                return
            end
            if ~isfield(se.userData, 'landmarkMeta')
                return
            end
            
            if isfield(ax.UserData, 'faceTfObj')
                NP.LM.PlotLandmark(ax, se, trialIdx, t, ax.UserData.faceTfObj);
            else
                NP.LM.PlotLandmark(ax, se, trialIdx, t);
            end
        end
        
        % Response
        function EmbeddedResp(ax, se)
            % Make a scatter plot of unit responses at a given time point
            
            trialIdx = ax.UserData.epoch;
            t0 = ax.UserData.time;
            if trialIdx > se.numEpochs || trialIdx < 1
                return
            end
            if ~ismember('spikeRate', se.tableNames)
                return
            end
            
            cTb = se.userData.ksMeta.clusTb;
            cc = LMV.Param.GetTaskPhaseColors(cTb.tId1);
            
            [t, mm] = se.GetArray('spikeRate', trialIdx);
            [~, I] = min(abs(t - t0));
            r = mm(I,:);
            rc = cc .* r';
            
            if ~isfield(ax.UserData, 'dots') || ~isvalid(ax.UserData.dots)
                xy = cTb.embedCoords;
                h = scatter(ax, xy(:,1), xy(:,2), 12, cc, 'MarkerFaceColor', 'flat');
                ax.UserData.dots = h;
                ax.Color = 'k';
                axis(ax, 'square');
                axis(ax, 'tight');
                ax.XTick = [];
                ax.YTick = [];
                ax.XLim = ax.XLim + [-1 1]*0.2;
                ax.YLim = ax.YLim + [-1 1]*0.2;
                ax.Title.String = sprintf("Single-unit activity on functional embedding (n = %i)", height(cTb));
            else
                h = ax.UserData.dots;
                h.CData = rc;
            end
        end
        
        function LabelEmbeddedUnits(ax, senTb)
            % 
            
            senIdx = ax.UserData.epoch;
            if senIdx > height(senTb) || senIdx < 1
                return;
            end
            
            se = senTb.se(senIdx);
            if ~isfield(se.userData, 'expClusTb')
                return
            end
            T = se.userData.expClusTb;
            xy = T.embedCoords;
            r = repmat(0.1, size(xy(:,1)));
            s = string(1:numel(T.clusId));
            
            if isfield(ax.UserData, 'circles')
                delete(ax.UserData.circles);
                delete(ax.UserData.labels);
            end
            ax.UserData.circles = viscircles(ax, xy, r, 'Color', 'white', 'DrawBackgroundCircle', false, 'LineWidth', 1);
            % ax.UserData.labels = text(ax, xy(:,1), xy(:,2)+r(1)*1.2, s, 'Color', 'white', 'FontSize', 12, ...
            %     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            ax.UserData.labels = text(ax, xy(:,1)+r(1)*1.2, xy(:,2), s, 'Color', 'white', 'FontSize', 10, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        end
        
        function PethHeatmap(ax, ce)
            % Plot PETHs of individual units across conditions
            
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > ce.numEpochs || trialIdx < 1
                return
            end
            cla(ax);
            
            % Get response array
            uTb = ce.clusTb;
            [uTb, I] = sortrows(uTb, {'tId1', 'tId2', 'tId3', 'tR3'}, {'ascend', 'ascend', 'ascend', 'descend'});
            [t, mm] = ce.GetArray('resp', [], I+1, 'Normalization', 'max');
            y = 1 : size(mm,2);
            
            % Color responses
            cc = LMV.Param.GetTaskPhaseColors(uTb.tId1);
            cm = permute(cc, [1 3 2]) .* mm';
            
            % Plot heatmap
            imagesc(ax, t, y, cm);
            hold(ax, 'on');
            
            ax.Box = 'off';
            ax.Color = 'k';
            ax.XLim = tWin;
            ax.YDir = 'reverse';
            ax.YLabel.String = "# of units";
            ax.XLabel.String = "Aligned time (s)";
        end
        
        function RasterStack(ax, senTb)
            % A wrapper for NP.UnitPlot.RasterStack
            
            senIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if senIdx > height(senTb) || senIdx < 1
                return;
            end
            
            se = senTb.se(senIdx);
            if ~isfield(se.userData, 'expClusTb')
                return
            end
            T = se.userData.expClusTb;
            clusInd = NP.Unit.ClusId2Ind(T.clusId, se);
            
            cla(ax);
            cc = LMV.Param.GetTaskPhaseColors(string(T.tId1));
            
            NP.UnitPlot.RasterStack(ax, se, [], tWin, clusInd, ...
                'Color', cc, 'LineWidth', 1.5, ...
                'PETH', false, 'PETHArgs', {'Color', cc}, ...
                'PeakSpikeRates', T.peakSpkRate);
            
            nUnits = numel(clusInd);
            xx = repmat(ax.XLim', 1, nUnits-1);
            yy = repmat((2:nUnits)-0.5, 2, 1);
            plot(ax, xx, yy, 'Color', [1 1 1]*0.5, 'LineWidth', 1);
            
            ax.YTickLabel = "u"+(1:nUnits);
            ax.Color = 'k';
        end
        
        % Linker
        function LinkerTrialEvents(ax, se)
            % A wrapper for NP.TaskBaseClass.PlotEventWindows
            
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > se.numEpochs || trialIdx < 1
                return;
            end
            
            if ~ismember('taskTime', se.tableNames)
                return;
            end
            
            colNames = {'stim', 'prod'};
            % cc = LMV.Param.GetTaskPhaseColors(se.userData.phaseName);
            cc = LMV.Param.GetTaskPhaseColors(colNames);
            
            NP.TaskBaseClass.PlotEventWindows(ax, se, trialIdx, tWin, colNames, ...
                'YRange', [2 3], 'Style', 'block', 'Color', cc, 'Alpha', .5, 'Text', false, 'TextArgs', {'FontSize', 12});
            
            ax.XLabel.String = [];
            ax.XTickLabel = [];
            ax.YLabel.String = [];
            ax.YTick = [];
            axis(ax, 'off');
        end
        
        function LinkerTrialPhone(ax, se)
            % A wrapper for NP.TaskBaseClass.PlotTGEHier
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > se.numEpochs || trialIdx < 1
                return;
            end
            if ~ismember('taskTime', se.tableNames)
                return;
            end
            NP.TaskBaseClass.PlotTGEHier(ax, se, trialIdx, tWin, [], 0.5, 1, 1, 'TierInd', 2, 'FontSize', 12);
            ax.Title.String = [];
            ax.XTickLabel = [];
            ax.XLabel.String = [];
            ax.XAxis.Visible = 'off';
            axis(ax, 'off');
        end
        
        function LinkerResponseStack(ax, ce)
            % A wrapper for MPlot.PlotHistStack and MPlot.PlotRasterStack
            
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > ce.numEpochs || trialIdx < 1
                return
            end
            
            % Configure by task phase
            % if ce.userData.phaseName == "stim"
            %     cla(ax);
            % end
            expClusTb = ce.userData.expClusTb;
            N = height(expClusTb);
            c = LMV.Param.GetTaskPhaseColors(repelem(ce.userData.phaseName, N));
            
            % Plot spike raster
            st = expClusTb.spikeTime(:,trialIdx)';
            y0 = (1:N)+0.25;
            MPlot.PlotRasterStack(st, y0, 'HeightScale', 0.4, 'Color', c, 'Parent', ax);
            
            % Get averaged responses
            [~, I] = MMath.SortLike(ce.clusTb.clusId, expClusTb.clusId);
            [t, mm] = ce.GetArray('resp', trialIdx, I+1);
            [~, sem] = ce.GetArray('sem', trialIdx, I+1);
            
            % Plot spike rates
            MPlot.PlotHistStack(t, mm+0.5, sem, 'Scaling', 0.4, 'Color', c, 'Parent', ax);
            ax.XLim = tWin;
            ax.YLabel.String = "Neuron #";
            ax.XLabel.String = "Aligned time (s)";
        end
        
    end
end

