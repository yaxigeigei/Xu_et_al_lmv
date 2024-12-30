classdef Mov
    % Mirro plays
    
    methods(Static)
        function SetupPlotter(mp, suffixName)
            
            import LMV.Linker.Mov
            
            mp.RemovePlot();
            
            figNum = sum(suffixName{:});
            rowDist = [1 2 3];
            colDist = 1;
            
            spInd = MPlot.FindSubplotInd(rowDist, colDist, 1, 1);
            spStr = MPlotter.SubplotArgs2Str(spInd{:});
            mp.AddPlot(figNum, spStr, @cla, '', 'epoch');
            mp.AddPlot(figNum, spStr, @(a,b) Mov.TrialPhone(a, b), 'ce'+suffixName, 'epoch');
            mp.AddPlot(figNum, spStr, @(a,b) Mov.TrialEvents(a, b), 'ce'+suffixName, 'epoch');
            mp.AddPlot(figNum, spStr, @(a,b) Mov.PhaseTitle(a, b), 'ce'+suffixName, 'epoch');
            mp.AddPlot(figNum, spStr);
            
            spInd = MPlot.FindSubplotInd(rowDist, colDist, 2, 1);
            spStr = MPlotter.SubplotArgs2Str(spInd{:});
            mp.AddPlot(figNum, spStr, @Mov.TrialMelSpectrogram, 'seTemp'+suffixName, 'epoch');
            mp.AddPlot(figNum, spStr);
            
            spInd = MPlot.FindSubplotInd(rowDist, colDist, 3, 1);
            spStr = MPlotter.SubplotArgs2Str(spInd{:});
            mp.AddPlot(figNum, spStr, @Mov.RasterStack, 'senTb'+suffixName, 'epoch');
            mp.AddPlot(figNum, spStr);
            
            mp.RefreshAll();
            f = figure(mp.plotTable.figureObj{1});
            % f.Position(1:2) = [0 50];
            f.Position(3:4) = [800 400];
            MPlot.Paperize(f, 'FontSize', 8);
        end
        
        function frames = MakeVideos(saveDir, mp, se, epInd)
            % 
            % 
            %   frames = MakeVideos(saveDir, mp, se, epInd)
            % 
            
            % Get variables
            f = figure(mp.plotTable.figureObj{1});
            recId = NP.SE.GetID(se);
            phaseName = se.userData.phaseName;
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            expClusId = se.userData.expClusTb.clusId;
            
            if ~exist("epInd", "var")
                epInd = 1 : se.numEpochs;
            end
            
            for i = epInd(:)'
                % Get time window
                tWin = [tt.speechMatchOn(i) tt.speechMatchOff(i)] + [-1 1]*0.3;
                
                % Set epoch and time
                mp.timeLimits = tWin;
                mp.epoch = i;
                mp.time = tWin(1);
                
                MPlot.Paperize(f, 'FontSize', 8);
                
                % Make movie
                vidFile = fullfile(saveDir, sprintf("u%i_%s_sent%i_trial%i.mp4", expClusId, phaseName, i, tv.trialNum(i)));
                fps = 60;
                frames = mp.MakeVideo(f, 1/fps, 'FrameRate', fps, 'FilePath', vidFile);
                
                % return
            end
        end
        
        function MakeSpeechAudio(saveDir, seTemp, epInd)
            % 
            % 
            %   MakeSpeechAudio(saveDir, seTemp, epInd)
            % 
            
            % Get variables
            recId = NP.SE.GetID(seTemp);
            phaseName = seTemp.userData.phaseName;
            [tt, tv] = seTemp.GetTable('taskTime', 'taskValue');
            
            for i = epInd(:)'
                % Get time window
                tWin = [tt.speechMatchOn(i) tt.speechMatchOff(i)] + [-1 1]*0.3;
                
                % Save speech audio
                if phaseName == "stim"
                    chan = "speaker1";
                else
                    chan = "mic";
                end
                micFile = sprintf("%s_%s_sent%i_trial%i.wav", recId, phaseName, i, tv.trialNum(i));
                NP.Audio.WriteTrial(seTemp, i, tWin, 'Channel', chan, 'FolderPath', saveDir, 'FileName', micFile);
                
                % return
            end
        end
        
        function MakeSpikeAudio(saveDir, seTemp, epInd, se, sr)
            % 
            % 
            %   MakeSpikeAudio(saveDir, seTemp, epInd, sr)
            % 
            
            % Get variables
            clusId = seTemp.userData.expClusTb.clusId;
            clusIdShort = clusId - NP.Unit.GetBaseClusId(clusId);
            phaseName = seTemp.userData.phaseName;
            [tt, tv] = seTemp.GetTable('taskTime', 'taskValue');
            
            for i = epInd(:)'
                % Get time window
                tWin = [tt.speechMatchOn(i) tt.speechMatchOff(i)] + [-1 1]*0.3;
                
                % Convert epoch time to recording time
                isTrial = se.GetTable('taskValue').trialNum == tv.trialNum(i);
                rt = se.GetReferenceTime;
                if phaseName == "stim"
                    tOn = se.GetTable("taskTime").stimOn(isTrial);
                else
                    tOn = se.GetTable("taskTime").prodOn{isTrial}(1);
                end
                tWin = tWin + tOn + rt(isTrial);
                
                % Generate spike audio
                spkFile = fullfile(saveDir, sprintf("u%i_%s_sent%i_trial%i.wav", clusId, phaseName, i, tv.trialNum(i)));
                sr.WriteSpikeAudio(spkFile, clusIdShort, tWin);
            end
        end
        
        function SentenceTitle(ax, se)
            % Add the text of sentence as axes title
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > se.numEpochs || trialIdx < 1
                return
            end
            
            if ~ismember('taskValue', se.tableNames)
                return
            end
            s = se.GetTable("taskValue").stimText{trialIdx};
            
            s(1) = upper(s(1));
            ax.Title.String = ['"' s '"'];
            ax.Title.FontAngle = "italic";
        end
        
        function PhaseTitle(ax, se)
            % Add the text of sentence as axes title
            trialIdx = ax.UserData.epoch;
            if trialIdx > se.numEpochs || trialIdx < 1
                return
            end
            if ~isfield(se.userData, 'phaseName')
                return
            end
            if se.userData.phaseName == "stim"
                lb = "Listen";
            else
                lb = "Speak";
            end
            ax.Title.String = lb;
        end
        
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
            
            colNames = {'stim', 'prod'};
            % cc = LMV.Param.GetTaskPhaseColors(se.userData.phaseName);
            cc = LMV.Param.GetTaskPhaseColors(colNames);
            
            NP.TaskBaseClass.PlotEventWindows(ax, se, trialIdx, tWin, colNames, ...
                'YRange', [2 3], 'Style', 'block', 'Color', cc, 'Alpha', .5, 'Text', false, 'TextArgs', {'FontSize', 12});
            
            ax.YLim = [0 3];
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
            NP.TaskBaseClass.PlotTGEHier(ax, se, trialIdx, tWin, [], 0.5, 1, 1, 'TierInd', 2, 'FontSize', 12);
            ax.Title.String = [];
            ax.XTickLabel = [];
            ax.XLabel.String = [];
            ax.XAxis.Visible = 'off';
            axis(ax, 'off');
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
            
            if isfield(se.userData, "phaseName") && se.userData.phaseName == "stim"
                chan = "speaker1";
            else
                chan = "mic";
            end
            
            NP.TaskBaseClass.PlotMelSpectrogram(ax, se, trialIdx, tWin, 'ChannelName', chan, 'WaveformHeight', 0);
            
            ax.XTick = [];
            ax.XTickLabel = [];
            ax.XLabel.String = [];
            ax.YLabel.String = "Spectrogram";
            ax.YTick = [];
        end
        
        function RasterStack(ax, senTb)
            % A wrapper for NP.UnitPlot.RasterStack
            
            senIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if senIdx > height(senTb) || senIdx < 1
                return;
            end
            
            se = senTb.se(senIdx);
            cTb = se.userData.expClusTb;
            clusInd = NP.Unit.ClusId2Ind(cTb.clusId, se);
            ce = senTb.ce(senIdx);
            
            cla(ax);
            NP.UnitPlot.RasterStack(ax, se, [], tWin, clusInd, 'LineWidth', 1.5, ...
                'PETH', true, ...
                'PeakSpikeRates', ce.clusTb.peakSpkRate, ...
                'PETHArgs', {'Color', LMV.Param.GetTaskPhaseColors(senTb.phase(senIdx)), 'Scaling', 0.8, 'LineWidth', 1.5});
            
            ax.XTickLabel = [];
            ax.XLabel.String = "Time";
            ax.YTick = [];
            ax.YTickLabel = [];
            ax.YLabel.String = "Trials";
            ax.XGrid = "off";
            % ax.YTickLabel = "u"+(1:nUnits);
            % ax.Color = 'k';
        end
        
        function ResponseStack(ax, ce)
            % A wrapper for MPlot.PlotHistStack and MPlot.PlotRasterStack
            
            trialIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if trialIdx > ce.numEpochs || trialIdx < 1
                return
            end
            
            % Configure by task phase
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
        
        % Bridge play
        function RealtimeRaster(ax, se)
            % A wrapper for NP.UnitPlot.RasterStack
            
            epIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if epIdx > se.numEpochs || epIdx < 1
                return
            end
            
            % Find unit
            cTb = se.userData.expClusTb;
            clusInd = NP.Unit.ClusId2Ind(cTb.clusId, se);
            
            % Slice out spike times data
            st = se.SliceEventTimes('spikeTime', tWin, [], clusInd, 'Fill', 'none');
            st = st.(1);
            
            % Plot raster
            % cla(ax);
            [t, y] = NP.UnitPlot.ConvertEventTimesForRasters(st);
            hh = MPlot.PlotPointAsLine(t, y, .6, 'Color', [0 0 0], 'LineWidth', 0.5, 'Parent', ax);
            hh.UserData.t = hh.XData;
            hh.UserData.y = hh.YData;
            
            ax.UserData.hh = hh;
            ax.XLim = tWin;
            % ax.YLim = [0 se.numEpochs+1];
            ax.YLim = [0 epIdx+1];
            ax.XTickLabel = [];
            ax.XLabel.String = "Time";
            ax.YTick = [];
            ax.YLabel.String = "Trials";
        end
        
        function UpdateRaster(ax, se)
            % Update the plot created by RealtimeRaster
            
            epIdx = ax.UserData.epoch;
            tEnd = ax.UserData.time;
            if epIdx > se.numEpochs || epIdx < 1
                return
            end
            if ~isfield(ax.UserData, "hh") || ~isvalid(ax.UserData.hh)
                return
            end
            
            hh = ax.UserData.hh;
            t = hh.UserData.t;
            y = hh.UserData.y;
            y(y>epIdx+0.5 | (y>epIdx-0.5 & t>tEnd)) = NaN;
            hh.YData = y;
        end
        
        function DelayBoundary(ax, se)
            % A wrapper of NP.TaskBaseClass.PlotEventBoundary
            
            epIdx = ax.UserData.epoch;
            tWin = ax.UserData.timeLimits;
            if epIdx > se.numEpochs || epIdx < 1
                return
            end
            
            % evtNames = ["stimMatchOff", "cue3On", "cue3Off", "prodMatchOn"];
            evtNames = ["stimMatchOff", "prodMatchOn"];
            NP.TaskBaseClass.PlotEventBoundary(ax, se, 1:epIdx, tWin, evtNames);
            
            ax.YDir = "normal";
        end
        
        function frames = MakeBridgePlayVideos(saveDir, mp, se, epInd)
            % 
            % 
            %   frames = MakeVideos(saveDir, mp, se, epInd)
            % 
            
            % Get variables
            f = figure(mp.plotTable.figureObj{1});
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            expClusId = se.userData.expClusTb.clusId;
            
            if ~exist("epInd", "var")
                epInd = 1 : se.numEpochs;
            end
            
            for i = epInd(:)'
                % Set epoch and time
                mp.epoch = i;
                mp.time = 0;
                
                MPlot.Paperize(f, 'FontSize', 8);
                
                % Make movie
                vidFile = fullfile(saveDir, sprintf("u%i_%s_rep%i_trial%i.mp4", expClusId, tv.stimId(i), i, tv.trialNum(i)));
                fps = 60;
                frames = mp.MakeVideo(f, 1/fps, 'FrameRate', fps, 'FilePath', vidFile);
                
                % return
            end
        end
        
        function MakeBridgePlaySpeechAudio(saveDir, se, epInd)
            % 
            % 
            %   MakeBridgePlaySpeechAudio(saveDir, se, epInd)
            % 
            
            % Get variables
            recId = NP.SE.GetID(se);
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            
            if ~exist("epInd", "var")
                epInd = 1 : se.numEpochs;
            end
            
            for i = epInd(:)'
                % Get time window
                tWin = [tt.trialOn(i) tt.prodMatchOff(i)] + [0 1]*0.3;
                
                % Save speech audio
                micFile = sprintf("%s_%s_rep%i_trial%i.wav", recId, tv.stimId(i), i, tv.trialNum(i));
                NP.Audio.WriteTrial(se, i, tWin, 'Channel', "mic", 'FolderPath', saveDir, 'FileName', micFile);
                
                % return
            end
        end
        
        function MakeBridgePlaySpikeAudio(saveDir, seSen, se, sr, epInd)
            % 
            % 
            %   MakeBridgePlaySpikeAudio(saveDir, seSen, se, sr, epInd)
            % 
            
            % Get variables
            clusId = seSen.userData.expClusTb.clusId;
            clusIdShort = clusId - NP.Unit.GetBaseClusId(clusId);
            [tt, tv] = seSen.GetTable('taskTime', 'taskValue');
            
            if ~exist("epInd", "var")
                epInd = 1 : seSen.numEpochs;
            end
            
            for i = epInd(:)'
                % Get time window
                tWin = [tt.trialOn(i) tt.prodMatchOff(i)] + [0 1]*0.3;
                
                % Convert epoch time to recording time
                isTrial = se.GetTable('taskValue').trialNum == tv.trialNum(i);
                rt = se.GetReferenceTime;
                tWin = tWin + rt(isTrial);
                
                % Generate spike audio
                spkFile = fullfile(saveDir, sprintf("u%i_%s_rep%i_trial%i.wav", clusId, tv.stimId(i), i, tv.trialNum(i)));
                sr.WriteSpikeAudio(spkFile, clusIdShort, tWin);
            end
        end
        
        
        
    end
    
end

