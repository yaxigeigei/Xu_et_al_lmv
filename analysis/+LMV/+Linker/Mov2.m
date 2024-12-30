classdef Mov2
    % Mirro plays
    
    methods(Static)
        function frames = MakeVideos(saveDir, mp, seTemp, epInd)
            % 
            % 
            %   frames = MakeVideos(saveDir, mp, seTemp, epInd)
            % 
            
            % Get variables
            f = figure(mp.plotTable.figureObj{1});
            recId = NP.SE.GetID(seTemp);
            phaseName = seTemp.userData.phaseName;
            [tt, tv] = seTemp.GetTable('taskTime', 'taskValue');
            expClusId = seTemp.userData.expClusTb.clusId;
            
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
        
        function MakeSpikeAudio(saveDir, seTemp, epInd ,se, sr)
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
        
    end
    
end

