classdef TIMIT < NP.TaskBaseClass
    
    methods(Static)
        function [tv, tt] = AddTaskValueTable(se)
            % Make and add taskValue table to se
            % 
            %   tv = AddTaskValueTable(se)
            % 
            
            % Continue with existing table or create new
            if ismember('taskValue', se.tableNames)
                tv = se.GetTable('taskValue');
            else
                tv = table;
            end
            tt = se.GetTable('taskTime');
            tv.trialNum = (1:height(tt))';
            tv.recId(:) = string(se.userData.expInfo.recId);
            se.SetTable('taskValue', tv, 'eventValues');
            
            % Set stim ID in taskValue
            NP.TaskBaseClass.MigrateStimId(se, '^.{5}_si.{3,}');
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            
            % Check required columns in taskTime table
            requiredCols = {'stim'};
            assert(all(ismember(requiredCols, tt.Properties.VariableNames)), ...
                "TIMIT task requires the 'stim' column in the taskTime table.");
            
            numTrials = 0;
            for i = 1 : height(tt)
                % Get task name
                eStim = tt.stim(i);
                if iscell(eStim)
                    eStim = eStim{1};
                end
                if isnan(eStim)
                    continue
                end
                taskName = lower(eStim(1).GetVfield('task'));
                if taskName ~= "timit"
                    continue
                end
                numTrials = numTrials + 1;
                
                % Get stim and prod texts
                stimText = eStim(1).GetParentLabel;
                
                tv.taskName(i) = taskName;
                tv.stimText(i) = stimText;
            end
            
            if ~numTrials
                fprintf("This se does not have any TIMIT trial\n");
                return
            end
            
            % Make unique numeric stim IDs
            m = tv.taskName == "timit";
            [~, ~, tv.stimNumId(m)] = unique(tv.stimText(m), 'stable');
            
            % Add taskValues table to se
            se.SetTable('taskValue', tv);
        end
        
        function seTb = SplitBySentence(se)
            % Split se by unique sentences
            
            % Make a lite se with only the necessary data
            seLite = se.Duplicate({'taskTime', 'taskValue', 'spikeRate'}, false);
            seLite.userData.experimentInfo = se.userData.experimentInfo;
            
            % Remove trials that are not TIMIT
            tv = seLite.GetTable('taskValue');
            seLite.RemoveEpochs(lower(tv.taskName) ~= "timit");
            
            % Split se by unique sentences
            ops = struct;
            ops.conditionVars = {'stimNumId'};
            seTb = NP.SE.SplitConditions(seLite, ops);
            
            % Add group info
            tv = seLite.GetTable('taskValue');
            seTb.stimText = unique(tv.stimText, 'stable');
            seTb.trialNum = arrayfun(@(x) x.GetColumn('taskValue', 'trialNum'), seTb.se, 'Uni', false);
        end
        
        function uTb = ComputePETH(seTb)
            % Compute PETHs from trials of each seTb row
            
            % Add PETH parameters to each se
            for i = 1 : height(seTb)
                tt = seTb.se(i).GetTable('taskTime');
                ops = NP.Param.Resample();
                ops.rsWin = [tt.stimOn(1) tt.stimOff(1)] + [-1 1]*0.5;
                ops.rsBinSize = 0.005;
                seTb.se(i).userData.ops = ops;
            end
            
            % Compute PETHs
            m = seTb.numTrial > 2;
            uTb = NP.PETH.ComputePETH(seTb.se(m));
        end
        
        function PlotTrialNP30B12(se, seTb, uTb, isSave)
            % 
            
            % Select sentences
            m = seTb.numTrial >= 5; % exclude sentences with too few repeats
            seTb = seTb(m,:);

            % Find the indices of the first trial of each group in se
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            trInd = cellfun(@(x) find(tv.trialNum == x(1), 1), seTb.trialNum);
            
            % Time window to plot
            tWin = [tt.stimOn(trInd) tt.stimOff(trInd)] + [-1 1]*0.5;
            
            % Units to plot
            uId = [1067 788 663 927 1244 314 321 1085 1276 1403 1279 562]; % from Sean
            uInd = NP.Unit.ClusId2Ind(uId, se);
            uInd = sort(uInd); % descending distance to tip
            uInd4Plot = repmat(uInd(:), [1 height(seTb)]); % replicate indices for all groups
            
            % Determine the number of figures and sentences per figure
            nSen = numel(trInd);
            nSenPerFig = 5;
            nFig = ceil(nSen / nSenPerFig);
            
            % Make plots
            f = gcf;
            f.WindowState = 'maximized';
            for k = 1 : nFig
                I = (k-1)*nSenPerFig+1 : k*nSenPerFig;
                
                clf;
                NP.TaskBaseClass.PlotTrial(se, tWin(I,:), trInd(I), uInd4Plot(:,I), uTb.peakSpkRate(uInd4Plot(:,I)), ...
                    'FeatureList', {'phone', 'mel'});
                
                tightfig;
                %             f.WindowState = 'maximized';
                %             drawnow();
                
                % Save figure
                if exist('isSave', 'var') && isSave
                    figDir = fullfile(NP.Data.GetAnalysisRoot, se.userData.experimentInfo.experimentId);
                    if ~exist(figDir, 'dir')
                        mkdir(figDir);
                    end
                    figName = "timit_" + "s" + I(1) + "-" + I(end) + "_"+ strjoin("u" + uId, "_");
                    saveas(f, fullfile(figDir, figName), 'png');
                end
            end
        end
        
        function PlotGroups(seTb, uTb, varargin)
            % Present each seTb row in a column of plots. Top are feature plots of the first trial. 
            % Bottom are raster and PETH overlay of individual units across trials in the group.
            % 
            %   PlotGroups(se, uTb)
            %   PlotGroups(se, uTb, unitInd)
            %   PlotGroups(se, uTb, unitInd, pkSpkRate)
            %   PlotGroups(..., 'PlotList', {'phone', 'acous', 'mel', 'landmark'})
            % 
            % Inputs
            %   se              A MSessionExplorer object.
            %   tWin            A n-by-2 numeric array of time windows. n is the number of trials.
            %                   If n == 1, this window will be applied to all trials.
            %   trialInd        A n-element vector of trial indices for the n trials. Each trial is plotted
            %                   in a column. One can index the same trial multiple times, e.g. [1 1 2 2].
            %   unitInd         The positional indices of units in data structures. Default is [], 
            %                   where no raster and PETH will be plotted. When it is a vector, these same 
            %                   units are plotted for each trial (column). When it is a m-by-n matrix, where 
            %                   n is the number of trials, different columns can have different units.
            %   pkSpkRate       A vector of the same size as unitInd that overwrite the peak spike 
            %                   rates the PETHs are normalized to.
            % 
            
            featList = {'phone', 'acous', 'mel', 'landmark', 'event'};
            
            % Parse user inputs
            p = inputParser();
            p.addOptional('unitInd', [], @isnumeric);
            p.addOptional('pkSpkRate', [], @isnumeric);
            p.addParameter('FeatureList', {'phone', 'acous', 'mel'}, @(x) all(ismember(x, featList)));
            p.parse(varargin{:});
            unitInd = p.Results.unitInd;
            pkSpkRate = p.Results.pkSpkRate;
            featList = p.Results.FeatureList;
            
            if isvector(unitInd)
                unitInd = repmat(unitInd(:), [1 numel(trialInd)]);
                pkSpkRate = repmat(pkSpkRate(:), [1 numel(trialInd)]);
            end
            
            % Plot parameters
            numRows = numel(featList) + round(size(unitInd,1) * 0.5);
            numCols = numel(trialInd);
            k = 0;
            
            % Features
            for i = 1 : numel(featList)
                for j = 1 : height(seTb)
                    
                    se = seTb.se(j);
%                     [tt, tv] = se.GetTable('taskTime', 'taskValue');
                    tr = 1;
                    w = se.userData.ops.rsWin;
                    
                    k = k + 1;
                    ax = MPlot.Axes(numRows, numCols, k); cla
                    
                    switch featList{i}
                        case 'phone'
                            NP.TaskBaseClass.PlotTrialPhone(ax, se, tr, w);
                        case 'acous'
                            NP.TaskBaseClass.PlotTrialAcous(ax, se, tr, w);
                        case 'mel'
                            NP.TaskBaseClass.PlotTrialMelSpectrogram(ax, se, tr, w);
                        case 'landmark'
                            NP.TaskBaseClass.PlotTrialFacialFeatures(ax, se, tr, w);
                        case 'event'
                            NP.TaskBaseClass.PlotTrialEvents(ax, se, tr, w);
                    end
                    
                    % Omit x-label if it's not the last row
                    if i ~= numel(featList) || ~isempty(unitInd)
                        ax.XLabel.String = [];
                    end
                end
            end
            
            % Rasters and PETHs
            if isempty(unitInd)
                return
            end
            for j = 1 : height(seTb)
                rowInd = numel(featList)+1 : numRows;
                colInd = ones(size(rowInd)) * j;
                k = sub2ind([numCols numRows], colInd, rowInd);
                ax = MPlot.Axes(numRows, numCols, k); cla
                
                se = seTb.se(j);
%                 w = tWin(j,:);
                w = se.userData.ops.rsWin;
                NP.TIMIT.PlotRasterPETH(ax, se, w, uTb(unitInd(:,j),:));
            end
            
        end
        
    end
end

