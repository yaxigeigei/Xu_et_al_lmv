classdef NBD < NP.TaskBaseClass
    
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
                "NBD task requires the 'stim' column in the taskTime table.");
            
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
                if taskName ~= "nbd"
                    continue
                end
                numTrials = numTrials + 1;
                
                % Get stim and prod texts
                stimText = eStim(1).GetParentLabel;
                
                tv.taskName(i) = taskName;
                tv.stimText(i) = stimText;
            end
            
            if ~numTrials
                fprintf("This se does not have any nbd trial\n");
                return
            end
            
            % Make unique numeric stim IDs
            m = tv.taskName == "nbd";
            [~, ~, tv.stimNumId(m)] = unique(tv.stimText(m), 'stable');
            
            % Add taskValues table to se
            se.SetTable('taskValue', tv);
        end

        function [uId, uIdUni] = GetSelectedClusId(recId)
            % Return handpicked cluster IDs for a given recording
            % 
            %   [uId, uIdUni] = GetSelectedClusId(recId)
            % 
            switch recId
                case 'NP102_B2', uId = [1020200143 1020200118 1020200113 1020200242 1020200273 1020200272 1020200137];
                case 'NP108_B1', uId = [1080100202 1080100179 1080100000 1080100315 1080100312 1080100305 1080100365 1080100200 1080100101];
                otherwise, uId = [];
            end
            uIdUni = NP.Unit.GetBaseClusId(recId) + uId;
        end

        function PlotOverview(se, varargin)
            % Plot session heatmap
            
            % Determine units to plot
            ud = se.userData;
            ut = NP.Unit.GetClusTb(se);
            if numel(varargin) == 1
                % By unit ID
                uId = varargin{1};
                uInd = NP.Unit.ClusId2Ind(uId, se);
                uInd = sort(uInd);
                p = 0;
                upp = 10;
            elseif numel(varargin) == 2
                % By page index and # of units per page
                [p, upp] = varargin{:};
                I = (1:upp)+upp*(p-1);
                I(I > height(ut)) = []; % remove out of range indices
                uInd = I;
            end
            uInd = uInd(:)';
            uId = NP.Unit.Ind2ClusId(uInd, se);
            
            % Configure plots
            nr = 2; % number of rows
            ppr = ceil(upp/nr) + 1; % plot per row = #unit / 3 rows + one label plot
            panels1 = repmat({'rate_map'}, [nr ppr]);
            panels1args = repmat({{}}, [nr ppr]);
            
            panels2 = repmat({'event_boundary'}, size(panels1));
            panels2args = repmat({{{'lastWordOn', 'lastWordOff'}, @(x) brewermap(x, 'Dark2')}}, size(panels1));
            
            panels3 = repmat({'block_boundary'}, size(panels1));
            % panels3args = repmat({{{'arithmetic_secondOperand'},'Label', false}}, size(panels1));
            panels3args = repmat({{{'surprisal_label'},'Label', false}}, size(panels1));
            
            k = 1;
            for i = 1 : size(panels1,1)
                % first column for labels
                panels1{i,1} = ''; % not plotting spikes in the first column
                panels2{i,1} = ''; % not plotting event_boundary in the first column
                panels3args{i,1}{3} = true; % show labels in the fir
                
                % following columns for rasters
                for j = 2 : size(panels1,2)
                    if k > numel(uInd)
                        panels1{i,j} = '';
                        panels2{i,j} = '';
                        panels3{i,j} = '';
                        continue
                    end
                    panels1args{i,j} = {uInd(k)};
                    k = k + 1;
                end
            end
            
            tWin = [-0.5, 2];
            
            % Plotting
            f = MPlot.Figure(38050301); clf
            f.WindowState = 'maximized';
            NP.PTBRead.PlotSE(se, panels1, 'PanelArgs', panels1args, 'TimeWindow', tWin);
            NP.PTBRead.PlotSE(se, panels2, 'PanelArgs', panels2args, 'TimeWindow', tWin);
            NP.PTBRead.PlotSE(se, panels3, 'PanelArgs', panels3args, 'TimeWindow', tWin);
            
            % Save figure
            recId = string(NP.SE.GetID(ud));
            region = NP.SE.GetRegion(ud(1));
            figDir = fullfile(NP.Data.GetAnalysisRoot, "se_"+region, recId);
            if ~exist(figDir, 'dir')
                mkdir(figDir);
            end
            % figName = strjoin([recId, "semsr", "page"+p, "u"+uId], "_");
            figName = strjoin([recId, "nbd", "page"+p], "_"); %filename is often too long with all the unitIds included
            if ~exist(figDir, 'dir')
                mkdir(figDir);
            end
            exportgraphics(f, fullfile(figDir, figName+".png"));
        end

        function PlotHistograms(se, varargin)
             % Plot session heatmap
            
            % Determine units to plot
            ud = se.userData;
            ut = NP.Unit.GetClusTb(se);
            if numel(varargin) == 1
                % By unit ID
                uId = varargin{1};
                uInd = NP.Unit.ClusId2Ind(uId, se);
                uInd = sort(uInd);
                p = 0;
                upp = 10;
            elseif numel(varargin) == 2
                % By page index and # of units per page
                [p, upp] = varargin{:};
                I = (1:upp)+upp*(p-1);
                I(I > height(ut)) = []; % remove out of range indices
                uInd = I;
            end
            uInd = uInd(:)';
            uId = NP.Unit.Ind2ClusId(uInd, se);
            
            % Configure plots
            nr = 2; % number of rows
            ppr = ceil(upp/nr) + 1; % plot per row = #unit / 3 rows + one label plot
            panels1 = repmat({'rate_map'}, [nr ppr]);
            panels1args = repmat({{}}, [nr ppr]);
            
            panels2 = repmat({'event_boundary'}, size(panels1));
            panels2args = repmat({{{'lastWordOn', 'lastWordOff'}, @(x) brewermap(x, 'Dark2')}}, size(panels1));
            
            panels3 = repmat({'block_boundary'}, size(panels1));
            % panels3args = repmat({{{'arithmetic_secondOperand'},'Label', false}}, size(panels1));
            panels3args = repmat({{{'surprisal_label'},'Label', false}}, size(panels1));
            
            k = 1;
            for i = 1 : size(panels1,1)
                % first column for labels
                panels1{i,1} = ''; % not plotting spikes in the first column
                panels2{i,1} = ''; % not plotting event_boundary in the first column
                panels3args{i,1}{3} = true; % show labels in the fir
                
                % following columns for rasters
                for j = 2 : size(panels1,2)
                    if k > numel(uInd)
                        panels1{i,j} = '';
                        panels2{i,j} = '';
                        panels3{i,j} = '';
                        continue
                    end
                    panels1args{i,j} = {uInd(k)};
                    k = k + 1;
                end
            end
            
            tWin = [-0.5, 2];
            
            % Plotting
            f = MPlot.Figure(38050301); clf
            f.WindowState = 'maximized';
            NP.PTBRead.PlotSE(se, panels1, 'PanelArgs', panels1args, 'TimeWindow', tWin);
            NP.PTBRead.PlotSE(se, panels2, 'PanelArgs', panels2args, 'TimeWindow', tWin);
            NP.PTBRead.PlotSE(se, panels3, 'PanelArgs', panels3args, 'TimeWindow', tWin);
            
            % Save figure
            recId = string(NP.SE.GetID(ud));
            region = NP.SE.GetRegion(ud(1));
            figDir = fullfile(NP.Data.GetAnalysisRoot, "se_"+region, recId);
            if ~exist(figDir, 'dir')
                mkdir(figDir);
            end
            % figName = strjoin([recId, "semsr", "page"+p, "u"+uId], "_");
            figName = strjoin([recId, "nbd", "page"+p], "_"); %filename is often too long with all the unitIds included
            if ~exist(figDir, 'dir')
                mkdir(figDir);
            end
            exportgraphics(f, fullfile(figDir, figName+".png"));
        end

        function tv = AddStimulusFeatures(tv)
            % Add word-level features to the taskTime table from the reference spreadsheet in /supporting_files/NBD_a0246a17-sentences-manual.csv

            % Get the directory of the current script
            scriptDir = fileparts(mfilename('fullpath'));

            % Go up one directory to the root of the project
            scriptDir = fullfile(scriptDir, '..');

            % Load the reference spreadsheet
            ref = readtable(fullfile(scriptDir, 'supporting_files', 'NBD_a0246a17-sentences-manual.csv'), 'TextType', 'string');

            features_to_copy = {'ortho_upoint', 'phono_upoint', 'sentence_length', 'word_length', 'ortho_n_dens_s', 'ortho_n_freq_s_m', 'phono_n_dens_s', 'phono_n_freq_s_m', 'sum_biphone', 'sum_bigram', 'n_phon', 'n_syll', 'pos'};
            % for i = 1 : length(features_to_copy)
            %     tv.(features_to_copy{i}) = NaN(height(tv), 1);
            % end
            % 
            % % Also add the new features
            % new_features = {'block', 'condition', 'OND', 'PND', 'surprisal', 'modality'};
            % for i = 1 : length(new_features)
            %     tv.(new_features{i}) = NaN(height(tv), 1);
            % end


            % Match the stimText in the taskTime table with the reference spreadsheet by comparing sentences (after removing all punctuation and capitalization)
            tv.stimText = lower(tv.stimText);
            ref.sentence = lower(ref.sentence);
            % But we want to keep apostrophes and hyphens
            tv.stimText = regexprep(tv.stimText, '[^\w\s''-]', '');
            ref.sentence = regexprep(ref.sentence, '[^\w\s''-]', '');
            tv.stimText = regexprep(tv.stimText, '\s+', ' ');
            ref.sentence = regexprep(ref.sentence, '\s+', ' ');

            for i = 1 : height(tv)
                if(i == 45)
                    tv.stimText(i) = "that's a no-entry sign"; %This sentence gets parsed incorrectly?
                    row = 55;
                elseif(i == 83)
                    tv.stimText(i) = "we're past the getting to know you phase"; %This sentence is different in the reference spreadsheet (hyphens are different)
                    row = 84;
                else
                    % Find the matching sentence in the reference spreadsheet
                    row = find(strcmp(tv.stimText(i), ref.sentence));
                end

                if isempty(row)
                    fprintf("No match for sentence: %s\n", tv.stimText(i));
                    continue
                end

                % tv.surprisal(i) = ref.label(row);
                tv.block(i) = ref.block(row);
                tv.condition{i} = ref.our_cond{row};
                % Condition is formatted as 'O+P- high', 'O-P- low', etc, where both +/- and high/low can vary independently. Parse these out to get the orthographic and phonological conditions.
                tv.OND(i) = contains(tv.condition(i), 'O+');
                tv.PND(i) = contains(tv.condition(i), 'P+');
                tv.surprisal(i) = contains(tv.condition(i), 'high');

                % Remove the surprisal part of the condition
                tv.condition(i) = regexprep(tv.condition(i), ' high| low', '');

                % Assign modality based on block membership: 0, 2 = auditory; 1, 3 = visual
                modalities = {'auditory', 'visual'};
                tv.modality{i} = modalities{1 + mod(tv.block(i), 2)};

                % Copy the rest of the features
                for j = 1 : length(features_to_copy)
                    tv.(features_to_copy{j})(i) = ref.(features_to_copy{j})(row);
                end
            end
        end

        function [tt, tv] = FindKeyTimes(tt, tv)
            % Get onset/offset times for the first, second, penultimate, and last words in each sentence, and add them to tt

            relevant_fields = {'firstWordOn', 'firstWordOff', 'secondWordOn', 'secondWordOff', 'penultimateWordOn', 'penultimateWordOff', 'lastWordOn', 'lastWordOff'};
            for i = 1:length(relevant_fields)
                tt.(relevant_fields{i}) = NaN(height(tt), 1);
            end
            
            for i = 1 : height(tt)
                estim = tt.stim(i);
                if iscell(estim)
                    estim = estim{1};
                end
                silent_intervals = false;
                if silent_intervals
                    % Note that for reading the even-index intervals in the TextGrid are empty strings corresponding to the transition periods between words.
                    tt.('firstWordOn')(i) = estim(1).T.T1(1);
                    tt.('firstWordOff')(i) = estim(1).T.T2(1);
                    tt.('secondWordOn')(i) = estim(1).T.T1(3);
                    tt.('secondWordOff')(i) = estim(1).T.T2(3);
                    tt.('penultimateWordOn')(i) = estim(1).T.T1(end-2);
                    tt.('penultimateWordOff')(i) = estim(1).T.T2(end-2);
                    tt.('lastWordOn')(i) = estim(1).T.T1(end);
                    tt.('lastWordOff')(i) = estim(1).T.T2(end);
                else
                    tt.('firstWordOn')(i) = estim(1).T.T1(1);
                    tt.('firstWordOff')(i) = estim(1).T.T2(1);
                    tt.('lastWordOn')(i) = estim(1).T.T1(end);
                    tt.('lastWordOff')(i) = estim(1).T.T2(end);
                end

            end
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

        function seSort = SortTrials(se)
            % Sort trials by modality and condition

            tv = se.GetTable('taskValue');

            % Make categorical version of modality
            modality_list = unique(tv.modality);
            tv.modality_categorical = categorical(tv.modality, modality_list, 'Ordinal', true);

            % Make categorical version of condition
            condition_list = unique(tv.condition);
            tv.condition_categorical = categorical(tv.condition, condition_list, 'Ordinal', true);

            % Make categorical version of surprisal
            surprisal_list = unique(tv.surprisal);
            tv.surprisal_categorical = categorical(tv.surprisal, surprisal_list, 'Ordinal', true);

            % And convert surprisal to "high" or "low" based on the whether the value is 0 or 1
            surprisals = {'low', 'high'};
            for i = 1 : length(surprisals)
                tv.surprisal_text(tv.surprisal == i) = surprisals(i);
            end
            % Make a label for each trial based on modality and condition.
            tv.condition_label = strcat(tv.modality, ' ',tv.condition);
            tv.surprisal_label = strcat(tv.modality, ' ',tv.surprisal_text);


            
            [tv, I] = sortrows(tv, {'modality_categorical', 'surprisal_categorical'});
            seSort = se.Duplicate;
            seSort.SortEpochs(I);
            seSort.SetTable('taskValue', tv);
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

