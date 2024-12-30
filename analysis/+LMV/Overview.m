classdef Overview
    % Make plots for overview
    
    methods(Static)
        function OneTrial(se, tr, tWin)
            % Plot speech labels, spectrogram, and some features for inspection
            % 
            %   OneTrial(se, tr)
            %   OneTrial(se, tr, tWin)
            % 
            % Inputs
            %   se          An MSeesionExplorer object.
            %   tr          Trial index.
            %   tWin        'full', 'stim', 'prod', or a two-element numeric vector of start and end time.
            % 
            
            tt = se.GetTable('taskTime');
            
            % Time window around production
            if ~exist('tWin', 'var') || isempty(tWin)
                tWin = 'full';
            end
            if ~isnumeric(tWin)
                switch tWin
                    case 'full'
                        tWin = [tt.cue1On(tr) max(tt.cue3Off(tr), tt.prod{tr}(end).T.tmax)] + [-1 1]*0.2;
                    case 'stim'
                        tWin = [tt.cue1On(tr) tt.stimOff(tr)] + [-1 1]*0.1;
                    case 'prod'
                        tWin = [tt.prod{tr}(1).T.tmin tt.prod{tr}(end).T.tmax] + [-1 1]*0.1;
                end
            end
            assert(~any(isnan(tWin)), "The time window contains NaN");
            
            % Plot
            tl = tiledlayout(5, 1);
            tl.Padding = 'compact';
            k = 1;
            
            ax = nexttile(k); cla; k = k + 1;
            NP.TaskBaseClass.PlotTGEHier(ax, se, tr, tWin, {'stim', 'prod'});
            ax.XLabel.String = [];
            ax.XTickLabel = [];
            
            ax = nexttile(k, [2 1]); cla; k = k + 2;
            NP.TaskBaseClass.PlotMelSpectrogram(ax, se, tr, tWin, 'Wave', 0.3, 'ShowIntensity', false);
            ax.XLabel.String = [];
            ax.XTickLabel = [];
            
            ax = nexttile(k, [2 1]); cla; k = k + 2;
            NP.TaskBaseClass.PlotPitch(ax, se, tr, tWin, 'Background', 'mel');
            
            % ax = nexttile(k, [2 1]); cla; k = k + 1;
            % NP.TaskBaseClass.PlotRelativePitch(ax, se, tr, tWin);
            
            % ax = nexttile(k); cla; k = k + 1;
            % NP.TaskBaseClass.PlotFujisaki(ax, se, tr, tWin);
        end
        
        function Session(se, varargin)
            % Plot session heatmap rasters
            % 
            %   Overview(se)
            %   Overview(se, ..., 'UnitIds', [])
            %   Overview(se, ..., 'UnitsPerPage', 10, 'Page', 0)
            %   Overview(se, ..., 'TaskPhase', 'full')
            %   Overview(se, ..., 'Folder', [])
            % 
            % Inputs
            %   se                  An MSessionExplorer object.
            %   
            %   There are two ways to specify which units to plot:
            %   1) Plot specified units using IDs.
            %       'UnitIDs'       A vector of cluster IDs. IDs should be consistent with the ID format (i.e. origianl 
            %                       KS short ID or unique long ID) in se.userData.ksMeta.clusTb.clusId. NaN IDs are 
            %                       legal and will be ignored.
            %   2) Plot a batch of units from all units on a single page/figure.
            %       'Page'          The index of the page (batch) to plot.
            %       'UnitsPerPage'  The number of units to be plotted on the page.
            %   
            %   'TaskPhase'         The task phase to set the time limits for. Support options include:
            %                         'full': from stimOn to last prodOff, with some padding.
            %   'Folder'            The path of a folder to save the plots as images. Default is empty [], and no plot
            %                       will be saved.
            %   
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('UnitIds', [], @isnumeric);
            p.addParameter('UnitsPerPage', 10, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Page', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('TaskPhase', 'full', @(x) ischar(x) || isstring(x));
            p.addParameter('PlotFun', 'rate_map', @(x) ischar(x) || isa(x, "function_handle"));
            p.addParameter('Folder', [], @(x) ischar(x) || isstring(x) || isempty(x));
            p.parse(varargin{:});
            uId = p.Results.UnitIds;
            upp = p.Results.UnitsPerPage;
            pg = p.Results.Page;
            taskPhase = p.Results.TaskPhase;
            fPlot = p.Results.PlotFun;
            figDir = p.Results.Folder;
            
            % Determine the units to plot
            ud = se.userData;
            ut = NP.Unit.GetClusTb(ud);
            if isempty(uId)
                % By # of units per page and page index
                uInd = (1:upp)+upp*(pg-1);
                uInd(uInd > height(ut)) = [];
            else
                % By unit ID
                uId(isnan(uId)) = [];
                uInd = NP.Unit.ClusId2Ind(uId, ud);
            end
            uInd(isnan(uInd)) = [];
            
            % Configure plots
            nr = 2; % number of rows
            ppr = ceil(upp/nr) + 1; % plot per row = #unit / 3 rows + one label plot
            panels1 = repmat({fPlot}, [nr ppr]);
            panels1args = repmat({{}}, [nr ppr]);
            
            ebArgs = { ...
                ["cue1On", "stimOn", "stimOff", "cue3On", "prodMatchOn", "prodMatchOff"], ...
                @(x) LMV.Param.GetTaskPhaseColors(["atten", "stim", "stim", "init", "prod", "prod"]) ...
                };
            panels2 = repmat({'event_boundary'}, size(panels1));
            panels2args = repmat({ebArgs}, size(panels1));
            
            bbArgs = {"stimText", 'Label', false, 'Color', [0 0 0 .3]};
            panels3 = repmat({'block_boundary'}, size(panels1));
            panels3args = repmat({bbArgs}, size(panels1));
            
            k = 1;
            for i = 1 : size(panels1,1)
                % first column for labels
                panels1{i,1} = ''; % not plotting spikes in the first column
                panels2{i,1} = ''; % not plotting event_boundary in the first column
                panels3args{i,1}{3} = true; % show labels in the first column
                
                % following columns for rasters
                for j = 2 : size(panels1,2)
                    if k > numel(uInd) || isnan(uInd(k))
                        panels1{i,j} = '';
                        panels2{i,j} = '';
                        panels3{i,j} = '';
                        continue
                    end
                    panels1args{i,j} = {uInd(k)};
                    k = k + 1;
                end
            end
            
            tWin = [-0.3 6];
            
            % Plotting
            f = gcf; clf
            f.WindowState = 'maximized';
            NP.TaskBaseClass.PlotSE(se, panels1, 'PanelArgs', panels1args, 'TimeWindow', tWin);
            NP.TaskBaseClass.PlotSE(se, panels2, 'PanelArgs', panels2args, 'TimeWindow', tWin);
            NP.TaskBaseClass.PlotSE(se, panels3, 'PanelArgs', panels3args, 'TimeWindow', tWin);
            
            % Save figure
            if isempty(figDir)
                return
            end
            if ~exist(figDir, 'dir')
                mkdir(figDir);
            end
            uId = NP.Unit.Ind2ClusId(uInd, se);
            uId = uId - NP.Unit.GetBaseClusId(se);
            figName = strjoin([NP.SE.GetID(se), "lmv", "page"+pg, "u"+uId(:)'], "_");
            exportgraphics(f, fullfile(figDir, figName+".png"));
        end
        
        function SessionFromCache(unitId, varargin)
            % Plot session rasters
            % 
            %   SessionFromCache(unitId)
            %   SessionFromCache(unitId, ...)
            % 
            % See also LMV.Overview.Session
            
            unitId = unique(unitId, 'stable');
            
            % Load cache
            sUnit = NP.Unit.LoadUnitCache(unitId, varargin{:});
            
            % Put data into se
            se = LMV.SE.UnitCache2SE(sUnit, varargin{:});
            
            % Add unit attributes to clusTb
            clusTb = NP.Unit.GetClusTb(se);
            
            % Plot
            if height(clusTb) > 16
                disp("Only plotting the first 16 selected units.");
            end
            upp = MMath.Bound(height(clusTb), [10 16]);
            LMV.Overview.Session(se, 'UnitId', clusTb.clusId, 'UnitsPerPage', upp, varargin{:});
        end
        
        function SessionOnBrush(fig, axesStruct)
            % 
            %   SessionOnBrush(fig, axesStruct)
            % 
            
            % Find brushed units
            hh = axesStruct.Axes.Children;
            b = [];
            for i = 1 : numel(hh)
                % Check if handle is the plot for brushing
                s = hh(i).UserData;
                if ~isfield(s, 'forBrush') || ~s.forBrush
                    continue
                end
                
                % Check if brushed any data
                b = logical(hh(i).BrushData);
            end
            if ~any(b)
                return
            end
            
            % Create figure if absent
            if isempty(fig.UserData) || ~isvalid(fig.UserData.cbFig)
                fig.UserData.cbFig = MPlot.Figure( ...
                    'Name', 'Selected Units', ...
                    'NumberTitle', 'off', ...
                    'Menubar', 'none', ...
                    'Toolbar', 'figure', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
            end
            
            b = find(b);
            if numel(b) > 50
                b = randsample(b, 50);
            end
            uTb = axesStruct.Axes.UserData.clusTb;
            clusId = uTb.clusId(b);
            
            figure(fig.UserData.cbFig);
            LMV.Overview.SessionFromCache(clusId, 'DataSource', 'm2', 'NumTrials', "mode");
        end
        
        function Sentences(senTb, varargin)
            % Plot sentence rasters aligned with task and speech features.
            % 
            %   Sentences(senTb)
            %   Sentences(..., 'UnitIds', [])
            %   Sentences(..., 'UnitsPerPage', [], 'Page', 0)
            %   Sentences(..., 'StimIdList', LMV.Param.stimIdList4)
            %   Sentences(..., 'Features', {'phone', 'mel'})
            %   Sentences(..., 'TaskPhase', 'full')
            %   Sentences(..., 'MinTaskHight', 10)
            %   Sentences(..., 'UnitLabels', ["depth"])
            %   Sentences(..., 'Folder', [])
            % 
            % Inputs
            %   senTb               A table of MSessionExplorer object.
            %   
            %   There are two ways to specify which units to plot:
            %   1) Plot specified units using IDs.
            %       'UnitIds'       A vector of cluster IDs. IDs should be consistent with the ID format (i.e. origianl 
            %                       KS short ID or unique long ID) in se.userData.ksMeta.clusTb.clusId. NaN IDs are 
            %                       legal and will be ignored.
            %   2) Plot a batch of units from all units on a single page/figure.
            %       'Page'          The index of the page (batch) to plot.
            %       'UnitsPerPage'  The number of units to be plotted on the page.
            %   
            %   'StimIdList'        A string array of stimuli to be included. Each stimulus plots a column.
            %   'Features'          
            %   'TaskPhase'         The task phase to set the time limits for. Support options include
            %                         'stim': 
            %                         'prod': 
            %                         'full': from stimOn to last prodOff, with some padding.
            %   'UnitLabels'        Any combination of "depth", "ttest", "zeta".
            %   'MinTaskHight'      The minimal percentage of height the task label plots should take.
            %   'TaskPhase'         The task phase to set the time limits for. Support options include:
            %   'Folder'            The path of a folder to save the plots as images. Default is empty [], and no plot
            %                       will be saved.
            % 
            
            % Process inputs
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('UnitIds', [], @isnumeric);
            p.addParameter('UnitsPerPage', 15, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Page', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('StimIdList', LMV.Param.stimIdList4, @isstring);
            p.addParameter('NColumns', 4, @isscalar);
            p.addParameter('Features', {'phone', 'mel'}, @(x) isstring(x) || iscellstr(x));
            p.addParameter('TaskPhase', 'full', @(x) ischar(x) || isstring(x));
            p.addParameter('MinTaskHight', 10, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('UnitLabels', ["depth"], @isstring);
            p.addParameter('Folder', [], @(x) ischar(x) || isstring(x));
            p.parse(varargin{:});
            uId = p.Results.UnitIds;
            upp = p.Results.UnitsPerPage;
            pg = p.Results.Page;
            featPanels = cellstr(p.Results.Features(:));
            stimIdList = p.Results.StimIdList;
            nCols = p.Results.NColumns;
            taskPhase = p.Results.TaskPhase;
            minFeatHight = p.Results.MinTaskHight;
            unitLabels = p.Results.UnitLabels;
            figDir = p.Results.Folder;
            
            % Determine the units to plot
            ud = senTb.se(1).userData;
            ut = NP.Unit.GetClusTb(ud);
            if isempty(uId)
                % By # of units per page and page index
                uInd = (1:upp)+upp*(pg-1);
                uInd(uInd > height(ut)) = NaN;
            else
                % By unit ID
                uInd = NP.Unit.ClusId2Ind(uId, ud);
            end
            
            % Configure time window
            switch taskPhase
                case 'full'
                    winEvts = {'trialOn', 'matchOff'};
                    winOffsets = [-.3 .3];
                case 'lmv'
                    winEvts = {'stimOn', 'matchOff'};
                    winOffsets = [-.3 .3];
                case 'stim'
                    winEvts = {'stimOn', 'stimOff'};
                    winOffsets = [-.3 .3];
                case 'prod'
                    winEvts = {'matchOn', 'matchOff'};
                    winOffsets = [-.3 .3];
                otherwise
                    error("'%s' is not a valid 'TaskPhase' parameter.", taskPhase);
            end
            
            % Select sentences
            if isempty(stimIdList)
                [~, ind] = sort(senTb.numTrial, 'descend');
                ind = ind(1:12);
                stimIdList = string(senTb.stimId(ind));
            else
                [~, ind] = MMath.SortLike(senTb.stimId, stimIdList, false);
            end
            senTb = senTb(ind,:);
            
            
            % Compute figure layout
            nStim = numel(stimIdList);
            nRowSets = ceil(nStim/nCols);
            
            nFeats = numel(featPanels);
            hFeats = repelem(round(100/nRowSets/(nFeats+numel(uInd))), nFeats);
            hFeats = max(minFeatHight, hFeats);
            rowDist = repmat([hFeats floor(100/nRowSets)-sum(hFeats)], 1, nRowSets);
            
            blockPanels = [featPanels; {'rasters'}];
            blockPanelArgs = [repmat({{}}, [numel(featPanels) 1]); {{uInd}}]; % TBD: specify unit labels
            nBlockPanels = numel(blockPanels);
            
            % Plot
            f = gcf;
            f.WindowState = 'maximized';
            pause(0.1);
            
            tl = tiledlayout(sum(rowDist), nCols);
            tl.Padding = 'compact';
            
            for i = 1 : nRowSets
                for j = 1 : nCols
                    % Find the current stim
                    k = (i-1)*nCols + j;
                    if k > nStim || k > height(senTb)
                        break
                    end
                    
                    % Initialize with empty panels
                    panels = cell(numel(rowDist), nCols);
                    panelArgs = cell(numel(rowDist), nCols);
                    
                    % Configure all figure panels for the current iteration
                    rowInd = (i-1)*nBlockPanels + (1:nBlockPanels);
                    panels(rowInd, j) = blockPanels;
                    panelArgs(rowInd, j) = blockPanelArgs;
                    
                    % Plot sentence block
                    NP.TaskBaseClass.PlotSE(senTb.se(k), panels, rowDist, 'PanelArgs', panelArgs, 'TimeWindow', winEvts, 'TimeWindowOffsets', winOffsets);
                end
            end
            
            % Save figure
            if isempty(figDir)
                return
            end
            if ~exist(figDir, 'dir')
                mkdir(figDir);
            end
            uId = NP.Unit.Ind2ClusId(uInd, ud);
            uId(isnan(uId)) = [];
            uId = uId - NP.Unit.GetBaseClusId(ud);
            figName = strjoin([NP.SE.GetID(ud), senTb.taskName(1), "page"+pg, "u"+uId(:)'], "_") + ".png";
            exportgraphics(f, fullfile(figDir, figName));
        end
        
        function SentencesFromCache(unitId, varargin)
            % Plot rasters of stim aligned to speech labels
            % 
            %   SentencesFromCache(unitId)
            %   SentencesFromCache(unitId, ...)
            % 
            % See also LMV.Overview.Sentences
            
            % Load cache
            sUnit = NP.Unit.LoadUnitCache(unitId, varargin{:});
            
            % Put data into se
            se = LMV.SE.UnitCache2SE(sUnit, varargin{:});
            
            % Add unit attributes to clusTb
            clusTb = NP.Unit.GetClusTb(se);
            tTb = struct2table(cellfun(@(s) s.tResp, sUnit), "AsArray", true);
            
            se.userData.ksMeta.clusTb = clusTb;
            
            % Make sentence table
            senTb = se.SplitConditions('stimId', 'taskValue');
            
            % Plot
            LMV.Overview.Sentences(senTb, 'Features', "phone", 'UnitId', clusTb.clusId, varargin{:});
        end
        
    end
    
end
