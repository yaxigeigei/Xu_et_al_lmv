classdef TaskBaseClass
    
    methods(Static)
        % Utilities
        function MigrateStimId(se, pattern)
            % Remove columns of NP.TGEvent objects from taskTime table whose names match the 'pattern', and 
            % add the names as stimId to the taskValue table.
            % Note that because the ID events often preceed the trimmed stim events. This function therefore 
            % matches ID events with stim events by finding the pair closest in time.
            % 
            %   MigrateStimId(se, pattern)
            % 
            % Inputs
            %   se          A MSessionExplorer object.
            %   pattern     A regular expression pattern.
            %               '^.{5}_si.{3,}' works with most TIMIT stim IDs.
            %               '^(ANT_|CAT_|SEMSR\d*_)' works for Semsr stim IDs.
            % 
            
            % Collect all valid events
            seLite = se.Duplicate({'taskTime'}, false);
            seLite.SliceSession(0, 'absolute');
            tt = seLite.GetTable('taskTime');
            varNames = string(tt.Properties.VariableNames);
            sid = cell(size(varNames));
            for i = numel(varNames) : -1 : 1
                vn = varNames(i);
                isRm(i) = ~isempty(regexp(vn, pattern, 'once'));
                if isRm(i)
                    if iscell(tt.(vn))
                        v = tt.(vn){1};
                    else
                        v = tt.(vn);
                    end
                    if class(v) == "MSessionExplorer.Event"
                        sid{i} = v;
                    end
                end
            end
            sid = cat(1, sid{:});
            sid(isnan(sid)) = [];
            
            % Remove ID related columns
            se.SetColumn('taskTime', isRm, []);
            
            % Check the presence of valid IDs
            if isempty(sid)
                return
            end
            
            % Match sid events to the closest stim events
            tv = se.GetTable('taskValue');
            stim = se.GetColumn('taskTime', 'stim');
            if iscell(stim)
                stim = cellfun(@(x) x(1), stim);
            end
            stim = stim + se.GetReferenceTime('taskTime');
            tStim = double(stim);
            tSID = double(sid);
            for i = 1 : numel(sid)
                [~, I] = min(abs(tStim - tSID(i)));
                tv.stimId(I) = sid(i).GetVfield('type');
            end
            
            se.SetTable('taskValue', tv);
        end
        
        function AddEventObjects(se, sDef, tbName)
            % Add task phase events to the taskTime table in se
            % 
            %   AddEventObjects(se, sDef)
            %   AddEventObjects(se, sDef, tbName)
            % 
            % Inputs
            %   se          MSessionsExplore object.
            %   sDef        A struct whose field names are the names of phase events and the field values 
            %               are 2-element cell arrays indicating the name of the onset and offset event.
            %   tbName      The table to read and save events. Default is 'taskTime'.
            % 
            
            % Get table
            if ~exist('tbName', 'var')
                tbName = 'taskTime';
            end
            tb = se.GetTable(tbName);
            
            % Make all columns cell array
            cellTb = tb;
            for i = 1 : width(tb)
                if ~iscell(cellTb.(i))
                    cellTb.(i) = num2cell(cellTb.(i));
                end
            end
            
            % Create events in taskTime table
            evtList = fieldnames(sDef);
            for i = 1 : height(tb)
                for j = 1 : numel(evtList)
                    % Get variables
                    en = evtList{j};
                    vn = sDef.(en);
                    t1 = cellTb.(vn{1}){i}; % onsets
                    t2 = cellTb.(vn{2}){i}; % offsets
                    
                    % Find the closest offset (before the next onset) for each onset time
                    t12 = [t1 NaN(size(t1))];
                    for n = 1 : numel(t1)
                        if n < numel(t1)
                            k = find(t1(n) < t2 & t2 < t1(n+1), 1);
                        else
                            k = find(t1(n) < t2 & t2 < t1(n)+10, 1);
                        end
                        if ~isempty(k)
                            t12(n,2) = t2(k);
                        end
                    end
                    
                    % Remove events where offset could not be found
                    isOrphan = any(isnan(t12), 2);
                    t12(isOrphan,:) = [];
                    if isempty(t12)
                        t12 = [NaN NaN];
                    end
                    
                    % Construct event objects
                    evt = MSessionExplorer.Event(t12(:,1));
                    evt = evt.SetTfield('tOn', t12(:,1));
                    evt = evt.SetTfield('tOff', t12(:,2));
                    tb.(en){i} = evt;
                end
            end
            se.SetTable(tbName, tb);
        end
        
        function [t1, t2, alignInfo] = FindMatchedSpeechTimes(e1, e2)
            % Find matching times from two sequences of TGEvent objects by global alignment
            % 
            %   [t1, t2, alignInfo] = FindMatchedSpeechTimes(e1, e2)
            % 
            % See also MLing.FindAlignedTokens
            
            % Get sequences of speech labels
            s1 = e1.GetParentLabel();
            s2 = e2.GetParentLabel();
            
            % Align labels
%             [k1, k2, alignInfo] = MLing.FindAlignedWords(s1, s2);
            [k1, k2, alignInfo] = MLing.FindAlignedTokens(s1, s2);
            alignInfo.kk = [k1 k2]';
            
            % Get onset and offset times of aligned events
            % offsets are subtracted by 1 millisec to avoid overlap with the following onsets
            e1 = e1(k1);
            t1 = e1.GetTfield('tmin');
            e2 = e2(k2);
            t2 = e2.GetTfield('tmin');
            
            % Delete duplicated times
            % e.g. if word A in e1 matches words X and Y in e2, A will be duplicated for each of X and Y
            [~, ind1] = unique(t1(:), 'stable');
            [~, ind2] = unique(t2(:), 'stable');
            ind = intersect(ind1, ind2);
            t1 = t1(ind);
            t2 = t2(ind);
        end
        
        function tt = FindMorphTimes(tt, eventNames)
            % Find matching key times between source trial and median sequence event temaplte
            % 
            %   tt = FindMorphTimes(tt, eventNames)
            % 
            
            T = tt{:,eventNames};
            dT = diff(T, 1, 2);
            
            % Exclude out-of-order event times by setting them to NaN
            for i = 1 : size(T,1)
                for j = 2 : size(T,2)
                    if any(T(i,j) <= T(i,1:j-1))
                        T(i,j) = NaN;
                        fprintf("Exclude the #%i event in epoch %i since it preceeds or overlap with a previous event.\n", j, i);
                    end
                end
            end
            
            % Construct template event times from median onset and intervals
            dT = median(dT, 1, 'omitnan');
            Tm = cumsum([0 dT]);
            
            % Get NaN mask based on both T and Tm
            M = ~isnan(T + Tm);
            
            for i = 1 : height(tt)
                tt.morphFrom{i} = T(i, M(i,:))';
                if(isempty(tt.morphFrom{i}))
                    disp(tt(i,:));
                end
                tt.morphTo{i} = tt.morphFrom{i}(1) + Tm(M(i,:))';
            end
        end
        
        function seTb = SplitBySentence(se)
            % Split se by unique sentences according to stimText in the taskValue table
            % 
            %   seTb = SplitBySentence(se)
            % 
            
            % Split
            seTb = NP.SE.SplitConditions(se, 'conditionVars', {'taskName', 'stimText', 'stimId'});
            
            % Add group info
            for i = 1 : height(seTb)
                tv = seTb.se(i).GetTable('taskValue');
                seTb.trialNum{i} = tv.trialNum;
                if ismember('tempTrialNum', tv.Properties.VariableNames)
                    seTb.tempTrialNum(i) = tv.tempTrialNum(1);
                else
                    seTb.tempTrialNum(i) = NaN;
                end
            end
        end
        
        function player = PlayAudio(se, trialIdx, tWin, colName)
            % Play mic or speaker audio
            
            if nargin < 4
                colName = 'mic';
            end
            if nargin < 3 || isempty(tWin)
                % Take data as is
                ni = se.GetTable('ni');
                ni = ni(trialIdx,:);
            else
                % Reslice data if time window is provided
                ni = se.SliceTimeSeries('ni', tWin, trialIdx, 'Fill', 'bleed');
            end
            
            % Play audio waveform
            fs = se.userData.niMeta.niSampRate;
            if ischar(fs)
                fs = str2double(fs);
            end
            wav = ni.(colName){1};
            player = audioplayer(wav, fs);
            player.play();
        end
        
        % General plots
        function axArray = PlotSE(se, panels, varargin)
            % Plot various speech features, unit raster, and PETHs aligned in time
            % 
            %   PlotSE(se, panels)
            %   PlotSE(se, panels, rowDist)
            %   PlotSE(se, panels, rowDist, colDist)
            %   PlotSE(..., 'PanelArgs', {{}})
            %   PlotSE(..., 'TimeWindow', {'trialOn', 'matchOff'})
            %   PlotSE(..., 'TimeWindowOffsets', [-0.3 0.3])
            %   PlotSE(..., 'ExampleTrial', [])
            % 
            % Inputs
            %   se              A MSessionExplorer object.
            %   panels          A m-by-n cell array. The shape specifies that the figure has m rows and n 
            %                   columns of plots (or Axes). The variable in each cell specifies what to 
            %                   plot, and can be one of the following:
            %                   1) A char string of plot name. Supported plots include 'words', 'mels', 
            %                      'saccades', 'reacts'.
            %                   2) A positional index of units in the spikeTime table.
            %                   3) Empty, then the Axes will be left empty.
            %                   The default panelList is {'words'}.
            %   rowDist         A m-element vector specifying integer height (i.e. row) ratios of the panels.
            %   colDist         A n-element vector specifying integer width (i.e. column) ratios of the panels.
            %   'PanelArgs'     A cell array of structs. The size of this array should match that of panels. 
            %                   Each struct contains the additional arguments for plotting. Therefore, the 
            %                   recepient function must support struct expansion of Parameter-Value pairs.
            %   'TimeWindow'        1) A n-by-2 numeric array of time windows. n is the number of trials.
            %                          If n == 1, this window will be applied to all trials.
            %   'TimeWindowOffsets' ...
            % 
            
            % Parse user inputs
            p = inputParser();
            p.addOptional('rowDist', [], @isnumeric);
            p.addOptional('colDist', [], @isnumeric);
            p.addParameter('TimeWindow', {'trialOn', 'matchOff'}, @(x) iscellstr(x) || isstring(x) || isnumeric(x));
            p.addParameter('TimeWindowOffsets', [-0.3 0.3], @isnumeric);
            p.addParameter('ExampleTrial', [], @isnumeric);
            p.addParameter('PanelArgs', {{}}, @iscell);
            p.parse(varargin{:});
            rowDist = p.Results.rowDist;
            colDist = p.Results.colDist;
            panelArgs = p.Results.PanelArgs;
            tWin = p.Results.TimeWindow;
            tWinOffsets = p.Results.TimeWindowOffsets;
            egTrial = p.Results.ExampleTrial;
            
            % Resolve the figure layout
            if isempty(rowDist)
                rowDist = ones(size(panels,1), 1);
            end
            if isempty(colDist)
                colDist = ones(size(panels,2), 1);
            end
            f = gcf;
            isTl = arrayfun(@(x) isa(x, 'matlab.graphics.layout.TiledChartLayout'), f.Children);
            if ~any(isTl)
                tl = tiledlayout(sum(rowDist), sum(colDist));
                tl.Padding = 'tight';
            else
                tl = f.Children(isTl);
            end
            
            % Match arguments size to panels
            if isscalar(panelArgs)
                panelArgs = repmat({{}}, size(panels));
            end
            
            % Find the index of the template trial
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            if isempty(egTrial) && ismember('tempTrialNum', tv.Properties.VariableNames)
                egTrial = find(tv.trialNum==tv.tempTrialNum(1), 1);
            end
            if isempty(egTrial)
                egTrial = 1;
            end
            
            % Find time windows
            if iscellstr(tWin) || isstring(tWin)
                tWin = cellstr(tWin);
                [onEvt, offEvt] = tWin{:};
                tWin = [tt.(onEvt) tt.(offEvt)] + tWinOffsets;
            end
            if size(tWin,1) < se.numEpochs
                tWin = repmat(tWin, [se.numEpochs 1]);
            end
            
            % Features
            axArray = cell(size(panels));
            for i = 1 : size(panels, 1)
                for j = 1 : size(panels, 2)
                    % Create Axes
                    ntArgs = MPlot.FindTileInd(rowDist, colDist, i, j);
                    ax = nexttile(tl, ntArgs{:});
                    hold(ax, 'on');
                    
                    % Plot a panel
                    p = panels{i,j};
                    args = panelArgs{i,j};
                    if isempty(p)
                        % plot nothing
                        
                    elseif isa(p, 'function_handle')
                        % Invoke custom function
                        p(ax, se, egTrial, tWin(egTrial,:), args{:});
                        
                    elseif ischar(p)
                        % Plot feature
                        switch p
                            case 'phone'
                                if isempty(args)
                                    args = {{'stim', 'prod'}};
                                end
                                NP.TaskBaseClass.PlotTGEHier(ax, se, egTrial, tWin(egTrial,:), args{:});
                            case 'phones'
                                if isempty(args)
                                    args = {"phone", {'stim', 'prod'}};
                                end
                                NP.TaskBaseClass.PlotTGE(ax, se, [], tWin, args{:});
                                
                            case 'mel'
                                NP.TaskBaseClass.PlotMelSpectrogram(ax, se, egTrial, tWin(egTrial,:), args{:});
                            case 'mels'
                                NP.TaskBaseClass.PlotMelSpectrograms(ax, se, [], tWin);
                            case 'acous'
                                NP.TaskBaseClass.PlotIntensity(ax, se, egTrial, tWin(egTrial,:));
                            case 'pitch'
                                NP.TaskBaseClass.PlotPitch(ax, se, egTrial, tWin(egTrial,:), 'Background', 'mel');
                            case 'fujisaki'
                                NP.TaskBaseClass.PlotFujisaki(ax, se, egTrial, tWin(egTrial,:));
                                
                            case 'saccades'
                                NP.TaskBaseClass.PlotTrialSaccade(ax, se, [], tWin);
                                
                            case 'event_boundary'
                                if isempty(args)
                                    args = {{'stimOff', 'cue3On', 'cue3Off'}};
                                end
                                NP.TaskBaseClass.PlotEventBoundary(ax, se, [], tWin(egTrial,:), args{:});
                            case 'block_boundary'
                                if isempty(args)
                                    args = {{'stimText'}};
                                end
                                NP.TaskBaseClass.PlotBlockBoundary(ax, se, [], tWin(egTrial,:), args{:});
                                
                            case 'raster'
                                NP.UnitPlot.Raster(ax, se, [], tWin(egTrial,:), args{:});
                            case 'rate_map'
                                NP.UnitPlot.Heatmap(ax, se, [], tWin(egTrial,:), args{:});
                            case 'rasters'
                                NP.UnitPlot.RasterStack(ax, se, [], tWin(egTrial,:), args{:});
                            case 'rate_maps'
                                NP.UnitPlot.HeatmapStack(ax, se, [], tWin(egTrial,:), args{:});
                                
                            case 'cla'
                                cla(ax);
                        end
                    end
                    
                    % Omit x-label if it's not the last row
                    if i ~= size(panels, 1)
                        ax.XLabel.String = [];
                    end
                    
                    axArray{i,j} = ax;
                end
            end
            
            axArray = [axArray{:}];
            axArray = reshape(axArray, size(panels));
        end
        
        function PlotEvents(ax, se, trialIdx, tWin, colNames)
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
        
        function PlotEventBoundary(ax, se, trialInd, tWin, colNames, cmapFunc)
            % Plot lines of the same events across trials, typically on top of raster
            % 
            %   PlotEventBoundary(ax, se, trialInd, tWin, colNames, cmapFunc)
            % 
            
            if ~exist('cmapFunc', 'var')
                cmapFunc = @lines;
            end
            
            colNames = cellstr(colNames);
            tt = se.SliceEventTimes('taskTime', tWin, trialInd, colNames); % not using 'bleed' otherwise TGEs will get converted to double
            
            if isempty(trialInd)
                y = (1 : se.numEpochs)';
            else
                y = (1 : numel(trialInd))';
            end
            if isscalar(y)
                return
            end
            y([1 end]) = y([1 end]) + [-1 1]'*0.5;
            cc = cmapFunc(numel(colNames));
            
            hold(ax, 'on');
            for k = 1 : numel(colNames)
                evt = tt.(colNames{k});
                if iscell(evt)
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
        
        function PlotBlockBoundary(ax, se, trialInd, tWin, colNames, varargin)
            % Plot seperator between different trial blocks
            % 
            %   PlotBlockBoundary(ax, se, trialInd, tWin, colNames)
            %   PlotBlockBoundary(ax, se, trialInd, tWin, colNames, ..., 'Label', true)
            %   PlotBlockBoundary(ax, se, trialInd, tWin, colNames, ..., 'Color', lines(numel(colNames)))
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
            ax.YLabel.String = 'Trials';
            MPlot.Axes(ax);
        end
        
        function PlotEventWindows(ax, se, trialIdx, tWin, colNames, varargin)
            % Plot task and speech event objects for a single trial
            % 
            %   PlotEventWindows(ax, se, trialIdx, tWin, colNames)
            %   PlotEventWindows(..., 'Style', 'block')
            %   PlotEventWindows(..., 'YRange', trialIdx+[-0.5 0.5])
            %   PlotEventWindows(..., 'Colors', lines(numel(colNames)))
            %   PlotEventWindows(..., 'Alpha', 0.2)
            %   PlotEventWindows(..., 'Text', false)
            %   PlotEventWindows(..., 'TextArgs', {})
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
            p.addParameter('TextArgs', {}, @iscell);
            p.parse(varargin{:});
            style = p.Results.Style;
            cc = p.Results.Colors;
            yRange = p.Results.YRange;
            alpha = p.Results.Alpha;
            isText = p.Results.Text;
            txtArgs = p.Results.TextArgs;
            
            if isscalar(alpha) && numel(colNames) > 1
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
                if isnan(evt)
                    continue
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
                    text(ax, x, y, ss, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', txtArgs{:});
                end
            end
            
            m = ~isinf(tWin);
            ax.XLim(m) = tWin(m);
            ax.XLabel.String = 'Time (s)';
            ax.YLim(1) = min([yRange(:); ax.YLim(1)]);
            ax.YLim(2) = max([yRange(:); ax.YLim(2)]);
            MPlot.Axes(ax);
        end
        
        % Phonetics plots
        function PlotTGE(ax, se, trialInd, tWin, varargin)
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
                            tge.Plot(k+[-.5 .5], 'Parent', ax, p.Unmatched);
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
        
        function PlotTGEHier(ax, se, trialIdx, tWin, colNames, varargin)
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
            tge.PlotTiers(varargin{:}, 'Parent', ax);
            
            if ~isempty(tWin)
                ax.XLim = tWin;
            end
            ax.XLabel.String = 'Time (s)';
            ax.YDir = 'reverse';
            ax.YTick = 1 : tgGetNumberOfTiers(tge(1))+1;
            ax.YTickLabel = {'Sent', 'Word', 'Phone'};
            MPlot.Axes(ax);
        end
        
        % Acoustics plots
        function PlotMelSpectrogram(ax, se, trialIdx, varargin)
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
            p.addParameter('ChannelName', "mic", @(x) ischar(x) || isstring(x));
            p.addParameter('WaveformHeight', 0.15, @isnumeric);
            p.addParameter('ShowIntensity', false, @islogical);
            p.parse(varargin{:});
            tWin = p.Results.tWin;
            chan = string(p.Results.ChannelName);
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
                    tbs{i} = se.SliceTimeSeries(tbNames{i}, tWin, trialIdx, 'Fill', 'none');
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
            S = mel.(chan){1};
            F = se.userData.melMeta.F;
            T = mel.time{1};
            NP.Audio.PlotMelSpectrogram(ax, S, F, T);
            hold(ax, 'on');
            melMax = max(hz2mel(F));
            
            ni = tbs{2};
            if ~isempty(ni)
                % Plot waveform
                tNI = ni.time{1};
                w = ni.(chan){1};
                tf = @(y) y/max(abs(w)) * hWf/2*melMax + (1-hWf/2)*melMax;
                w = tf(w);
                plot(ax, tNI, w, 'k');
                
                it = tbs{3};
                if ~isempty(it)
                    % Plot envelope, its peaks, and the peaks of d(env)/dt
                    tIt = it.time{1};
                    plot(ax, tIt, tf(it.env{1}), 'r', 'LineWidth', 2);
                    plot(ax, tIt, tf(it.peakEnv{1}/2), 'g', 'LineWidth', 2);
                    plot(ax, tIt, tf(it.peakRate{1}*3), 'y', 'LineWidth', 2);
                end
            end
            
            if ~isempty(tWin)
                ax.XLim = tWin;
            end
%             ax.YLim = [0 max(yNI)];
            ax.XLabel.String = 'Time (s)';
            MPlot.Axes(ax);
        end
        
        function PlotIntensity(ax, se, trialIdx, tWin)
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
            cc = lines;
            plot(ax, tIt, it.env{1}, 'Color', cc(2,:), 'LineWidth', 2);
            plot(ax, tIt, it.dEnv{1}/10, 'Color', cc(3,:), 'LineWidth', 2);
            plot(ax, tIt, it.peakRate{1}*4, 'Color', cc(5,:), 'LineWidth', 2);
            
            ax.XLim = tWin;
            ax.XLabel.String = 'Time (s)';
            ax.YLim = [-maxVal maxVal] * 1.2;
            ax.YLabel.String = 'AU';
            MPlot.Axes(ax);
        end
        
        function PlotPitch(ax, se, trialIdx, varargin)
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
            
            if ~ismember('pitch', se.tableNames)
                return
            end
            
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
        
        function PlotRelativePitch(ax, se, trialIdx, varargin)
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
        
        function PlotBinnedPitch(ax, se, trialIdx, varargin)
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
        
        function PlotFujisaki(ax, se, trialIdx, varargin)
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
            ax.YLim = yRange;
%             ax.YLim = [0 10];
            ax.YLabel.String = "Relative F0";
            MPlot.Axes(ax);
        end
        
        function PlotMelSpectrograms(ax, se, trialInd, varargin)
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
        
        % Facial tracking plots
        function PlotFacialFeatures(ax, se, trialIdx, tWin)
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
        
        function PlotSaccade(ax, se, trialInd, tWin)
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
        
        % Time morphing plots
        function PlotAlignment(seTb, varargin)
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
        
        function PlotAlignmentScores(ax, se, trialInd)
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
        
        function PlotAlignedWords(ax, se, trialInd, tWin, colName)
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
        
        function PlotSpeedMap(ax, se, trialInd, tWin)
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

