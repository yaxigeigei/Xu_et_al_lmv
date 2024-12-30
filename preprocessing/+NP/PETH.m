classdef PETH
    
    methods(Static)
        function seOut = ComputePETH(se, ops)
            % Compute PETHs for each unit in each se object
            % 
            %   seOut = ComputePETH(se)
            %   seOut = ComputePETH(se, ops)
            % 
            % Inputs
            %   se          MSessionExplorer object(s) containing the following data:
            %                   se.userData.ksMeta.clusTb
            %                   'spikeRate' or 'spikeTime' table (will use 'spikeRate' table if both are present)
            %   ops         Option struct(s) containing the following fields:
            %                   rsWin           Time window to compute PETH, in a two-element vector.
            %                   rsBinSize       PETH bin size in second.
            %               If ops is not provided or left empty, the function will find it at se.userData(1).ops. 
            %               If only one struct is provided but se includes more than one object, the same struct 
            %               will be applied to all se objects.
            % Output
            %   seOut       MSessionExplorer object(s) where each epoch is computed from an input se object.
            %               It contains the following tables:
            %                   'mean'          Timeseries of PETH means.
            %                   'sem'           Timeseries of standard deviation of the means.
            %                   'taskValue'     Task attributes of the first or the template (if available) trial.
            %               And seOut.userData.clusTb
            %                   It is a copy of se.userData.keMeta.clusTb with the following additional columns:
            %                       peakSpkCount, peakSpkRate, peakSpkProb
            %                   If seOut has multiple epochs (i.e. PETHs from multiple se objects), these columns 
            %                   store the peak values across all PETHs.
            % 
            
            if nargin < 2 || isempty(ops)
                ops = arrayfun(@(x) x.userData(1).ops, se);
            end
            
            if numel(se) > 1
                % Compute recursively
                seArray = se;
                if numel(ops) == 1
                    ops = repmat(ops, size(seArray));
                end
                seArray = arrayfun(@(x,y) NP.Unit.ComputePETH(x,y), seArray, ops);
                
                % Find peak stats from all results
                for i = numel(seArray) : -1 : 1
                    cTb = seArray(i).userData.clusTb;
                    peakSpkCount(:,i) = cTb.peakSpkCount;
                    peakSpkRate(:,i) = cTb.peakSpkRate;
                    peakSpkProb(:,i) = cTb.peakSpkProb;
                end
                cTb.peakSpkCount = max(peakSpkCount, [], 2);
                cTb.peakSpkRate = max(peakSpkRate, [], 2);
                cTb.peakSpkProb = max(peakSpkProb, [], 2);
                
                % Merge results
                seOut = Merge(seArray);
                seOut.userData = [];
                seOut.userData.clusTb = cTb;
                return
            end
            
            fprintf("Compute PETHs from %s\n", NP.SE.GetID(se));
            
            % Compute histograms
            tWin = ops.rsWin;
            tEdges = (tWin(1) : ops.rsBinSize : tWin(2))';
            if ismember('spikeRate', se.tableNames)
                respTb = se.SliceTimeSeries('spikeRate', tWin, 'Fill', 'none');
                respTb = se.ResampleTimeSeries(respTb, tEdges);
                [mm, ee, stats] = MNeuro.MeanTimeSeries(respTb{:,2:end});
            else
                respTb = se.SliceEventTimes('spikeTime', tWin, 'Fill', 'none');
                [mm, ee, stats] = MNeuro.MeanEventRate(respTb, tEdges);
            end
            
            % Organize PETH timeseries into tables
            rowDist = size(mm,1);
            colDist = ones(1, size(mm,2));
            mm = mat2cell(mm, rowDist, colDist);
            ee = mat2cell(ee, rowDist, colDist);
            
            tCenters = MMath.BinEdges2Centers(tEdges);
            clusIds = respTb.Properties.VariableNames;
            mm = MSessionExplorer.MakeTimeSeriesTable(tCenters, mm, 'VariableNames', clusIds, 'Verbose', false);
            ee = MSessionExplorer.MakeTimeSeriesTable(tCenters, ee, 'VariableNames', clusIds, 'Verbose', false);
            
            % Add stats to clusTb
            cTb = NP.Unit.GetClusTb(se);
            cTb.peakSpkCount = stats.pkVal * ops.rsBinSize;
            cTb.peakSpkRate = stats.pkVal;
            cTb.peakSpkProb = stats.pkProb;
            
            % Get task attributes from taskValue table
            tv = se.GetTable('taskValue');
            k = [];
            if ismember('tempTrialNum', tv.Properties.VariableNames)
                k = find(tv.trialNum == tv.tempTrialNum(1), 1);
            end
            if isempty(k)
                k = 1;
            end
            tv = tv(k,:);
            
            % Construct se for output
            seOut = MSessionExplorer();
            seOut.SetTable('mean', mm, 'timeSeries');
            seOut.SetTable('sem', ee, 'timeSeries');
            seOut.SetTable('taskValue', tv, 'eventValues');
            seOut.userData.clusTb = cTb;
        end
        
        function mi = ComputeModulationIndex(mm, ee)
            % Compute error-discounted modulation index
            %
            %   mi = ComputeModulationIndex(mm)
            %   mi = ComputeModulationIndex(mm, ee)
            % 
            % Inputs
            %   mm          1) A t-by-r matrix of mean responses where t is the number of time points and 
            %                  r is the number of responses.
            %               2) A cell array of mean response vectors.
            %   ee          Error of mm, with matching dimensions. Will use zero error if not provided.
            % Output
            %   mi          A numeric matrix of modulation indices.
            % 
            
            if nargin < 2
                ee = mm;
            end
            
            if istable(mm)
                mm = mm{:,:};
                ee = ee{:,:};
            end
            
            if isnumeric(mm)
                if isvector(mm)
                    mm = mm(:);
                    ee = ee(:);
                end
                rowDist = size(mm,1);
                colDist = ones(size(mm,2), 1);
                mm = mat2cell(mm, rowDist, colDist);
                ee = mat2cell(ee, rowDist, colDist);
            end
            
            mi = zeros(size(mm));
            for i = 1 : size(mm,1)
                for j = 1 : size(mm,2)
                    m = mm{i,j};
                    e = ee{i,j};
                    if nargin < 2
                        e(:) = 0;
                    end
                    mPlus = m+e;
                    mMinus = m-e;
                    amp = max(mMinus) - min(mPlus);
                    amp = max(amp, 0); % avoid negative
                    mi(i,j) = amp ./ max(mPlus);
                end
            end
        end
        
        function tb = ComputePatternStats(ce, histInd, isCat)
            % Compute properties stats of PETH patterns
            % 
            %   tb = ComputePatternStats(ce)
            %   tb = ComputePatternStats(ce, histInd)
            %   tb = ComputePatternStats(ce, histInd, isCat)
            % 
            % Inputs
            %   histInd         
            %   isCat           
            % Output
            %   tb              A table with the following variables
            %     pkIdx           Index of the time bin where each trace peaks
            %     pkVal           The value at pkIdx
            %     entropy         Shannon's entropy of each trace
            %     mi              Modulation index
            
            if nargin < 3
                isCat = false;
            end
            if nargin < 2 || isempty(histInd)
               histInd = 1 : ce.numEpochs;
            end
            
            [mTb, eTb] = ce.GetTable('resp', 'sem');
            mTb = mTb(histInd, 2:end);
            eTb = eTb(histInd, 2:end);
            
            if isCat
                mTb = CatRows(mTb);
                eTb = CatRows(eTb);
            end
            
            function tbOut = CatRows(tbIn)
                tbOut = table;
                for c = 1 : width(tbIn)
                    vn = tbIn.Properties.VariableNames{c};
                    tbOut.(vn) = {cat(1, tbIn.(vn){:})};
                end
            end
            
            tb = table;
            for i = width(mTb) : -1 : 1
                for j = height(mTb) : -1 : 1
                    % Get a PETH and errors
                    m = mTb.(i){j};
                    e = eTb.(i){j};
                    
                    % Compute peak info
                    [tb.pkVal(i,j), tb.pkIdx(i,j)] = max(m);
                    
                    % Compute entropy
                    tb.entropy(i,j) = MMath.Entropy(m/sum(m));
                    
                    % Compute modulation index
                    rel = max(m-e, 0);
                    tb.mi(i,j) = (max(rel) - min(rel)) ./ max([m; NP.Param.normAddMax]);
                end
            end
        end
        
        function [sSpk, clusTb] = ExtractPeriEventSpikes(se, ops)
            % Extract peri-event spike times and compute PETHs
            % 
            %   [sSpk, clusTb] = ExtractPeriEventSpikes(se, ops)
            % 
            
            sSpk = struct;
            
            for i = 1 : numel(ops.rsVars)
                % Unpack struct
                s = ops.rsVars(i);
                tn = s.tableName;
                vn = s.varNames;
                rd = s.readout;
                if isempty(rd)
                    rd = "default";
                end
                fprintf("Extract spike times using %s from the '%s' table using %s readout.\n", strjoin(vn, ', '), tn, rd);
                
                % Get trigger table
                k = strcmp(tn, se.tableNames);
                if ~se.isEventTimesTable(k) && ~se.isTimesSeriesTable(k)
                    fprintf("The requested table '%s' is not an applicable tyep.\n", tn);
                    continue
                end
                trigTb = se.GetTable(tn);
                
                % Iterate through each variable of interest in the trigger table
                for j = 1 : numel(vn)
                    % Get spikeTime table
                    stTb = se.GetTable('spikeTime');
                    vnUni = vn{j};
                    if ismember(vnUni, fieldnames(sSpk))
                        vnUni = vnUni + "_" + tn;
                    end
                    
                    % Make sure epoch data is of type double and in cell array
                    trig = trigTb.(vn{j});
                    if ~iscell(trig)
                        trig = num2cell(trig);
                    end
                    trig = cellfun(@double, trig, 'Uni', false);
                    
                    % Convert trigger timeseries to event times
                    if se.isTimesSeriesTable(k)
                        isShort = cellfun(@numel, trig) < 3; % findpeaks requires at least 3 date points
                        trig(isShort) = {zeros(3,1)}; % this will make findpeaks return empty event time
                        trigTb.time(isShort) = {(1:3)'};
                        [~, trig] = cellfun(@(x,t) findpeaks(x,t), trig, trigTb.time, 'Uni', false);
                    end
                    
                    % Convert trigger times to slicing windows
                    tWins = cellfun(@(x) x + ops.rsWin, trig, 'Uni', false);
                    
                    % Remove empty epochs
                    isNone = cellfun(@(x) isempty(x) || all(isnan(x)), trig);
                    trig(isNone) = [];
                    tWins(isNone) = [];
                    stTb(isNone,:) = [];
                    if all(isNone)
                        stTb = array2table(num2cell(NaN(1, width(stTb))), 'VariableNames', stTb.Properties.VariableNames);
                        sSpk.(vnUni) = stTb;
                        continue
                    end
                    
                    % Slicing
                    stTb = se.SliceEventTimes(stTb, tWins, [], ops.unitInd);
                    
                    % Reference to trigger times
                    trig = cat(1, trig{:});
%                     trig(isnan(trig)) = 0; % this line shouldn't be necessary
                    seTemp = MSessionExplorer;
                    seTemp.SetTable('spikeTime', stTb, 'eventTimes');
                    seTemp.AlignTime(trig);
                    seTemp.Column2Cell('spikeTime');
                    stTb = seTemp.GetTable('spikeTime');
                    
                    sSpk.(vnUni) = stTb;
                end
            end
            
            % Compute PETHs
            fprintf("Compute PETHs.\n");
            tEdges = (ops.rsWin(1) : 0.025 : ops.rsWin(2))';
            fns = fieldnames(sSpk);
            peth = cell(numel(fns), 3);
            for e = 1 : numel(fns)
                st = sSpk.(fns{e}){:,:};
                [rm, re] = MNeuro.MeanEventRate(st, tEdges);
                peth{e,1} = rm;
                peth{e,2} = re;
                peth{e,3} = size(st, 1);
            end
            
            rMax = max(cat(1,peth{:,1}), [], 1);
            peth(:,1:2) = cellfun(@(x) x./rMax, peth(:,1:2), 'Uni', false);
            
            tCenters = MMath.BinEdges2Centers(tEdges);
            units = stTb.Properties.VariableNames;
            sPeth = cell(numel(units), numel(fns));
            for u = 1 : numel(units)
                for e = 1 : numel(fns)
                    s = struct;
                    s.time = tCenters;
                    s.mean = peth{e,1}(:,u);
                    s.sem = peth{e,2}(:,u);
                    s.nEvents = peth{e,3};
                    s.scale = rMax(u);
                    sPeth{u,e} = s;
                end
            end
            pethTb = cell2table(sPeth, 'VariableNames', fns);
            
            clusTb = NP.Unit.GetClusTb(se);
            if ~isempty(ops.unitInd)
                clusTb = clusTb(ops.unitInd,:);
            end
            clusTb = [clusTb pethTb];
        end
        
        function PlotFeatureResponse(sSpk, pethTb, uInd, varargin)
            % Plot peri-event spike rasters and histograms
            % 
            %   PlotFeatureResponse(sSpk, pethTb, uInd)
            % 
            
            p = inputParser;
            p.addParameter('FeatureNames', fieldnames(sSpk), @(x) isstring(x) || iscellstr(x) || ischar(x));
            p.addParameter('Raster', true, @isscalar);
            p.addParameter('PETH', true, @isscalar);
            p.parse(varargin{:});
            featNames = p.Results.FeatureNames;
            
            nFeat = numel(featNames);
            nRow = 7;
            nCol = ceil(nFeat/(nRow-1));
            tl = tiledlayout(nRow, nCol);
            tl.Padding = 'compact';
            
            for u = 1 : numel(uInd)
                ax = nexttile;
                NP.UnitPlot.Waveform(ax, pethTb, uInd(u));
                ax.Title.String = sprintf('u%i@%ium', pethTb.clusId(uInd(u)), pethTb.depth(uInd(u)));
            end
            for i = (numel(uInd)+1) : nCol
                ax = nexttile;
                ax.Visible = 'off';
            end
            
            for i = 1 : numel(featNames)
                fn = featNames{i};
                st = sSpk.(fn){:,uInd};
                nEvts = size(st,1);
                rsInd = unique(round(linspace(1, nEvts, 50)));
                st = st(rsInd,:);
                
%                 stClean = cell(1, size(st,2));
%                 isRm = cellfun(@(x) isempty(x) || all(isnan(x)), st);
%                 for j = 1 : numel(stClean)
%                     stUnit = st(:,j);
%                     stUnit(isRm(:,j)) = [];
%                     stClean{j} = stUnit;
%                 end
                
                sPeth = pethTb.(fn);
                t = sPeth(uInd(1)).time;
                hh = cat(2, sPeth(uInd).mean);
                ee = cat(2, sPeth(uInd).sem);
                mi = NP.PETH.ComputeModulationIndex(hh, ee);
                
                ax = nexttile;
                if numel(t) == size(hh,1)
                    MPlot.PlotHistStack(t, hh, ee);
                    MPlot.PlotRasterStack(st);
                end
                plot(ax, [0 0]', ax.YLim, 'Color', [0 0 0 .3]);
                ax.XLim = t([1 end]);
                ax.XGrid = 'off';
                ax.Box = 'off';
                ax.Title.String = sprintf("%s (n=%i)", featNames{i}, nEvts);
                ax.YTickLabel = "MI:"+round(mi,2);
                ax.YTickLabelRotation = 90;
            end
        end
        
    end
    
end