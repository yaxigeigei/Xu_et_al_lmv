classdef Unit
    
    methods(Static)
        % Response
        function RemoveUnits(se, uInd)
            % Remove units from se
            % 
            %   RemoveUnits(se, uInd)
            % 
            
            if islogical(uInd)
                uInd = find(uInd);
            end
            
            se.userData.ksMeta.clusTb(uInd,:) = [];
            
            for i = 1 : numel(se.tableNames)
                tn = se.tableNames{i};
                if ~startsWith(tn, 'spike')
                    continue
                end
                tb = se.GetTable(tn);
                if se.isTimesSeriesTable(i)
                    tb(:,uInd+1) = [];
                elseif se.isEventTimesTable(i)
                    tb(:,uInd) = [];
                else
                    error("'%s' has invalid table type.", tn);
                end
                se.SetTable(tn, tb);
            end
        end
        
        function AddSpikeRateTable(se, ops)
            % Convert spike times to spike rates and add to se
            %
            %   AddSpikeRateTable(se, ops)
            % 
            
            % Vectorize spike times
            seLite = se.Duplicate({'spikeTime'}, false);
            seLite.SliceSession(0, 'absolute');
            
            % Clean ISI violated spikes
            spk = seLite.GetTable('spikeTime');
            spk = NP.Unit.CleanSpikes(spk, NP.Param.RP);
            
            % Add time lag
            for i = 1 : width(spk)
                if isnumeric(spk.(i))
                    spk.(i) = {spk.(i)};
                end
                spk.(i){1} = spk.(i){1} + ops.spkLagInSec(1);
            end
            seLite.SetTable('spikeTime', spk);
            
            % Binning
            tRef = se.GetReferenceTime('spikeTime');
            if ismember('ni', se.tableNames)
                tNI = se.GetColumn('ni', 'time');
                tWin = [tNI{1}(1) tNI{end}(end)] + tRef([1 end]')';
            else
                ksOps = se.userData.ksMeta.ops;
                tWin = ksOps.trange * ksOps.fs / se.userData.apMeta.imSampRate;
            end
            tEdges = tWin(1) : ops.spkBinSize : tWin(2);
            r = seLite.ResampleEventTimes('spikeTime', tEdges, 'Normalization', 'countdensity');
            
            % Smoothing
            if ops.spkKerSize
                for i = 2 : width(r)
                    r.(i){1} = MNeuro.Filter1(r.(i){1}, 1/ops.spkBinSize, 'gaussian', ops.spkKerSize);
                end
            end
            seLite.SetTable('spikeRate', r, 'timeSeries', 0);
            seLite.RemoveTable('spikeTime');
            
            % Reslice
            tSlice = tRef;
            tSlice(1) = 0;
            seLite.SliceSession(tSlice, 'absolute');
            seLite.AlignTime(tRef-tSlice);
            r = seLite.GetTable('spikeRate');
            se.SetTable('spikeRate', r, 'timeSeries', tRef);
        end
        
        function AddSimSpikeTimeTable(se, tbName)
            % Add a table of spike time simulated or reconstructed from spike rates
            %
            %   AddSimSpikeTimeTable(se)
            %   AddSimSpikeTimeTable(se, tbName)
            % 
            if ~exist('tbName', 'var') || isempty(tbName)
                tbName = 'spikeTimeSim';
            end
            sr = se.GetTable('spikeRate');
            rt = se.GetReferenceTime('spikeRate');
            rng(61);
            st = NP.Unit.SimSpiking(sr);
            se.SetTable(tbName, st, 'eventTimes', rt);
        end
        
        function st = SimSpiking(sr)
            % Simulate spike times from spike rates
            % 
            %   spikeTimeTable = SimSpiking(spikeRateTable)
            %
            st = sr(:,2:end);
            rsF = 10e3;
            for i = 1 : height(sr)
                t = sr.time{i};
                tRS = (t(1) : 1/rsF : t(end))';
                for j = 2 : width(sr)
                    r = sr.(j){i};
                    if ~any(r > 0) % for speeding up
                        st.(j-1){i} = NaN;
                        continue
                    end
                    rRS = interp1(t, r, tRS, 'linear');
                    p = rRS / rsF; % convert spikes per sec to spk prob per time bin
                    rn = rand(size(rRS));
                    isSpk = rn < p;
                    spk = find(isSpk);
                    if isempty(spk)
                        spk = NaN;
                    else
                        spk = tRS(spk);
                    end
                    st.(j-1){i} = spk;
                end
            end
            st = NP.Unit.CleanSpikes(st, NP.Param.RP);
        end
        
        function DetrendSpikeRates(se, varargin)
            % Detrend spike rate for each unit
            % 
            %   DetrendSpikeRates(se)
            % 
            
            rt = se.GetReferenceTime("spikeRate");
            
            if ismember("spikeRate0", se.tableNames)
                srTb = se.GetTable("spikeRate0");
            else
                srTb = se.GetTable("spikeRate");
                se.SetTable("spikeRate0", srTb, 'eventTimes', rt);
            end
            
            t = cellfun(@(x,y) x+y, srTb.time, num2cell(rt), 'Uni', false);
            t = cat(1, t{:});
            nSp = cellfun(@numel, srTb.time);
            for i = 2 : width(srTb)
                r0 = cat(1, srTb.(i){:});
                m = r0 > 0;
                mdl = fitlm(t(m), r0(m));
                r = r0 ./ predict(mdl, t);
                r = r / sum(r, "omitmissing") * sum(r0, "omitmissing");
                srTb.(i) = mat2cell(r, nSp);
            end
            
            se.SetTable("spikeRate", srTb);
        end
        
        % Quality control
        function clusTb = ComputeQualityTable(se, taskNames)
            % Return cluster table with updated quality metrics
            % 
            %   clusTb = ComputeQualityTable(se)
            %   clusTb = ComputeQualityTable(se, taskNames)
            % 
            
            % Start with the existing cluster table
            NP.Unit.MergeMeta(se);
            clusTb = NP.Unit.GetClusTb(se);
            clusTb = NP.Unit.AddRecMeta(se, clusTb);
            
            % Correct contamination stats
            clusTb = NP.Unit.NoDupContamStats(clusTb);
            
            % Add span metric
            if ~exist('taskNames', 'var') || isempty(taskNames)
                clusTb.fracSpan = NP.Unit.ComputeTaskSpan(se);
            else
                clusTb.fracSpan = NP.Unit.ComputeTaskSpan(se, taskNames);
            end
            
            % Add quality labels
            clusTb.isRPV = clusTb.RPV_ND > NP.Param.maxRPV;
            clusTb.isContam = clusTb.contamND > NP.Param.maxContam;
            clusTb.isShort = clusTb.fracSpan < NP.Param.minFracSpan;
            clusTb.isMulti = clusTb.isRPV | clusTb.isContam;
            clusTb.isSingle = ~clusTb.isMulti;
        end
        
        function RemoveDuplicateSpikes(se)
            % Remove spikes following an ISI less than the threshold defined by NP.Param.minISI.
            % This is a wrapper of NP.Unit.CleanSpikes that takes se object as input.
            % 
            %   NP.Unit.RemoveDuplicateSpikes(se)
            % 
            
            st = se.GetTable('spikeTime');
            st = NP.Unit.CleanSpikes(st, NP.Param.minISI);
            se.SetTable('spikeTime', st);
        end
        
        function RemoveRPVSpikes(se)
            % Remove spikes following an ISI less than the refractory period defined by NP.Param.RP.
            % This is a wrapper of NP.Unit.CleanSpikes that takes se object as input.
            %
            %   NP.Unit.RemoveRPVSpikes(se)
            %
            st = se.GetTable('spikeTime');
            st = NP.Unit.CleanSpikes(st, NP.Param.RP);
            se.SetTable('spikeTime', st);
        end
        
        function [spkTb, vioTb] = CleanSpikes(spkTb, isiLimit)
            % Remove the spikes following ISI less than the specified limit in sec
            % This can be used to fix double-counting, e.g. isiLimit = NP.Param.minISI (0.5e-3)
            % or to clean RPVs, e.g. isi_limit = NP.Param.RP (1.5e-3)
            %
            %   [spkTb, vioTb] = NP.Unit.CleanSpikes(spkTb, isiLimit)
            %
            
            spkCell = table2cell(spkTb);
            vioCell = cell(size(spkCell));
            
            for i = 1:numel(spkCell)
                isVio = diff(spkCell{i}) < isiLimit;
                isVio = [false; isVio(:)];
                if any(isVio)
                    vioCell{i} = spkCell{i}(isVio);
                    spkCell{i} = spkCell{i}(~isVio);
                end
            end
            
            spkTb = cell2table(spkCell, 'VariableNames', spkTb.Properties.VariableNames);
            vioTb = cell2table(vioCell, 'VariableNames', spkTb.Properties.VariableNames);
            
            spkTb = col2cell(spkTb);
            vioTb = col2cell(vioTb);
            
            function tb = col2cell(tb)
                for k = 1 : width(tb)
                    C = tb.(k);
                    if ~iscell(C)
                        C = num2cell(C);
                    end
                    tb.(k) = C;
                end
            end
        end
        
        function StandardizeSpikeTimeTable(se)
            % Make sure spike time table is in cell array
            tb = se.GetTable("spikeTime");
            for k = 1 : width(tb)
                C = tb.(k);
                if ~iscell(C)
                    C = num2cell(C);
                end
                tb.(k) = C;
            end
            se.SetTable("spikeTime", tb);
        end
        
        function s = ComputeContamStats(tSpk)
            % Compute unit quality metrics about contamination
            % 
            %   s = ComputeContamStats(tSpk)
            % 
            % Input
            %   tSpk                A vector of spike times.
            % 
            % Output
            %   Struct s with the following fields
            %   s.isiEdges          Bin edges of inter-spike interval (ISI) histogram in seconds.
            %   s.isiCount          Spike count of ISI histogram.
            %   s.RPV               Refractory period violation rate (%). The RP threshold is NP.Param.RP.
            %   s.meanActiveRate    Mean spike rate when unit is active.
            %   s.activeThreshold   The threshold for active firing, computed from NP.Param.activePrct.
            %   s.contam            Contamination rate (%) (Hill, Mehta and Kleinfeld, J Neuro, 2011).
            %
            
            % Compute inter-spike intervals
            tSpk = unique(tSpk); % remove merging artifacts - Kilosort or Phy's bug?
            isi = diff(tSpk);
            
            % ISI histogram and the rate of refractory period violation (RPV)
            tEdges = 0 : 0.5e-3 : 0.02;
            nISI = histcounts(isi, [tEdges Inf]);
            RPV = sum(nISI(tEdges < NP.Param.RP)) / numel(isi);
            s.isiEdges = tEdges;
            s.isiCount = nISI(1:end-1);
            s.RPV = RPV * 100;
            
            % Contamination rate
            r = 1 ./ isi;
            th = 1 / prctile(isi, NP.Param.activePrct);
            isActive = r > th;
            rMeanActive = sum(isActive) / sum(isi(isActive));
            c = MNeuro.ClusterContamination(RPV, rMeanActive, NP.Param.RP);
            s.meanActiveRate = rMeanActive;
            s.activeThreshold = th;
            s.contam = c * 100;
        end
        
        function clusTb = NoDupContamStats(clusTb)
            % Estimate RPV and contam rates without duplicated spikes
            % 
            %   clusTb = NoDupContamStats(clusTb)
            % 
            % Inputs require the following columns in clusTb
            %   isiEdges, isiCount, rMeanActive
            % 
            % Outputs add to the following columns in clusTb
            %   isiCountND, RPV_ND, contamND
            % 
            % See also NP.Unit.ComputeContamStats
            
            for i = 1 : height(clusTb)
                % Replace values in duplication bins with the mean of non-dup RPV bins
                tEdges = clusTb.isiEdges{i};
                nISI = clusTb.isiCount{i};
                isDup = tEdges < NP.Param.minISI; % mask of bins for duplications and the subset of RPVs
                isVio = tEdges < NP.Param.RP; % mask of all RPVs bins
                nISI(isDup) = mean(tEdges(isVio & ~isDup));
                clusTb.isiCountND{i} = nISI;
                
                % Compute the rate of refractory period violation (RPV)
                d = sum(nISI) - sum(clusTb.isiCount{i});
                RPV = sum(nISI(isVio)) / (clusTb.numSpikes(i) - 1 + d);
                clusTb.RPV_ND(i) = RPV * 100;
                
                % Contamination rate
                c = MNeuro.ClusterContamination(RPV, clusTb.meanActiveRate(i), NP.Param.RP);
                clusTb.contamND(i) = c * 100;
            end
        end
        
        function [clusTb, isSingle] = FindSingleUnits(clusTb, varargin)
            % Return clusTb with only single units.
            % 
            %   [clusTb, isSingle] = FindSingleUnits(clusTb, varargin)
            % 
            % Inputs
            %   clusTb              aaa
            %   'MaxRPV'            aaa
            %   'MaxContam'         aaa
            %   'MinNumSpikes'      aaa
            %   'MinSNR'            aaa
            % Outputs
            %   clusTb              aaa
            %   isSingle            aaa
            %   
            
            p = inputParser();
            p.addParameter('MaxRPV', NP.Param.maxRPV, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('MaxContam', NP.Param.maxContam, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('MinNumSpikes', NP.Param.minNumSpk, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('MinSNR', NP.Param.minSNR, @(x) isscalar(x) && isnumeric(x));
            p.parse(varargin{:});
            
            isSingle = ...
                clusTb.RPV < p.Results.MaxRPV & ...
                clusTb.contam < p.Results.MaxContam & ...
                clusTb.numSpikes > p.Results.MinNumSpikes & ...
                clusTb.SNR > p.Results.MinSNR;
            
            clusTb = clusTb(isSingle,:);
        end
        
        function [fracSpan, taskNames] = ComputeTaskSpan(se, taskNames, capTrialDur)
            % Compute the fraction of task duration that a unit is present
            % 
            %   [fracSpan, taskNames] = ComputeTaskSpan(se)
            %   [fracSpan, taskNames] = ComputeTaskSpan(se, taskNames)
            %   [fracSpan, taskNames] = ComputeTaskSpan(se, taskNames, capTrialDur)
            % 
            % Inputs
            %   se              An MSessionExplorer object with trials as epochs.
            %   taskNames       One or more task names to compute span for. The function finds all unique tasks 
            %                   from the taskValue table if the input is empty [].
            %   capTrialDur     The maximal duration to be included in computing span. Default is 10 seconds.
            %                   This helps prevent underestimation of span when there are large gaps in the 
            %                   task and/or extra waiting time after the task ends.
            % Output
            %   fracSpan        A unit-by-task array indicating the number of rounds a unit appeared for each task.
            %   taskNames       Same as the input taskNames. Useful when the input was not provided.
            % 
            
            if ~ismember('spikeSpan', se.tableNames)
                NP.Unit.AddSpikeSpanTable(se);
            end
            
            rt = se.GetReferenceTime('spikeSpan');
            dur = [diff(rt); Inf];
            tWin = zeros(numel(rt), 2);
            if nargin < 3 || isempty(capTrialDur)
                capTrialDur = 10;
            end
            tWin(:,2) = min(dur, capTrialDur);
            
            ssTb = se.SliceTimeSeries('spikeSpan', tWin);
            
            if ~ismember('taskValue', se.tableNames)
                % Compute span over the full recording
                taskNames = [];
                fracSpan = mean(cell2mat(ssTb{:,2:end}), 1)';
                return
            end
            
            % Find tasks to compute for
            tv = se.GetTable('taskValue');
            if nargin < 2 || isempty(taskNames)
                taskNames = unique(tv.taskName, 'stable');
            end
            taskNames = lower(string(taskNames));
            
            % Compute span for each task
            for i = numel(taskNames) : -1 : 1
                % Find task epochs
                isTask = strcmpi(tv.taskName, taskNames(i));
                if ~any(isTask)
                    error("Cannot compute unit span for '%s' since the task is not present in this recording.", taskNames(i));
                end
                
                % Compute fraction of time
                fracSpan(:,i) = mean(cell2mat(ssTb{isTask,2:end}), 1)';
                
                % Convert fraction of time to fraction of target task length
                switch taskNames(i)
                    case "lmv"
                        n = 72;
                    case "timit"
                        n = 200;
                    otherwise
                        n = sum(isTask);
                end
                fracSpan(:,i) = fracSpan(:,i) * sum(isTask) / n;
            end
        end
        
        function AddSpikeSpanTable(se, binSize)
            % Compute 'spikeSpan' timeseries table and add it to se
            %
            %   AddSpikeSpanTable(se)
            %   AddSpikeSpanTable(se, binSize)
            % 
            
            % Vectorize spike times
            seLite = se.Duplicate({'spikeTime'}, false);
            seLite.SliceSession(0, 'absolute');
            
            % Clean ISI violated spikes
            spk = seLite.GetTable('spikeTime');
            spk = NP.Unit.CleanSpikes(spk, NP.Param.RP);
            
            % Determine time bins
            tRef = se.GetReferenceTime('spikeTime');
            if ismember('ni', se.tableNames)
                tNI = se.GetColumn('ni', 'time');
                tWin = [tNI{1}(1) tNI{end}(end)] + tRef([1 end]')';
            else
                ksOps = se.userData.ksMeta.ops;
                tWin = ksOps.trange * ksOps.fs / se.userData.apMeta(1).imSampRate;
            end
            if nargin < 2
                binSize = 0.1; % default resolution 0.1 sec
            end
            tEdges = tWin(1) : binSize : tWin(2);
            
            % Find temporal spans
            ssTb = table;
            ssTb.time{1} = MMath.BinEdges2Centers(tEdges);
            for i = 1 : width(spk)
                uName = spk.Properties.VariableNames{i};
                try
                    ssTb.(uName){1} = NP.Unit.FindTemporalSpan(spk.(uName){1}, tEdges);
                catch
                    disp("?");
                end
            end
            seLite.SetTable('spikeSpan', ssTb, 'timeSeries', 0);
            seLite.RemoveTable('spikeTime');
            
            % Reslice
            tSlice = tRef;
            tSlice(1) = 0;
            seLite.SliceSession(tSlice, 'absolute');
            seLite.AlignTime(tRef-tSlice);
            ssTb = seLite.GetTable('spikeSpan');
            se.SetTable('spikeSpan', ssTb, 'timeSeries', tRef);
        end
        
        function [isActive, tBin] = FindTemporalSpan(tSpk, tEdges)
            % Find periods where interspike interval (ISI) is shorter than NP.Param.activePrct percentile of all ISIs.
            % This measures whether a unit is present or minimally active in a recording.
            %
            %   [isActive, tBin] = FindTemporalSpan(tSpk, tEdges)
            % 
            
            % Compute inter-spike intervals
            tSpk = unique(tSpk); % remove duplicates and sort
            isi = diff(tSpk);
            
            % Thresholding
            tBin = MMath.BinEdges2Centers(tEdges(:));
            hBin = diff(tEdges(1:2)) / 2;
            th = prctile(isi, NP.Param.activePrct);
            isActive = false(size(tBin));
            for i = 1 : numel(isi)
                if isi(i) < th
                    mask = tBin > (tSpk(i)-hBin) & tBin < (tSpk(i+1)+hBin); % extend half bins to make sure every ISI passing threshold contributes
                    isActive(mask) = true;
                end
            end
            
            % Fill gaps less than 120 seconds (roughly 18 trials or one round)
            isInactive = ~isActive;
            isInactive([1 end]) = true; % capture the trimming
            gapInd = MMath.Logical2Bounds(isInactive);
            gapDur = diff(tBin(gapInd), 1, 2);
            for i = 2 : numel(gapDur)-1 % exclude the beginning and the end
                if gapDur(i) < 120
                    isActive(gapInd(i,1):gapInd(i,2)) = true;
                end
            end
        end
        
        % Utilities
        function MergeMeta(se)
            % Merge clusTb from multiple ksMeta structs
            % 
            %   MergeMeta(se)
            % 
            
            % Get metadata struct
            fprintf("\nUse metadata from se.usreData.ksMeta.\n");
            meta = se.userData.ksMeta;
            
            if isscalar(meta)
                fprintf("Metadata contains only one struct. No need to merge.\n");
                return
            end
            
            % Back up original metadata
            fprintf("Backup saved at se.userData.ksMeta0.\n");
            se.userData.ksMeta0 = meta;
            
            % Concatenate clusTbs
            clusTbCat = cat(1, meta.clusTb);
            meta = meta(1);
            meta.clusTb = clusTbCat;
            se.userData.ksMeta = meta;
        end
        
        function [uId, uIdUni] = GetSelectedClusId(recId)
            % Return handpicked cluster IDs for a given recording
            % 
            %   [uId, uIdUni] = GetSelectedClusId(recId)
            % 
            uId{1} = LMV.Param.GetSelectedClusId(recId);
            uId{2} = Semsr.Param.GetSelectedClusId(recId);
            uId = cat(2, uId{:});
            uId = unique(uId);
            uIdUni = uId + NP.Unit.GetBaseClusId(recId);
        end
        
        function SetUniqueClusId(se)
            % Make cluster IDs unique across experiments and subjects
            % New ID has 9 digits, SSBBPUUUU, where SS is subject number, BB is block number,
            % P is probe or imec number (added in preprocessing), UUUU is Kilosort cluster ID.
            % 
            %   SetUniqueClusId(se)
            % 
            % See also NP.Unit.GetRecBaseId
            
            id = se.userData.ksMeta.clusTb.clusId(1);
            if id > 1e7
                fprintf("\nThis se already uses unique clusters IDs, e.g. the first cluster is %i.\n", id);
                return
            end
            
            % Get the base ID
            baseId = NP.Unit.GetBaseClusId(se);
            fprintf("\nCluster IDs will be incremented by %i\n", baseId);
            
            % Change uId in metadata
            ksMeta = se.userData.ksMeta;
            ksMeta.clusTb.clusId = ksMeta.clusTb.clusId + baseId;
            ksMeta.spkTb.clusId = ksMeta.spkTb.clusId + baseId;
            se.userData.ksMeta = ksMeta;
            fprintf("Changed cluster IDs in ksMeta\n");
            
            % change uId in table variable names
            tbNames = intersect({'spikeTime', 'spikeRate', 'spikeSpan', 'resp'}, se.tableNames);
            for i = 1 : numel(tbNames)
                tn = tbNames{i};
                tb = se.GetTable(tn);
                vn = tb.Properties.VariableNames;
                m = startsWith(vn, 'u');
                vn(m) = cellstr(cellfun(@(x) "u"+(str2double(x(2:end))+baseId), vn(m)));
                tb.Properties.VariableNames = vn;
                se.SetTable(tn, tb);
                fprintf("Changed cluster IDs in '%s' table\n", tn);
            end
        end
        
        function baseId = GetBaseClusId(src)
            % Get the base that makes cluster IDs unique across experiments and subjects
            % New ID has 9 digits, SSBBPUUUU, where SS is subject number, BB is block number,
            % P is probe or imec number (added in preprocessing), UUUU is Kilosort cluster ID.
            % Therefore, the increment or base is SSBBP0000.
            % 
            %   baseId = GetBaseClusId(src)
            % 
            % See also NP.Unit.SetUniqueClusId
            % 
            if isnumeric(src)
                src = num2str(src);
                subjNum = str2double(src(1:end-7));
                blkNum = str2double(src(end-6:end-5));
            else
                [~, subjId, blkId] = NP.SE.GetID(src);
                subjNum = str2double(erase(subjId, 'NP'));
                blkNum = str2double(erase(blkId, 'B'));
            end
            baseId = subjNum*1e7 + blkNum*1e5;
        end
        
        function probeId = GetProbeId(clusId)
            % Get probe ID (imec number) from unique cluster IDs
            % 
            %   probeId = GetProbeId(clusId)
            % 
            clusId = cellstr(string(clusId));
            probeId = cellfun(@(x) x(end-4), clusId, 'Uni', false);
            probeId = str2double(probeId);
        end
        
        function clusId = ToUniqueClusId(clusId, src)
            % Convert raw clusId to unique clusId
            % 
            %   clusId = ToUniqueClusId(clusId, src)
            % 
            
            % Convert name IDs to numeric IDs
            isText = isstring(clusId) || ischar(clusId) || iscellstr(clusId);
            if isText
                clusId = erase(clusId, 'u');
                clusId = str2double(clusId);
            end
            
            % Check ID type
            isRaw = clusId < 1e6;
            if any(~isRaw)
                fprintf("\n%i clusId are already unique, e.g. %i.\n", sum(~isRaw), clusId(find(~isRaw,1)));
            end
            
            % Conversion
            if ischar(src) || iscellstr(src)
                src = string(src);
            end
            if isscalar(src)
                clusId(isRaw) = clusId(isRaw) + NP.Unit.GetBaseClusId(src);
            else
                clusId = clusId + arrayfun(@NP.Unit.GetBaseClusId, src) .* isRaw;
            end
            
            % Convert back to name IDs
            if isText
                clusId = "u" + clusId;
            end
        end
        
        function numId = GetAllClusId(src)
            % Get all the cluster IDs from the source data structure
            % src can be se, userData, ksMeta, a spike table, or anything compatible
            % 
            %   numId = GetAllClusId(se)
            %   numId = GetAllClusId(se.userData)
            %   numId = GetAllClusId(se.userData.ksMeta)
            %   numId = GetAllClusId(se.userData.ksMeta.clusTb)
            %   numId = GetAllClusId(dataTable)
            % 
            
            % Denest or concatenate cell array
            if iscell(src)
                src = cat(1, src{:});
            end
            
            % Get clusTb if src is MSessionExplorer object or metadata struct(s)
            if isa(src, 'MSessionExplorer') || isstruct(src)
                src = NP.Unit.GetClusTb(src);
            end
            
            % Get clusId from table
            colNames = src.Properties.VariableNames;
            if ismember('clusId', colNames)
                % src is a clusTb
                numId = src.clusId;
            else
                % src is a data table such as spikeTime or spikeRate
                nameId = colNames;
                if strcmp(nameId{1}, 'time')
                    nameId = nameId(2:end);
                end
                numId = str2double(erase(nameId, 'u'));
            end
        end
        
        function numId = Ind2ClusId(ind, src)
            % Convert positional indices of units to cluster IDs (in numbers, not strings)
            % 
            %   numId = Ind2ClusId(ind, src)
            % 
            allNumId = NP.Unit.GetAllClusId(src);
            numId = NaN(size(ind));
            m = ~isnan(ind) | ind > numel(allNumId);
            numId(m) = allNumId(ind(m));
        end
        
        function ind = ClusId2Ind(id, src)
            % Convert cluster IDs to positional indices
            % 
            %   ind = ClusId2Ind(id, src)
            % 
            
            % Standardize ID datatype to numbers
            if ischar(id) || iscellstr(id) || isstring(id)
                nameId = string(cellstr(id));
                numId = str2double(erase(nameId, 'u'));
            else
                numId = id;
            end
            
            % Find indices
            allNumId = NP.Unit.GetAllClusId(src);
            ind = NaN(size(numId));
            for k = 1 : numel(numId)
                idx = find(allNumId == numId(k), 1);
                if ~isempty(idx)
                    ind(k) = idx;
                end
            end
        end
        
        function clusTb = GetClusTb(src)
            % Get concatenated clusTb from all ksMeta in the source data structure
            % src can be se, userData, ksMeta, a spike table, or anything compatible
            % 
            %   clusTb = GetClusTb(se)
            %   clusTb = GetClusTb(se.userData)
            %   clusTb = GetClusTb(se.userData.ksMeta)
            % 
            
            if iscell(src)
                % Denest or concatenate cell array inputs
                src = cat(1, src{:});
            end
            
            if isa(src, 'MSessionExplorer')
                % Take out userData struct(s) from se object
                assert(numel(src)==1, "MSessionExplorer input must be a scalar");
                src = src.userData;
            end
            
            % Collect and concatenate clusTb from all metadata struct(s)
            clusTb = cell(size(src));
            for i = 1 : numel(src)
                if isfield(src(i), 'ksMeta')
                    clusTb{i} = src(i).ksMeta.clusTb;
                elseif isfield(src, 'clusTb')
                    clusTb{i} = src(i).clusTb;
                else
                    error('Unknown input struct. Expect it to be userData or ksMeta.');
                end
            end
            clusTb = cat(1, clusTb{:});
        end
        
        function clusTb = AddRecMeta(src, clusTb)
            % Add recId, subjectId, region to clusTb
            % 
            %   clusTb = AddRecMeta(src, clusTb)
            % 
            [recId, subjectId] = NP.SE.GetID(src);
            region = NP.SE.GetRegion(src);
            clusTb.recId(:) = string(recId);
            clusTb.subjectId(:) = string(subjectId);
            clusTb.region(:) = string(region);
        end
        
        function ConvertDepth(se)
            % Make y-coordinates of all units and channels relative to surface and downward as positive
            % 
            %   NP.Unit.ConvertDepth(se)
            % 
            
            if ~isfield(se.userData, 'recMeta')
                warning("No conversion was made: spreadsheet metadata struct 'recMeta' is not found in se.userData");
                return
            end
            if isfield(se.userData, 'SurfaceCh_depth_nCh_0')
                warning("No conversion was made: found original surface info SurfaceCh_depth_nCh_0, indicating that the conversion has been done already.");
                return
            end
            S = se.userData.recMeta.SurfaceCh_depth_nCh;
            if isempty(S) || isnan(S(2))
                warning("No conversion was made: cortical surface information is missing from the SurfaceCh_depth_nCh field.");
                return
            end
            
            y0 = S(2); % surface coordinate in um is at the second position
            if ischar(y0)
                warning("No conversion was made: SurfaceCh_depth_nCh field is likely a note instead of numeric values.");
                return
            end
            if ~y0
                disp("No conversion was made: surface position is already at zero.");
                return
            end
            
            f = @(y) -(y-y0); % zero surface and flip sign such that downward is positive
            
            for i = 1 : numel(se.userData.ksMeta)
                s = se.userData.ksMeta(i);
                s.chanMapTb.ycoords = f(s.chanMapTb.ycoords);
                s.clusTb.depth = f(s.clusTb.depth);
                s.spkTb.centCoords(:,2) = f(s.spkTb.centCoords(:,2));
                s.tempTb.chanY = cellfun(f, s.tempTb.chanY, 'Uni', false);
                se.userData.ksMeta(i) = s;
            end
            
            for i = 1 : numel(se.userData.lfMeta)
                s = se.userData.lfMeta(i);
                s.chanTb.ycoords = f(s.chanTb.ycoords);
                s.chanTbKS.ycoords = f(s.chanTbKS.ycoords);
                se.userData.lfMeta(i) = s;
            end
            
            se.userData.recMeta.SurfaceCh_depth_nCh_0 = S; % backup original surface info
            se.userData.recMeta.SurfaceCh_depth_nCh(2) = 0;
        end
        
        function WriteTrial(se, sr, clusId, varargin)
            % Write waveform from trial(s) to individual audio file(s)
            %
            %   NP.Unit.WriteTrial(se, sr, clusId)
            %   NP.Unit.WriteTrial(se, sr, clusId, trialInd)
            %   NP.Unit.WriteTrial(se, sr, clusId, trialInd, tWin)
            %   NP.Unit.WriteTrial(..., 'FolderPath', '')
            %   NP.Unit.WriteTrial(..., 'FileNames', '')
            % 
            
            p = inputParser();
            p.addOptional('trialInd', 1:se.numEpochs, @isnumeric);
            p.addOptional('tWin', [-Inf Inf], @isnumeric);
            p.addParameter('FolderPath', '', @(x) isstring(x) || ischar(x));
            p.addParameter('FileNames', '', @(x) isstring(x) || ischar(x) || iscellstr(x));
            p.parse(varargin{:});
            trialInd = p.Results.trialInd;
            tWin = p.Results.tWin;
            folderPath = p.Results.FolderPath;
            fileNames = cellstr(p.Results.FileNames);
            
            if ~isempty(folderPath) && ~exist(folderPath, 'dir')
                mkdir(folderPath);
            end
            
            tRef = se.GetReferenceTime('spikeTime');
            tv = se.GetTable('taskValue');
            tv = tv(trialInd,:);
            
            for i = trialInd(:)'
                % Get file path and name
                if i > numel(fileNames) || isempty(fileNames{i})
                    fileName = strjoin([recId "trial"+tv.trialNum(i) tWinName(1)+"-"+tWinName(2)+"s.wav"], '_');
                else
                    fileName = fileNames{i};
                end
                filePath = fullfile(folderPath, fileName);
                
                % Make spike audio
                tWinGlb = tWin + tRef(i);
                sr.WriteSpikeAudio(filePath, clusId, tWinGlb);
            end
        end
        
        function clusTb2 = AlignClusTb(clusTb1, clusTb2, isRmRddCols)
            % Sort clusTb2 in the same order as clusTb1 and optioanlly remove redundant columns from clusTb2
            % 
            %   clusTb2 = AlignClusTb(clusTb1, clusTb2)
            %   clusTb2 = AlignClusTb(clusTb1, clusTb2, isRmRddCols)
            % 
            
            % Remove units from clusTb2 that do not exist in clusTb1
            m = ismember(clusTb2.clusId, clusTb1.clusId);
            clusTb2(~m,:) = [];
            
            % Sort clusTb2
            [~, I] = MMath.SortLike(clusTb2.clusId, clusTb1.clusId);
            clusTb2 = clusTb2(I,:);
            
            % Remove redundant table columns
            if nargin < 3
                isRmRddCols = false;
            end
            if isRmRddCols
                m = ismember(clusTb2.Properties.VariableNames, clusTb1.Properties.VariableNames);
                clusTb2(:,m) = [];
            end
        end
        
        % Unit cache
        function CreateUnitCache(seSrc, cacheDir)
            % Extract and save data of each unit to file
            % 
            %   CreateUnitCache(sePaths)
            %   CreateUnitCache(sePaths, cacheDir)
            %   CreateUnitCache(seArray, cacheDir)
            % 
            % Input
            %   sePaths         One or more paths of saved se objects.
            %   seArray         One or more se objects.
            %   cacheDir        The directory path to save cache in. The seArray input requires cacheDir 
            %                   to be provided. If using the sePaths input, by default, cache will be 
            %                   saved in a subfolder named 'unit_cache' in the folder where the se files live. 
            %                   This default location can be overwritten if cacheDir is provided. 
            % Cache
            %   The cache file is named as "u<clusId>.mat", and the following variables are saved or 
            %   updated (if the cache file already exists).
            %   unitInfo        A row of se.userData.ksMeta.clusTb table converted to struct.
            %   tt              A 'taskTime' table where each row is an example trial of a unique stim.
            %   tv              A 'taskValue' table where each row is an example trial of a unique stim.
            %   st              Spike times saved in a #stim-element cell array where each element is 
            %                   again a #trial-element cell array.
            % 
            % See also NP.Unit.AddData2UnitCache
            
            for i = 1 : numel(seSrc)
                % Get se object
                if iscellstr(seSrc) || isstring(seSrc)
                    % From file
                    seSrc = cellstr(seSrc);
                    se = NP.SE.LoadSession(seSrc{i});
                    seDir = fileparts(seSrc{i});
                else
                    % Use input
                    se = seSrc(i);
                    assert(~isempty(cacheDir), "When using seArray as input, cacheDir must be provided to specify the location to save the cache.");
                end
%                 fprintf("\nCache unit data from %s\n", NP.SE.GetID(se));
                
                % Determine folder to save cache
                if ~exist('cacheDir', 'var') || isempty(cacheDir)
                    cacheDir = fullfile(seDir, "unit_cache");
                end
                if ~exist(cacheDir, 'dir')
                    mkdir(cacheDir);
                end
                
                % Compute sentence PETH
                senTb = NP.TaskBaseClass.SplitBySentence(se);
                [ce, senTb] = LMV.SE.ComputeSentencePETH(senTb);
                
                % Extract and save unit data
                clusTb = NP.Unit.GetClusTb(se);
                [tt, tv] = ce.GetTable('taskTime', 'taskValue');
                stSen = arrayfun(@(x) x.GetTable('spikeTime'), senTb.se, 'Uni', false);
                stUnit = cell(height(clusTb), 1);
                for u = 1 : height(clusTb)
                    stUnit{u} = cellfun(@(x) x.(u), stSen, 'Uni', false);
                end
                
                s.unitInfo = table2struct(clusTb);
                s.tt = tt;
                s.tv = tv;
                s.st = stUnit;
                NP.Unit.AddData2UnitCache(cacheDir, clusTb.clusId, s);
            end
        end
        
        function AddData2UnitCache(cacheDir, clusId, dataStruct)
            % Add data variables to unit cache
            % 
            %   AddData2UnitCache(cacheDir, clusId, dataStruct)
            % 
            % Input
            %   cacheDir        The directory path for the cache.
            %   clusId          A vector of cluster IDs.
            %   dataStruct      A struct of data to be cached. Field names determines the variable names 
            %                   used in unit cache. Two forms of field values are supported:
            %                   1) When the number of elements of a field value equals the number of units, 
            %                      each element is saved to its corresponding unit's cache.
            %                   2) When unequal, the field value is copied as a whole and saved to every 
            %                      unit caches.
            % Cache
            %   The cache file is named as "u<clusId>.mat" with the specified variables saved or updated 
            %   (if the cache file already exists).
            % 
            % See also NP.Unit.CreateUnitCache
            
            fn = fieldnames(dataStruct);
            
            for u = 1 : numel(clusId)
                % Put data for this unit in a struct
                s = struct;
                for n = 1 : numel(fn)
                    d = dataStruct.(fn{n});
                    if numel(d) == numel(clusId)
                        % Take an element from the variable
                        if iscell(d)
                            s.(fn{n}) = d{u};
                        else
                            s.(fn{n}) = d(u);
                        end
                    else
                        % Use the variable as a whole
                        s.(fn{n}) = d;
                    end
                end
                
                % Add variables to cache
                uCachePath = fullfile(cacheDir, "u"+clusId(u)+".mat");
                if exist(uCachePath, 'file')
                    fprintf("u%i: update existing unit cache\n", clusId(u));
                    save(uCachePath, '-struct', 's', '-append');
                else
                    fprintf("u%i: create unit cache\n", clusId(u));
                    save(uCachePath, '-struct', 's');
                end
            end
        end
        
        function sUnit = LoadUnitCache(clusId, varargin)
            % Plot rasters of stim aligned to speech labels
            % 
            %   sUnit = LoadUnitCache(clusId)
            %   sUnit = LoadUnitCache(clusId, 'DataSource', 'm1')
            % 
            
            % Parse inputs
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('DataSource', 'm1', @(x) ischar(x) || isstring(x) || isempty(x));
            p.parse(varargin{:});
            dataSrc = p.Results.DataSource;
            
            % Check cluster IDs
            clusId(isnan(clusId)) = [];
            assert(~isempty(clusId), "clusId must contain at least one valid ID.\n");
            
            % Resolve the path of data source folder
            if isempty(dataSrc)
                dataSrc = "m1";
            end
            dataSrc = string(dataSrc);
            if any(strcmpi(dataSrc, ["m1", "m2"]))
                srcFolder = LMV.Data.GetAnalysisDir("data", "se_"+lower(dataSrc), "unit_cache");
            elseif strcmpi(dataSrc, "pe")
                srcFolder = LMV.Data.GetAnalysisDir("peri_event", "unit_cache");
            else
                srcFolder = dataSrc;
            end
            
            % Load cache
            srcFiles = fullfile(srcFolder, "u"+clusId+".mat");
            sUnit = arrayfun(@load, srcFiles, 'Uni', false);
        end
        
        function AddSTA2UnitCache(seSrc, cacheDir)
            % Extract and save data of each unit to file
            % 
            %   CreateUnitCache(sePaths)
            %   CreateUnitCache(sePaths, cacheDir)
            %   CreateUnitCache(seArray, cacheDir)
            % 
            % Input
            %   sePaths         One or more paths of saved se objects.
            %   seArray         One or more se objects.
            %   cacheDir        The directory path to save cache in. The seArray input requires cacheDir 
            %                   to be provided. If using the sePaths input, by default, cache will be 
            %                   saved in a subfolder named 'unit_cache' in the folder where the se files live. 
            %                   This default location can be overwritten if cacheDir is provided. 
            % Cache
            %   The cache file is named as "u<clusId>.mat", and the following variables are saved or 
            %   updated (if the cache file already exists).
            %   unitInfo        A row of se.userData.ksMeta.clusTb table converted to struct.
            %   tt              A 'taskTime' table where each row is an example trial of a unique stim.
            %   tv              A 'taskValue' table where each row is an example trial of a unique stim.
            %   st              Spike times saved in a #stim-element cell array where each element is 
            %                   again a #trial-element cell array.
            % 
            % See also NP.Unit.AddData2UnitCache
            
            for i = 1 : numel(seSrc)
                % Get se object
                if iscellstr(seSrc) || isstring(seSrc)
                    % From file
                    seSrc = cellstr(seSrc);
                    se = NP.SE.LoadSession(seSrc{i});
                    seDir = fileparts(seSrc{i});
                else
                    % Use input
                    se = seSrc(i);
                    assert(~isempty(cacheDir), "When using seArray as input, cacheDir must be provided to specify the location to save the cache.");
                end
%                 fprintf("\nCache unit data from %s\n", NP.SE.GetID(se));
                
                % Determine folder to save cache
                if ~exist('cacheDir', 'var') || isempty(cacheDir)
                    cacheDir = fullfile(seDir, "unit_cache");
                end
                if ~exist(cacheDir, 'dir')
                    mkdir(cacheDir);
                end
                
                % Compute STA
                senTb = NP.STA.SplitBySentence(se);
                [ce, senTb] = LMV.SE.ComputeSentencePETH(senTb);
                
                % Compute FTA
                
                
                % Extract and save unit data
                clusTb = NP.Unit.GetClusTb(se);
                [tt, tv] = ce.GetTable('taskTime', 'taskValue');
                stSen = arrayfun(@(x) x.GetTable('spikeTime'), senTb.se, 'Uni', false);
                stUnit = cell(height(clusTb), 1);
                for u = 1 : height(clusTb)
                    stUnit{u} = cellfun(@(x) x.(u), stSen, 'Uni', false);
                end
                
                s.unitInfo = table2struct(clusTb);
                s.tt = tt;
                s.tv = tv;
                s.st = stUnit;
                NP.Unit.AddData2UnitCache(cacheDir, clusTb.clusId, s);
            end
        end
        
    end
    
end
