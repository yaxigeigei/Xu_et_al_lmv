classdef PE
    % 
    % Main Data Structures
    % 
    %   seqData     A level-by-phase table where each element is a ssTb (see below).
    % 
    %   ssTb        An m-by-n table where each element is a seedStruct (see below).
    %               m is the number of different seed features (e.g. "T").
    %               n is the number of differently indexed sequences.
    %   
    %   seedStruct  A structrue of sequence data with a given seed.
    %               seed        The name of the seed feature, e.g. "AH".
    %               level       The name of sequence level, e.g. "phone" or "word".
    %               positions   The indices of sequence elements relative to the seed position.
    %               nElem       The numnber of elements in each sequence.
    %               seqTb       See below.
    %   
    %   seqTb       A table where rows are unique sequences and columns are various data 
    %               or labels associated with these sequences.
    %               seqStr      String representation of the sequences, e.g. "T AA".
    %               seqTge      NP.TGEvent object vectors.
    %               nSample     Number of samples (i.e. repeats) for each sequence.
    %               time0       Original seed timestamp in each sample.
    %               spikeTime, tResp, hh, ee, stats
    %                           These are added by LMV.PE.AddResp2SeqTb. In stats, the fields 
    %                           xrPkVal, xrPkLag can be further added by LMV.Linker.ComputeSeqPethXC
    %               tFeat, (featname)
    %                           These are added by LMV.PE.AddFeat2SeqTb.
    % 
    
    methods(Static)
        % Sequence extraction
        function [seedTb, senTb, seV] = ExtractPeriPhoneData(se)
            % Extract relationship between seed phones and their sentence positions
            % 
            %   [seedTb, senTb, seV] = ExtractPeriPhoneData(se)
            % 
            % Inputs
            %   se              A MSessionExplorer object.
            % 
            % Outputs
            %   seedTb          A table of seed info. Each row is a single seed event in the recording (e.g. a phone).
            %                   Columns are the follwoing:
            %                   source      The source of the event, "stim" or "prod".
            %                   name        The name of the event, e.g. "N", "OW".
            %                   time        Event onset time.
            %                   senIdx      The row index in senTb where the event is from.
            %                   wordIdx     Which word in the sentence that this event is at.
            %                   syllIdx     Which syllable in the sentence that this event is at.
            %                   phoneIdx    Which phone in the sentence that this event is at.
            % 
            %   senTb           A table of sentence info. Each row is a single sentence repeat in the recording.
            %                   Columns are the follwoing:
            %                   source      The source of the sentence, "stim" or "prod".
            %                   stimId      TIMIT style sentence ID.
            %                   text        Sentence transcript.
            %                   tge         NP.TGEvent object of the sentence.
            %                   word        Vector of NP.TGEvent objects cut to word level.
            %                   syll        Vector of NP.TGEvent objects cut to syllable level.
            %                   phone       Vector of NP.TGEvent objects cut to phone level.
            % 
            %   seV             Vectorized se with a single epoch.
            
            seV = se.Duplicate;
            
            % Add stimId to V field of sentence TGE objects
            [tt, tv] = se.GetTable("taskTime", "taskValue");
            for i = 1 : height(tt)
                stimId = tv.stimId(i);
                tt.stim(i) = tt.stim(i).SetVfield('stimId', stimId);
                tt.prod{i} = tt.prod{i}.SetVfield('stimId', repelem(stimId, numel(tt.prod{i})));
            end
            seV.SetTable('taskTime', tt);
            
            % Vectorize se
            seV.SliceSession(0, 'absolute');
            
            % Get sentence TGE objects
            tt = seV.GetTable('taskTime');
            sen = cat(1, tt.stim{:}, tt.prod{:});
            
            % Add syllables to sentence objects
            lookupPath = "lmv_syllabification_lookup.csv";
            lookup = readtable(lookupPath);
            sen = sen.AddSyllableTier(lookup);
            
            % Create sentence table
            senTb = table;
            senTb.source = sen.GetVfield('source');
            senTb.stimId = sen.GetVfield('stimId');
            senTb.text = sen.GetParentLabel;
            senTb.tge = sen;
            
            % Cache objects cut to all tiers
            senTb.word = arrayfun(@(x) Cut(x), senTb.tge, 'Uni', false);
            senTb.syll = cellfun(@(x) Cut(x), senTb.word, 'Uni', false);
            senTb.phone = cellfun(@(x) Cut(x), senTb.syll, 'Uni', false);
            
            % Create phoneme table
            seedCell = cell(0);
            for i = 1 : numel(sen)
                nSyl = 0;
                nPhn = 0;
                
                wrd = Cut(sen(i));
                for j = 1 : numel(wrd)
                    syl = Cut(wrd(j));
                    for k = 1 : numel(syl)
                        nSyl = nSyl + 1;
                        
                        phn = Cut(syl(k));
                        for m = 1 : numel(phn)
                            nPhn = nPhn + 1;
                            
                            s = struct;
                            s.source = senTb.source(i);
                            s.name = erase(phn(m).GetParentLabel, digitsPattern);
                            s.time = double(phn(m));
                            s.senIdx = i;
                            s.wordIdx = j;
                            s.syllIdx = nSyl;
                            s.phoneIdx = nPhn;
                            seedCell{end+1} = s;
                        end
                    end
                end
            end
            seedTb = struct2table([seedCell{:}]);
        end
        
        function ssTb = ComputeSeq(seedTb, senTb, varargin)
            % Find unique sequences with different starting features and at different levels
            % 
            %   ssTb = ComputeSeq(seedTb, senTb)
            %   ssTb = ComputeSeq(seedTb, senTb, ..., 'SeedFeatures', [])
            %   ssTb = ComputeSeq(seedTb, senTb, ..., 'Source', "prod")
            %   ssTb = ComputeSeq(seedTb, senTb, ..., 'Level', "phone")
            %   ssTb = ComputeSeq(seedTb, senTb, ..., 'Positions', 0)
            %   ssTb = ComputeSeq(seedTb, senTb, ..., 'IncludePartial', false)
            % 
            % Inputs
            %   seedTb, senTb       See the outputs of LMV.PE.ExtractPeriPhoneData.
            %                       senTb should have variables named by the level of interest (see 'Level' below) 
            %                       which contain NP.TGEvent objects cut to the right tier.
            %                       See the output of LMV.PE.ExtractPeriPhoneData for more details.
            %   'SeedFeatures'      Names of seed features.
            %   'Level'             "phone", "syll", or "word".
            %   'Positions'         1) A numeric vector of position indices relative to the seed feature. 
            %                          Default is 0, indicating just the seed feature itself.
            %                       2) A cell array of the vectors in 1).
            %   'Source'            "stim" or "prod" (default).
            %   'IncludePartial'    Whether or not to include sequences that partially exceed sentence boundaries. 
            %                       Default is false.
            % 
            % Output
            %   ssTb                See help LMV.PE
            % 
            % See also LMV.PE, LMV.PE.ExtractPeriPhoneData
            
            p = inputParser;
            p.addParameter('SeedFeatures', [], @(x) isstring(x) || iscellstr(x) || ischar(x) || isempty(x));
            p.addParameter('Source', "prod", @(x) any(strcmpi(x, ["stim", "prod"])));
            p.addParameter('Level', "phone", @(x) any(strcmpi(x, ["phone", "syll", "word"])));
            p.addParameter('Positions', 0, @(x) isnumeric(x) || iscell(x));
            p.addParameter('IncludePartial', false, @islogical);
            p.parse(varargin{:});
            seeds = string(p.Results.SeedFeatures);
            src = string(lower(p.Results.Source));
            lvl = string(lower(p.Results.Level));
            posInd = p.Results.Positions;
            isIncPartial = p.Results.IncludePartial;
            
            if isnumeric(posInd)
                posInd = {posInd};
            end
            
            % Get features of interest
            if isempty(seeds)
                % Find unique features and rank them by the number of occurence
                [N, C] = histcounts(categorical(seedTb.name));
                [~, I] = sort(N, 'descend');
                seeds = string(C(I));
            end
            seedTb.cat = categorical(seedTb.name, seeds, 'Ordinal', true);
            
            % Remove irrelevant instances
            %   Those associated with sentences of other sources
            %   Those outside the features of interest
            m = ismember(seedTb.senIdx, find(senTb.source==src)) & ~ismissing(seedTb.cat);
            seedTb = seedTb(m,:);
            
            % 
            nFeat = numel(seeds);
            nSeq = numel(posInd);
            outCell = cell(nFeat, nSeq);
            for i = 1 : nSeq
                % Extract sequence context
                seedTb = LMV.PE.ExtractSeq(seedTb, senTb, lvl, posInd{i}, isIncPartial);
                hasSeq = ~ismissing(seedTb.seqStr);
                
                for j = 1 : nFeat
                    % Find unique sequences for this feature
                    isFeat = seedTb.name == seeds(j);
                    seqList = unique(seedTb.seqStr(hasSeq & isFeat));
                    
                    % Initialize seedStruct
                    s = struct;
                    s.seed = seeds(j);
                    s.level = lvl;
                    s.positions = posInd{i};
                    s.nElem = numel(posInd{i});
                    
                    tb = table;
                    tb.seqStr = seqList;
                    for k = 1 : numel(seqList)
                        % Find sequences
                        isSeq = isFeat & seedTb.seqStr==seqList(k);
                        I = find(isSeq, 1);
                        tge = seedTb.seqTge{I};
                        tb.seqTge{k} = tge;
                        tb.nSample(k) = sum(isSeq);
                        tb.time0{k} = seedTb.time(isSeq);
                    end
                    s.seqTb = tb;
                    
                    outCell{j,i} = s;
                end
            end
            
            ssTb = cell2table(num2cell(outCell), 'RowNames', seeds, 'VariableNames', "seq"+(1:nSeq));
        end
        
        function seedTb = ExtractSeq(seedTb, senTb, lvl, posInd, isIncludePartial)
            % Extract sequences around events
            % 
            %   seedTb = ExtractSeq(seedTb, senTb, lvl, posInd, isIncludePartial)
            % 
            
            for i = 1 : height(seedTb)
                % Get the full object sequence from the sentence table
                k = seedTb.senIdx(i);
                tge = senTb.(lvl){k};
                
                % Get the indices of peri-event elements
                k = seedTb.(lvl+"Idx")(i);
                ind = k + posInd;
                
                % Handle out-of-range indices
                isOut = ind < 1 | ind > numel(tge);
                if (~isIncludePartial && any(isOut)) || all(isOut)
                    seedTb.seqStr(i) = string(missing);
                    seedTb.seqTge{i} = tge([]);
                    continue
                end
                ind(isOut) = [];
                
                % Extract the peri-event sequence
                seedTb.seqStr(i) = erase(join(tge(ind).GetParentLabel, " "), digitsPattern);
                seedTb.seqTge{i} = tge(ind) - seedTb.time(i);
            end
        end
        
        % Adding data to seqTb
        function seqTb = AddResp2SeqTb(seqTb, se, tWin, uInd, mdls)
            % Add peri-event spike times and PETHs to seqTb
            % 
            %   seqTb = AddResponses2SeqTb(seqTb, se, tWin, uInd, mdls)
            % 
            % Inputs
            %   seqTb       See help LMV.PE
            %   seV         Vectorized se with a single epoch.
            %   tWin        A peri-event time window to compute PETHs.
            % 
            % Output
            %   seqTb       Same as input but with the following columns added:
            %               spikeTime
            %               t
            %               hh
            %               ee
            %               stats
            % 
            % See also LMV.PE
            
            clusTb = NP.Unit.GetClusTb(se);
            if exist('uInd', 'var') && ~isempty(uInd)
                clusTb = clusTb(uInd,:);
            else
                uInd = [];
            end
            
            for k = 1 : height(seqTb)
                % Resampling parameters
                tBin = 0.01;
                
                % Slice out perievent spike times and rates
                t0 = seqTb.time0{k};
                st = se.SliceEventTimes('spikeTime', {t0 + tWin}, [], uInd);
                % sr = se.SliceTimeSeries('spikeRate', {t0 + tWin + [-tBin tBin]}, [], uInd+1);
                
                % Re-reference timestamps to seed event times
                for i = 1 : width(st)
                    st.(i) = cellfun(@(x,y) x-y, st.(i), num2cell(t0), 'Uni', false);
                end
                % sr.time = cellfun(@(x,y) x-y, sr.time, num2cell(t0), 'Uni', false);
                
                % Resample perievent spike rates
                tEdges = (tWin(1) : tBin : tWin(2))';
                tCenters = MMath.BinEdges2Centers(tEdges);
                sr = se.ResampleEventTimes(st, tEdges, 'Normalization', 'countdensity');
                rSeq = cellfun(@(x) MNeuro.Filter1(x, 1/tBin, 'Gaussian', 0.015), sr{:,2:end}, 'Uni', false);
                
                % Compute PETHs
                [hh, ee, stats] = MNeuro.MeanTimeSeries(rSeq);
                stats.clusId = clusTb.clusId;
                
                % Find response at optimal encoding time
                seqTb.optiResp0{k} = NaN(size(st));
                stats.optiTime(:) = 0.15; % default
                stats.optiResp(:) = NaN;
                if exist('mdls', 'var')
                    for u = 1 : height(stats)
                        if iscell(mdls)
                            mdl = mdls{u};
                        else
                            mdl = mdls;
                        end
                        if ~isempty(mdl)
                            stats.optiTime(u) = mdl.r2t;
                        end
                        [~, iOpti] = min(abs(tCenters - stats.optiTime(u)));
                        seqTb.optiResp0{k}(:,u) = cellfun(@(x) x(iOpti), rSeq(:,u));
                        stats.optiResp(u) = hh(iOpti, u);
                    end
                end
                
                seqTb.spikeTime{k} = st{:,:};
                seqTb.tResp{k} = tCenters;
                seqTb.hh{k} = hh;
                seqTb.ee{k} = ee;
                seqTb.stats{k} = stats;
            end
        end
        
        function seqTb = AddFeat2SeqTb(seqTb, se, tWin, tbName, colNames, featName)
            % Add peri-event feature arrays to seqTb
            % 
            %   seqTb = AddFeat2SeqTb(seqTb, se, tWin, tbName, colNames, featName)
            % 
            % Inputs
            %   seqTb       See help LMV.PE
            %   seV         Vectorized se with a single epoch.
            %   tWin        A peri-event time window to compute PETHs.
            % 
            % Output
            %   seqTb       Same as input but with the following columns added:
            %               tFeat       A vector of feature timestamps for each sequence.
            %               (featName)  A s-by-1 cell array where s is the number of samples for each sequence.
            %                           Each element in the cell array is a t-by-d numeric array. t is the number 
            %                           of time points, d is the number of sub-features, e.g. bins of spectrogram.
            % 
            
            if ~exist('featName', 'var')
                featName = tbName;
            end
            if ~exist('colNames', 'var')
                colNames = [];
            end
            
            for k = 1 : height(seqTb)
                % Slice out perievent feature timeseries
                t0 = seqTb.time0{k};
                ts = se.SliceTimeSeries(tbName, {t0 + tWin}, [], colNames);
                
                % Re-reference spike times to event times
                ts.time = cellfun(@(x,y) x-y, ts.time, num2cell(t0), 'Uni', false);
                
                % Resample timeseries to have consistent samples
                fs = 100;
                tEdges = (tWin(1) : 1/fs : tWin(2))';
                ts = se.ResampleTimeSeries(ts, tEdges);
                
                % Use a NP.CodingExplorer to organize data arrays
                ce = NP.CodingExplorer;
                ce.SetTable(tbName, ts, 'timeSeries');
                [t, arr] = ce.GetArray(tbName, 'DimCat', 0);
                
                seqTb.tFeat{k} = t;
                seqTb.(featName){k} = arr;
            end
        end
        
        function seqTb = AddWaveform2SeqTb(seqTb, se, tWin, colNames)
            % Add peri-event waveforms to seqTb
            % 
            %   seqTb = AddWaveform2SeqTb(seqTb, se, tWin, tbName, colNames)
            % 
            % Inputs
            %   seqTb       See help LMV.PE
            %   seV         Vectorized se with a single epoch.
            %   tWin        A peri-event time window to compute PETHs.
            % 
            % Output
            %   seqTb       Same as input but with the following columns added:
            %               tFeat       A vector of feature timestamps for each sequence.
            %               (featName)  A s-by-1 cell array where s is the number of samples for each sequence.
            %                           Each element in the cell array is a t-by-d numeric array. t is the number 
            %                           of time points, d is the number of sub-features, e.g. bins of spectrogram.
            % 
            
            tbName = "ni";
            if ~exist('colNames', 'var') || isempty(colNames)
                colNames = ["mic", "speaker1"];
            end
            colNames = string(colNames);
            
            for k = 1 : height(seqTb)
                % Slice out perievent feature timeseries
                t0 = seqTb.time0{k};
                ts = se.SliceTimeSeries(tbName, {t0 + tWin}, [], colNames);
                
                % Re-reference timestamps to event times
                ts.time = cellfun(@(x,y) x-y, ts.time, num2cell(t0), 'Uni', false);
                
                % Use a NP.CodingExplorer to organize data arrays
                ce = NP.CodingExplorer;
                ce.SetTable(tbName, ts, 'timeSeries');
                [T, W] = ce.GetArray(tbName, 'DimCat', 0);
                nSp = min(cellfun(@numel, T));
                
                seqTb.tWaveform{k} = T{1}(1:nSp);
                for i = 1 : numel(colNames)
                    arr = cellfun(@(x) x(1:nSp,i), W, 'Uni', false);
                    seqTb.(colNames(i)){k} = cat(2, arr{:});
                end
            end
        end
        
        % 
        function clusTbs = ExtractSeqResp0(s, targets)
            % 
            %   clusTbs = ExtractSeqResp(s, targets)
            % 
            
            clusTbs = cell(size(targets));
            
            % Unpack input
            se = s.se;
            recId = NP.SE.GetID(se);
            clusTb = NP.Unit.GetClusTb(se);
            
            % Load phone stRFs
            mdlTb = LMV.TRF.LoadModels("phone_"+targets, recId);
            [~, trfInd] = MMath.SortLike(mdlTb.clusId, clusTb.clusId);
            mdlTb = mdlTb(trfInd,:);
            
            % Loop through targets
            for tIdx = 1 : numel(targets)
                % Keep segmental models
                target = targets(tIdx);
                mdls = mdlTb.("phone_"+target);
                mdls = LMV.RF.UpdateModelR2(mdls);
                isSeg = ~cellfun(@(x) isempty(x) || x.r2 < 0.01 || x.null.r2Pval > 0.05, mdls);
                if ~any(isSeg)
                    fprintf("%s does not have any segmental phone_%s models\n", recId, target);
                    continue
                end
                segTb = clusTb(isSeg,:);
                segTb.mdls = mdls(isSeg);
                
                % Extract responses and features from se
                switch target
                    case 'stim'
                        tWin = [-0.1 0.4];
                    case 'prod'
                        tWin = [-0.4 0.1];
                    case 'feedback'
                        tWin = [-0.1 0.4];
                    otherwise
                        error("'%s' is not a supported target name.", target);
                end
                
                S = s.seqData.(target){"phone"}{:,:};
                
                for k = 1 : numel(S)
                    seqTb = S{k}.seqTb;
                    if isempty(seqTb)
                        continue
                    end
                    seqTb(seqTb.nSample<3,:) = [];
                    seqTb = LMV.PE.AddResp2SeqTb(seqTb, se, tWin, find(isSeg), mdls(isSeg));
                    S{k}.seqTb = seqTb;
                end
                
                % Mask out sequences with no further branching
                nSeq = cellfun(@(x) height(x.seqTb), S);
                for k = size(nSeq,2) : -1 : 2
                    m = nSeq(:,k) <= nSeq(:,k-1);
                    nSeq(m,k) = NaN;
                end
                
                % Loop through units
                for uIdx = 1 : height(segTb)
                    % Get the optimal encoding time
                    mdl = segTb.mdls{uIdx};
                    tResp = S{1}.seqTb.tResp{1};
                    [~, iOpt] = min(abs(tResp - mdl.r2t));
                    
                    % Find the spike rate at optimal encoding time
                    R = cell(size(S));
                    for i = 1 : numel(S)
                        seqTb = S{i}.seqTb;
                        if isempty(seqTb) || isnan(nSeq(i))
                            R{i} = NaN;
                            continue
                        end
                        % hh = cellfun(@(x) x(:,isSeg), seqTb.hh, 'Uni', false);
                        % R{i} = cellfun(@(y) y(iOpt,uIdx), hh);
                        R{i} = cellfun(@(x) x(:,uIdx), seqTb.optiResp0, 'Uni', false);
                    end
                    rTb = s.seqData.(target){"phone"}; % just to borrow the headers
                    rTb{:,:} = R;
                    segTb.seqR{uIdx} = rTb;
                end
                
                clusTbs{tIdx} = segTb;
            end
            
        end
        
        function clusTbs = ExtractSeqResp(s, targets)
            % 
            %   clusTbs = ExtractSeqResp(s, targets)
            % 
            
            clusTbs = cell(size(targets));
            
            % Unpack input
            se = s.se;
            recId = NP.SE.GetID(se);
            clusTb = NP.Unit.GetClusTb(se);
            
            % Load phone stRFs
            mdlTb = LMV.TRF.LoadModels("phone_"+targets, recId);
            [~, trfInd] = MMath.SortLike(mdlTb.clusId, clusTb.clusId);
            mdlTb = mdlTb(trfInd,:);
            
            % Loop through targets
            for tIdx = 1 : numel(targets)
                % Keep segmental models
                target = targets(tIdx);
                mdls = mdlTb.("phone_"+target);
                mdls = LMV.RF.UpdateModelR2(mdls);
                isSeg = ~cellfun(@(x) isempty(x) || x.r2 < 0.01 || x.null.r2Pval > 0.05, mdls);
                if ~any(isSeg)
                    fprintf("%s does not have any segmental phone_%s models\n", recId, target);
                    continue
                end
                segTb = clusTb(isSeg,:);
                segTb.mdls = mdls(isSeg);
                
                % Extract responses and features from se
                switch target
                    case 'stim'
                        tWin = [-0.1 0.4];
                    case 'prod'
                        tWin = [-0.4 0.1];
                    case 'feedback'
                        tWin = [-0.1 0.4];
                    otherwise
                        error("'%s' is not a supported target name.", target);
                end
                
                S = s.seqData.(target){"phone"}{:,:};
                
                for k = 1 : numel(S)
                    seqTb = S{k}.seqTb;
                    if isempty(seqTb)
                        continue
                    end
                    seqTb(seqTb.nSample<3,:) = [];
                    seqTb = LMV.PE.AddResp2SeqTb(seqTb, se, tWin, find(isSeg), mdls(isSeg));
                    S{k}.seqTb = seqTb;
                end
                
                % Mask out sequences with no further branching
                nSeq = cellfun(@(x) height(x.seqTb), S);
                for k = size(nSeq,2) : -1 : 2
                    m = nSeq(:,k) <= nSeq(:,k-1);
                    nSeq(m,k) = NaN;
                end
                
                % Extract signle sequence responses for each unit
                for uIdx = 1 : height(segTb)
                    T = cell(size(S));
                    for i = 1 : size(S,1)
                        maxLen = find(~isnan(nSeq(i,:)), 1, 'last');
                        for j = maxLen
                            seqTb = S{i,j}.seqTb;
                            if isempty(seqTb) || isnan(nSeq(i,j))
                                T{i} = [];
                                continue
                            end
                            rTb = table;
                            rTb.clusId = repelem(segTb.clusId(uIdx), sum(seqTb.nSample))';
                            rTb.seqStr = repelem(seqTb.seqStr', seqTb.nSample)';
                            for k = 1 : size(S,2)
                                if k <= maxLen
                                    seqStrK = cellfun(@(x) erase(x(1:min(end,k)).GetAllParentLabel, digitsPattern), seqTb.seqTge);
                                else
                                    seqStrK = string(NaN(size(seqTb.seqTge)));
                                end
                                rTb.("seqStr"+k) = repelem(seqStrK', seqTb.nSample)';
                            end
                            rTb.time0 = cat(1, seqTb.time0{:});
                            rTb.optiResp = cell2mat(cellfun(@(x) x(:,uIdx), seqTb.optiResp0, 'Uni', false));
                            T{i,j} = rTb;
                        end
                    end
                    segTb.seqR{uIdx} = cat(1, T{:});
                end
                
                clusTbs{tIdx} = segTb;
            end
            
        end
        
        function m = ComputeSeqMod(rTb, isShuffle)
            % 
            %   m = ComputeSeqMod(rTb)
            %   m = ComputeSeqMod(rTb, isShuffle)
            % 
            
            if ~exist('isShuffle', 'var')
                isShuffle = false;
            end
            
            m = [0 NaN NaN];
            maxLen = numel(split(rTb.seqStr(1), " "));
            for i = 1 : maxLen-1
                parentSeqList = unique(rTb.("seqStr"+i), 'stable');
                d = NaN(size(parentSeqList));
                for j = 1 : numel(parentSeqList)
                    isSeq = rTb.("seqStr"+i)==parentSeqList(j);
                    subTb = rTb(isSeq,:);
                    if isShuffle
                        subTb.optiResp = randsample(subTb.optiResp, height(subTb));
                    end
                    chTb = groupsummary(subTb, "seqStr"+(i+1), "mean", "optiResp");
                    d(j) = (max(chTb.mean_optiResp) - min(chTb.mean_optiResp)) / mean(subTb.optiResp);
                end
                m(i+1) = mean(d, 'omitmissing');
            end
        end
        
        function [fTb, rTb] = MakeFlattenedTables(C, stRFs)
            % 
            %
            %   [fTb, rTb] = MakeFlattenedTables(C, stRFs)
            % 
            
            % Construct feature table
            nSeed = numel(C);
            tbs = cell(nSeed, 1);
            for i = 1 : numel(C)
                s = C{i};
                nSeq = height(s.seqTb);
                if ~nSeq
                    continue
                end
                tb = table;
                tb.seed = repmat(s.seed, [nSeq 1]);
                tb.level = repmat(s.level, [nSeq 1]);
                tb.positions = repmat({s.positions}, [nSeq 1]);
                tb.nElem = repmat(s.nElem, [nSeq 1]);
                tb.seqStr = s.seqTb.seqStr;
                tb.nRep = s.seqTb.nSample;
                if ismember("mel", s.seqTb.Properties.VariableNames)
                    tb.tFeat = cellfun(@(x) x{1}, s.seqTb.tFeat, 'Uni', false);
                    tb.mel = cellfun(@(x) mean(cat(3,x{:}), 3, 'omitmissing'), s.seqTb.mel, 'Uni', false);
                end
                if ismember("tWaveform", s.seqTb.Properties.VariableNames)
                    tb.tWaveform = s.seqTb.tWaveform;
                    tb.speaker1 = s.seqTb.speaker1;
                    tb.mic = s.seqTb.mic;
                end
                
                tbs{i} = tb;
            end
            fTb = cat(1, tbs{:});
            
            
            % Construct response table
            nUnit = height(s.seqTb.stats{1});
            tbs = cell(nUnit, nSeed);
            for i = 1 : numel(C)
                s = C{i};
                nSeq = height(s.seqTb);
                if ~nSeq
                    continue
                end
                for j = 1 : nUnit
                    tb = table;
                    tb.seqStr = s.seqTb.seqStr;
                    tb.clusId = cellfun(@(x) x.clusId(j), s.seqTb.stats);
                    tb.tResp = s.seqTb.tResp;
                    tb.resp = cellfun(@(x) x(:,j), s.seqTb.hh, 'Uni', false);
                    tbs{j,i} = tb;
                end
            end
            rTb = cat(1, tbs{:});
            
            
            % Find the indices of optimal coding time
            tResp = rTb.tResp{1};
            iDefaut = round(numel(tResp)/2);
            rTb.opResp(:) = NaN;
            for uIdx = 1 : nUnit
                mdl = stRFs{uIdx};
                if isempty(mdl)
                    iOpt = iDefaut;
                else
                    [mdl.r2, mdl.r2t, mdl.r2Idx] = LMV.RF.FindPeakR2(mdl.r2each, mdl.dt);
                    [~, iOpt] = min(abs(tResp - mdl.r2t));
                    
                    % Use a default +-0.15s if too little variance explained or optimal model is not significant
                    if mdl.r2 < 0.01 || mdl.null.r2Pval > 0.05
                        iOpt = iDefaut;
                    end
                end
                isUnit = rTb.clusId == tbs{uIdx,1}.clusId(1);
                rTb.iOpti(isUnit) = iOpt;
                rTb.opResp(isUnit) = cellfun(@(x) x(iOpt), rTb.resp(isUnit));
            end
            
        end
        
        % Plotting
        function PlotPeriEventSeqArray(C, varargin)
            % Plot an array of panels each contains rasters and PETHs separated by different speech sequences
            % 
            %   PlotPeriEventSeqArray(C)
            % 
            
            p = inputParser;
            p.addParameter('UnitIdx', 1, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('MaxSpikeRate', [], @(x) isscalar(x) && isnumeric(x));
            p.addParameter('FeatureName', [], @(x) ischar(x) || isstring(x) || iscellstr(x) || isempty(x));
            p.parse(varargin{:});
            uIdx = p.Results.UnitIdx;
            maxR = p.Results.MaxSpikeRate;
            feat = p.Results.FeatureName;
            
            if isempty(maxR) && isempty(feat)
                % Find the maximum PETH peak spike rate
                maxR = zeros(size(C));
                for i = 1 : numel(C)
                    if isempty(C{i}.seqTb)
                        continue
                    end
                    maxR(i) = max(cellfun(@(x) x.pkVal(uIdx), C{i}.seqTb.stats));
                end
                maxR = max(maxR(:));
            end
            
            % Plot
            [nRow, nCol] = size(C);
            tl = tiledlayout(nRow, nCol);
            tl.Padding = "compact";
            for i = 1 : nRow
                for j = 1 : nCol
                    ax = nexttile;
                    if isempty(feat)
                        LMV.PE.PlotSeqResp(ax, C{i,j}, 'MaxSpikeRate', maxR, 'UnitIdx', uIdx);
                    elseif strcmpi(feat, "tgtiers")
                        LMV.PE.PlotSeqTGTiers(ax, C{i,j});
                    else
                        LMV.PE.PlotSeqFeat(ax, C{i,j}, feat, 'UnitIdx', uIdx);
                    end
                end
            end
        end
        
        function PlotSeqResp(ax, s, varargin)
            % Plot rasters and PETHs separated by different speech sequences
            % 
            %   PlotSeqResp(ax, s)
            % 
            
            p = inputParser;
            p.addParameter('UnitIdx', 1, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('MaxSpikeRate', [], @(x) isscalar(x) && isnumeric(x));
            p.parse(varargin{:});
            uIdx = p.Results.UnitIdx;
            maxR = p.Results.MaxSpikeRate;
            
            % Organize table
            tb = s.seqTb;
            if isempty(tb)
                return
            end
            tb(tb.nSample<3,:) = []; % remove seqs that have too few samples
            if isempty(tb)
                return
            end
            r = round(cellfun(@(x) x.optiResp(uIdx), tb.stats));
            [~, I] = sort(r, 'descend');
            tb = tb(I,:);
            
            % PETH
            hh = cat(3, tb.hh{:});
            hh = squeeze(hh(:,uIdx,:));
            ee = cat(3, tb.ee{:});
            ee = squeeze(ee(:,uIdx,:));
            if isempty(maxR)
                maxR = max(hh(:));
            end
            hh = hh ./ maxR;
            ee = ee ./ maxR;
            MPlot.PlotHistStack(tb.tResp{1}, hh, ee, 'Scaling', 0.6);
            
            % Raster
            st = cellfun(@(x) x(:,uIdx), tb.spikeTime, 'Uni', false);
            MPlot.PlotRasterStack(st, 1:height(tb), 'HeightScale', 0.6, 'Parent', ax);
            
            % Labels
            for i = 1 : height(tb)
                tb.seqTge{i}.Plot(i-0.5+[0 .2], 'FontSize', 8);
                hold(ax, 'on');
            end
            
            r = round(cellfun(@(x) x.optiResp(uIdx), tb.stats));
            labelArray = r'+"Hz";
            labelArray(ismissing(labelArray)) = "";
            tickLabels = sprintf('%s\n', labelArray(:));
            ax.YTickLabel = strtrim(tickLabels);
            
            ax.Title.String = sprintf("/%s/ %+i %s", MLing.ARPA2IPA(s.seed), s.nElem, s.level);
            MPlot.Axes(ax);
        end
        
        function PlotSeqFeat(ax, s, feat, varargin)
            % Plot feature timeseries separated by different speech sequences
            % 
            %   PlotSeqFeat(ax, s)
            % 
            
            p = inputParser;
            p.addParameter('UnitIdx', 1, @(x) isscalar(x) && isnumeric(x));
            p.parse(varargin{:});
            uIdx = p.Results.UnitIdx;
            
            % Organize table
            tb = s.seqTb;
            if isempty(tb)
                return
            end
            tb(tb.nSample<3,:) = []; % remove seqs that have too few samples
            if isempty(tb)
                return
            end
            pkSpkRate = round(cellfun(@(x) x.pkVal(uIdx), tb.stats));
            [~, I] = sort(pkSpkRate, 'descend');
            tb = tb(I,:);
            
            % Compute mean feature time series
            for i = height(tb) : -1 : 1
                MM{i} = mean(cat(3, tb.(feat){i}{:}), 3, 'omitnan')';
            end
            MPlot.PlotHeatmapStack(tb.tFeat{i}{1}, MM, 'Scaling', 0.6, 'Parent', ax);
            
            % Labels
            for i = 1 : height(tb)
                tb.seqTge{i}.Plot(i-0.5+[0 .2], 'FontSize', 8);
                hold(ax, 'on');
            end
            
            pkSpkRate = round(cellfun(@(x) x.pkVal(uIdx), tb.stats));
            labelArray = pkSpkRate'+"Hz";
            labelArray(ismissing(labelArray)) = "";
            tickLabels = sprintf('%s\n', labelArray(:));
            ax.YTickLabel = strtrim(tickLabels);
            
%             ax.Title.String = sprintf("/%s/ %+i %s, %i seqs", MLing.ARPA2IPA(s.feat), s.nElem, s.level, height(tb));
            ax.Title.String = sprintf("/%s/ %+i %s", MLing.ARPA2IPA(s.feat), s.nElem, s.level);
            MPlot.Axes(ax);
        end
        
        function PlotSeqTGTiers(ax, s, varargin)
            % Plot TextGrid tiers for different speech sequences
            % 
            %   PlotSeqTGTiers(ax, s)
            % 
            
            p = inputParser;
            p.addParameter('UnitIdx', 1, @(x) isscalar(x) && isnumeric(x));
            p.parse(varargin{:});
            uIdx = p.Results.UnitIdx;
            
            % Organize table
            tb = s.seqTb;
            if isempty(tb)
                return
            end
            tb(tb.nSample<3,:) = []; % remove seqs that have too few samples
            if isempty(tb)
                return
            end
            
            % Labels
            for i = 1 : height(tb)
                tb.seqTge{i}.PlotTiers(i-0.45, 0.9/3);
                hold(ax, 'on');
            end
            
            labelArray = "n=" + tb.nSample';
            labelArray(ismissing(labelArray)) = "";
            tickLabels = sprintf('%s\n', labelArray(:));
            ax.YTick = 1 : height(tb);
            ax.YTickLabel = strtrim(tickLabels);
            ax.YLim = [0.5 height(tb)+0.5];
            
            ax.Title.String = sprintf("/%s/ %+i %s", MLing.ARPA2IPA(s.feat), s.nElem, s.level);
            MPlot.Axes(ax);
        end
        
        function PlotArticResp(mdl, seqTb, featName, varargin)
            % Plot spike raster with top artic trajectories
            % 
            %   PlotArticResp(mdl, seqTb, featName)
            %   PlotArticResp(mdl, seqTb, featName, ..., 'Prctile', [75 99])
            %   PlotArticResp(mdl, seqTb, featName, ..., 'Interval', 1)
            %   PlotArticResp(mdl, seqTb, featName, ..., 'Parent', nexttile)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('Prctile', 75, @(x) isnumeric(x) && numel(x)<3 && ~isempty(x));
            p.addParameter('Interval', 1, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Parent', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
            prc = p.Results.Prctile;
            itvl = p.Results.Interval;
            ax = p.Results.Parent;
            
            % Get feature timeseries
            tFeat = seqTb.tFeat{1}{1};
            A = cat(1, seqTb.(featName){:});
            A = cat(2, A{:})';
            
            % Slice out spike times and re-reference to the optimal coding time
            st = cat(1, seqTb.spikeTime{:});
            st = cellfun(@(x) x - mdl.r2t, st, 'Uni', false);
            
            % Sort samples by values at time zero
            [~, I] = min(abs(tFeat - 0));
            [a, I] = sort(A(:,I), 1, "ascend");
            A = A(I,:);
            st = st(I);
            
            % Remove outliers
            % m = ~isoutlier(a, "median", "ThresholdFactor", 2);
            m = ~isoutlier(a);
            a = a(m);
            A = A(m,:);
            st = st(m);
            
            % Find samples in percentile range
            prcVals = prctile(a, prc);
            if isscalar(prcVals)
                m = a > prcVals;
            else
                m = a > prcVals(1) & a < prcVals(2);
            end
            a = a(m);
            A = A(m,:);
            st = st(m);
            
            % Subsampling
            iSub = 1 : itvl : numel(st);
            N = numel(iSub);
            
            % Compute mean feature timeseries
            [Am, Asd, Ase] = MMath.MeanStats(A, 1);
            k = max(Am+Ase);
            Am = Am / k * N;
            Asd = Asd / k * N;
            Ase = Ase / k * N;
            
            % Compute mean spike rate
            tWin = seqTb.tResp{1}([1 end]) - mdl.r2t;
            % tWin = [-0.1 0.1];
            binEdges = tWin(1) : 0.01 : tWin(2);
            tResp = MMath.BinEdges2Centers(binEdges);
            R = cellfun(@(x) MNeuro.Filter1(x, [], 'bin', binEdges), st, 'Uni', false);
            R = cellfun(@(x) MNeuro.Filter1(x, 100, 'Gaussian', 0.015), R, 'Uni', false);
            [Rm, Rse] = MNeuro.MeanTimeSeries(R);
            k = max(Rm+Rse);
            Rm = Rm / k * N;
            Rse = Rse / k * N;
            
            % Plotting
            if isempty(ax)
                ax = nexttile;
            end
            MPlot.PlotRaster(st(iSub), [], N/100, 'ColorArray', [0 0 0], 'Parent', ax);
            MPlot.ErrorShade(tFeat, Am, Ase, 'Color', 'b', 'Parent', ax);
            plot(ax, tFeat, Am, 'Color', 'b');
            MPlot.ErrorShade(tResp, Rm, Rse, 'Color', 'k', 'Parent', ax);
            plot(ax, tResp, Rm, 'Color', 'k');
            
            ax.XLim = [-1 1]*0.1; % + mdl.r2t;
            ax.YLim = [0 N+1];
            ax.YTick = [];
            ax.XLabel.String = "Time (s)"; %  + \Deltat
            ax.YLabel.String = sprintf("%.1f Hz", k);
            ax.Title.String = sprintf("%s (%s pct, n=%i)", NP.Artic.GetLabels(featName), strjoin(string(prc), '-'), N);
            MPlot.Axes(ax);
        end
        
    end
    
end
