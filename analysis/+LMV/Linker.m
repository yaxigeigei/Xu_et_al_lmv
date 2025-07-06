classdef Linker
    % 
    % Main Data Structures
    % 
    %   triTb       An m-by-n table of seedStructs. m is the number of different seeds (e.g. "T"). 
    %               n is the number of different targets including "stim", "feedback", and "prod".
    %   
    %   
    % See also LMV.PE
    
    properties(Constant)
        types = ["mirror", "bridge", "feedback"];
        bridgeFN = [410100450]; % due to false negative responsiveness; 410100441
        currentModel = "smooth_lm";
        scoreTh = 0.1;
        respTimeWin = [-0.4 0.4];
    end
    
    methods(Static)
        % Utilities
        function [cidUni, cid] = GetSelectedClusId(type, recId)
            % Return handpicked cluster IDs for a given recording
            % 
            %   [cidUni, cid] = GetSelectedClusId(type, recId)
            % 
            
            type = string(type);
            recId = string(recId);
            
            if numel(recId) > 1
                [cidUni, cid] = arrayfun(@(x) LMV.Linker.GetSelectedClusId(type, x), recId, 'Uni', false);
                cidUni = [cidUni{:}];
                cid = [cid{:}];
                return
            end
            
            if type == "mirror"
                switch recId
                    % IFG
                    % STG
                    case 'NP35_B2', cid = []; % v2
                    case 'NP43_B1', cid = []; % v2
                    case 'NP50_B3', cid = []; % v2
                    case 'NP53_B1', cid = []; % v2
                    % vPrCG
                    case 'NP38_B5', cid = []; % none; v2
                    case 'NP38_B6', cid = []; % none; v2
                    case 'NP52_B1', cid = [94]; % v2
                    case 'NP52_B2', cid = [153]; % v2
                    case 'NP54_B1', cid = []; % v2
                    % mPrCG
                    case 'NP41_B1', cid = [170 224 400 416 424]; % v1.5
                    case 'NP44_B2', cid = [462 479]; % v2
                    case 'NP44_B3', cid = []; % v2
                    case 'NP45_B1', cid = [4 165 286 314]; % v2
                    case 'NP45_B2', cid = []; % v2
                    case 'NP46_B1', cid = [0 3 5 13 30 98 236 248]; % v2
                    case 'NP47_B3', cid = []; % v1.5
                    case 'NP56_B1', cid = [10714]; % v2
                    otherwise, cid = [];
                end
            elseif type == "bridge"
                switch recId
                    % IFG
                    % vPrCG
                    case 'NP38_B5', cid = []; % none; v2
                    case 'NP38_B6', cid = []; % none; v2
                    case 'NP52_B1', cid = [16]; % v2
                    case 'NP52_B2', cid = []; % v2
                    case 'NP54_B1', cid = []; % v2
                    % mPrCG
                    case 'NP41_B1', cid = [155 211 450 463 454]; % 441; v1.5
                    case 'NP44_B2', cid = []; % v2
                    case 'NP44_B3', cid = [575 594]; % v2
                    case 'NP45_B1', cid = []; % v2
                    case 'NP45_B2', cid = []; % v2
                    case 'NP46_B1', cid = [250 271]; % v2
                    case 'NP47_B3', cid = []; % v1.5
                    case 'NP56_B1', cid = []; % v2
                    otherwise, cid = [];
                end
            elseif type == "feedback"
                switch recId
                    % IFG
                    % STG
                    case 'NP35_B2', cid = [795 829]; % v2
                    case 'NP43_B1', cid = [161 210 307 319 381 413 428 473 544 549 555]; % v2
                    case 'NP50_B3', cid = [452 456]; % v2
                    case 'NP53_B1', cid = []; % v2
                    % vPrCG
                    case 'NP38_B5', cid = [266]; % none; v2
                    case 'NP38_B6', cid = [105 177]; % none; v2
                    case 'NP52_B1', cid = []; % v2
                    case 'NP52_B2', cid = []; % v2
                    case 'NP54_B1', cid = []; % v2
                    % mPrCG
                    case 'NP41_B1', cid = []; % v1.5
                    case 'NP44_B2', cid = []; % v2
                    case 'NP44_B3', cid = [585 587 595]; % v2
                    case 'NP45_B1', cid = []; % v2
                    case 'NP45_B2', cid = [52]; % v2
                    case 'NP46_B1', cid = [128]; % v2
                    case 'NP47_B3', cid = []; % v1.5
                    case 'NP56_B1', cid = []; % v2
                    otherwise, cid = [];
                end
            else
                cid = [];
            end
            
            cidUni = NP.Unit.GetBaseClusId(recId) + cid;
        end
        
        function uu = GetExampleUnitInfo(type)
            % Get info of example units
            % 
            %   uu = GetExampleUnitInfo(type)
            % 
            uu = cell(0);
            switch char(type)
                case 'mirror'
                    uu{end+1} = struct('clusId', 460100003, 'seqStr', "earth's the matter", 'seed', "DH");
                    uu{end+1} = struct('clusId', 520200153, 'seqStr', "will you tell", 'seed', "Y");
                    uu{end+1} = struct('clusId', 460100000, 'seqStr', "got enough blankets", 'seed', "AH");
                    uu{end+1} = struct('clusId', 460100248, 'seqStr', "get a tax", 'seed', "AH");
                    uu{end+1} = struct('clusId', 460100030, 'seqStr', "plenty of time", 'seed', "AH"); % two activations
                    uu{end+1} = struct('clusId', 460100271, 'seqStr', "sometime after eight", 'seed', "AE"); % two activations
                    % uu{end+1} = struct('clusId', 520200153, 'seqStr', "pulled my leg", 'seed', "AY");
                    % uu{end+1} = struct('clusId', 450100004, 'seqStr', "in love with", 'seed', "AH"); % a bit noisy in stim
                    % uu{end+1} = struct('clusId', 460100098, 'seqStr', "earth's the matter", 'seed', "DH");
                    % uu{end+1} = struct('clusId', 430100161, 'seqStr', "me by surprise", 'seed', "AY"); % from STG
                case 'bridge'
                    uu{end+1} = struct('clusId', 410100463, 'seqStr', "you got enough", 'seed', "G"); % similar to 410100155
                    uu{end+1} = struct('clusId', 450100286, 'seqStr', "me by surprise", 'seed', "B");
                    uu{end+1} = struct('clusId', 410100155, 'seqStr', "have you got", 'seed', "Y");
                    uu{end+1} = struct('clusId', 410100454, 'seqStr', "something pulled", 'seed', "NG"); % tuned position is "UH"
                    uu{end+1} = struct('clusId', 450100165, 'seqStr', "got plenty of", 'seed', "N");
                    uu{end+1} = struct('clusId', 410100450, 'seqStr', "they even pay", 'seed', "V"); % mixture of feedback and bridge
                    % uu{end+1} = struct('clusId', 410100441, 'seqStr', "we've got", 'seed', "G");
                case 'feedback'
                    uu{end+1} = struct('clusId', 440300585, 'seqStr', "search this house", 'seed', "S");
                    uu{end+1} = struct('clusId', 500300452, 'seqStr', "girl nodded understandingly", 'seed', "N");
                    uu{end+1} = struct('clusId', 430100549, 'seqStr', "of time to", 'seed', "AY");
                    uu{end+1} = struct('clusId', 430100582, 'seqStr', "enough blankets", 'seed', "K");
                    uu{end+1} = struct('clusId', 440200479, 'seqStr', "it was nobody's", 'seed', "Z"); % tuned to "OW"
                    % uu{end+1} = struct('clusId', 430100210, 'seqStr', "me by surprise", 'seed', "AY");
                case 'mirror_movie'
                    uu{end+1} = struct('clusId', 460100003, 'seqStr', "earth's the matter", 'seed', "DH");
                    % uu{end+1} = struct('clusId', 460100271, 'seqStr', "sometime after eight", 'seed', "AE"); % two activations
                    uu{end+1} = struct('clusId', 520200153, 'seqStr', "will you tell", 'seed', "Y");
            end
        end
        
        function [clusTb, clusTbPath] = LoadClusTb(mdlName)
            % Load cached linker clusTb
            % 
            %   [clusTb, clusTbPath] = LoadClusTb()
            %   [clusTb, clusTbPath] = LoadClusTb(mdlName)
            % 
            if nargin < 1
                mdlName = "smooth_lm";
            end
            clusTbPath = fullfile(LMV.Data.GetAnalysisDir, "linker", mdlName, "linker_clusTb.mat");
            load(clusTbPath, 'clusTb');
        end
        
        % Unit responses
        function cTb = ComputeSentenceResponseFromCache(clus)
            % Compute M2 sentence responses for each unit
            % 
            %   cTb = ComputeSentenceResponseFromCache(clus)
            % 
            
            % Initialize cluster table
            if ~istable(clus)
                cTb = table;
                cTb.clusId = clus;
            else
                cTb = clus;
            end
            
            % Define resampling parameters
            % tStart = cTb.uCache{1}.tt.stimMatchOn(1);
            tStart = 0;
            tEnd = cTb.uCache{1}.tt.prodMatchOff(1);
            rsBinSize = 0.02;
            rsKerSize = 0.05;
            rsPad = [-1 1];
            rsEdges = MMath.BinCenters2Edges(tStart+rsPad(1) : rsBinSize : tEnd+rsPad(2));
            
            plotPad = [-1 1]*0.3;
            t = MMath.BinEdges2Centers(rsEdges)';
            isShow = t > tStart+plotPad(1) & t < tEnd+plotPad(2);
            t = t(isShow);
            
            % Iterate through units
            for i = 1 : height(cTb)
                fprintf("%i/%i\tu%i\n", i, height(cTb), cTb.clusId(i));
                s = cTb.uCache{i};
                
                % Exclude sentences with less than 3 repeats
                nRep = cellfun(@numel, s.st);
                isSen = nRep > 2;
                nRep = nRep(isSen);
                stCell = s.st(isSen);
                tt = s.tt(isSen,:);
                tv = s.tv(isSen,:);
                
                % Compute spike rates of all trials
                stAll = cat(1, stCell{:});
                srAll = cellfun(@(x) MNeuro.Filter1(x, [], 'bin', rsEdges), stAll, 'Uni', false);
                srAll = cellfun(@(x) MNeuro.Filter1(x, 1/rsBinSize, 'gaussian', rsKerSize)', srAll, 'Uni', false);
                srAll = cat(1, srAll{:});
                
                nTrials = cellfun(@numel, stCell);
                srSen = mat2cell(srAll, nTrials);
                
                % Iterate through sentences
                nSen = numel(stCell);
                meanCell = cell(nSen,1);
                sigCell = cell(nSen,1);
                for j = 1 : nSen
                    % Compute sentence mean spike rates
                    meanCell{j} = mean(srSen{j}(:,isShow))';
                    
                    % Bootstrap mean spike rates
                    nBoot = 100;
                    M = cell(1,nBoot);
                    for k = 1 : nBoot
                        ind = randsample(size(srAll,1), size(srSen{j},1));
                        M{k} = mean(srAll(ind, isShow))';
                    end
                    M = cat(2, M{:});
                    sigCell{j} = meanCell{j} > prctile(M, 95, 2);
                end
                meanResp = cat(2, meanCell{:});
                isSigResp = cat(2, sigCell{:});
                
                sigResp = meanResp;
                sigResp(~isSigResp) = NaN;
                
                % Check significant activation in three phases
                isStim = tt.stim(1).MaskTimestamps(t);
                isProd = tt.prod{1}.MaskTimestamps(t);
                isDelay = t > tt.stim(1).T.tmax & t < tt.cue3(1).T.tOn;
                isInit = t > tt.cue3(1).T.tOn & t < tt.prod{1}.T.tmin;
                durSig = [ ...
                    sum(isStim & isSigResp); ...
                    sum(isDelay & isSigResp); ...
                    sum(isInit & isSigResp); ...
                    sum(isProd & isSigResp); ...
                    ] * rsBinSize;
                isPhase = durSig > 0.2;
                
                mirrorResp = sigResp;
                mirrorResp(:, ~all(isPhase([1 4],:),1)) = NaN;
                
                bridgeResp = sigResp;
                bridgeResp(:, ~all(isPhase,1)) = NaN;
                
                cTb.stimId{i} = tv.stimId;
                cTb.stim{i} = tt.stim;
                cTb.prod{i} = cellfun(@(x) x(1), tt.prod);
                cTb.nRep{i} = nRep;
                cTb.time{i} = t;
                cTb.resp{i} = meanResp;
                cTb.sigResp{i} = sigResp;
                cTb.sigPhaseDur{i} = durSig;
                cTb.mirrorResp{i} = mirrorResp;
                cTb.bridgeResp{i} = bridgeResp;
            end
            
        end
        
        % Time morphing
        function seMorph = MorphSession(se)
            % 
            %   seMorph = MorphSession(se)
            % 
            
            % Operate on hard copy
            se = se.Duplicate;
            
            % Match prod to stim
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            rt = se.GetReferenceTime('taskTime');
            [tt, tv] = LMV.SE.MatchProd2Stim(tt, tv, rt);
            se.SetTable('taskTime', tt);
            
            % Cut sentences to phones
            %   direct cutting may leave sentence objects into wrong adjacent epochs
            if ~iscell(tt.prod)
                tt.prod = num2cell(tt.prod);
            end
            for i = 1 : height(tt)
                tt.stimPhn{i} = Cut(Cut(tt.stim(i)));
                tt.prodPhn{i} = Cut(Cut(tt.prod{i}));
            end
            se.SetTable('taskTime', tt);
            
            % Reslice se by stim and prod phases
            tCue = sort([tt.cue1On+rt; tt.cue3On+rt]);
            se.SliceSession(tCue, 'absolute');
            se.SetColumn('taskTime', 'cueOn', zeros(se.numEpochs, 1));
            
            % Add expanded taskValue table
            tv2 = tv(repelem(1:height(tv), 2), :);
            tv2.phase = repmat(["stim"; "prod"], [height(tv) 1]);
            tv2.alignScore(tv2.phase=="stim") = Inf;
            tv2.tReact(tv2.phase=="stim") = 0;
            se.SetTable('taskValue', tv2, 'eventValues');
            
            % Find morph times
            tt2 = se.GetTable('taskTime');
            for i = 1 : 2 : height(tt2)
                tt2.combPhn{i} = tt2.stimPhn{i};
                tt2.combPhn{i+1} = tt2.prodPhn{i+1};
            end
            [tt2, tv2] = LMV.Linker.FindMorphTimes(tt2, tv2, [], 'best');
            se.SetTable('taskTime', tt2);
            se.SetTable('taskValue', tv2);
            
            % Morph se
            se.SetColumn('taskTime', 'trialOn', zeros(se.numEpochs, 1));
            seMorph = NP.SE.MorphSession(se);
            
            % Add prodMatchOn, prodMatchOff to stim epochs
            tt2 = seMorph.GetTable('taskTime');
            for i = 1 : 2 : height(tt2)
                tt2.prodMatchOn(i) = tt2.prodMatchOn(i+1);
                tt2.prodMatchOff(i) = tt2.prodMatchOff(i+1);
            end
            seMorph.SetTable('taskTime', tt2);
        end
        
        function [tt, tv] = FindMorphTimes(tt, tv, rt, targetType)
            % Find matching key times between source and temaplte trials
            % 
            %   [tt, tv] = FindMorphTimes(tt, tv)
            %   [tt, tv] = FindMorphTimes(tt, tv, rt)
            %   [tt, tv] = FindMorphTimes(tt, tv, rt, targetType)
            % 
            
            if nargin < 4
                targetType = 'median';
            end
            targetType = lower(targetType);
            
            % Preallocate columns
            tt.morphFrom = num2cell(zeros(size(tt.cue1On)));
            tt.morphTo = num2cell(zeros(size(tt.cue1On)));
            tv.tempTrialNum = NaN(size(tt.cue1On));
            
            % Select phones for each epoch
            for i = height(tt) : -1 : 1
                phn = tt.combPhn{i};
                phn = phn(phn > tt.cueOn(i)); % exclude events before cueOn or are NaN
                if exist('rt', 'var') && ~isempty(rt) && i < height(tt)
                    phn = phn(phn+rt(i) < tt.cueOn(i+1)+rt(i+1)); % exclude events that bleed into the next trial
                end
                phones{i,1} = phn;
                hasPhn(i,1) = ~isempty(phn);
            end
            
            % Find unique stimuli
            [stimList, ~, trialStimId] = unique(tv.stimText);
            
            for s = 1 : numel(stimList)
                % Find non-empty trials of the same sentence
                trialInd = find(trialStimId==s & hasPhn);
                
                % Find the trial(s) that has the highest alignment score
                sc = tv.alignScore(trialInd);
                bestTrials = trialInd(sc == max(sc));
                
                % Break ties by choosing the trial with reaction time closest to median
                tReact = tv.tReact(bestTrials);
                [~, I] = min(abs(tReact - median(tReact)));
                bestTrial = bestTrials(I);
                bestRep = trialInd == bestTrial;
                
                % Get the template speech events
                tp = phones{bestTrial};
                assert(~isempty(tp), "Template for '%s' is empty", stimList(s));
                
                % Align repeats and collect all the aligned times
                preEvt = {'cueOn'}; % events preceeding speech
                nPre = numel(preEvt);
                T = NaN(numel(trialInd), nPre+numel(tp));
                for j = 1 : numel(trialInd)
                    % Get the source speech events to be morphed
                    k = trialInd(j);
                    src = phones{k};
                    if isempty(src)
                        continue
                    end
                    
                    % Find timestamps of matched words
                    [~, tSrc, info] = NP.TaskBaseClass.FindMatchedSpeechTimes(tp, src);
                    
                    % Store full source event times
                    kTp = [1:nPre nPre+info.kk(1,:)];
                    tSrc = [tt{k,preEvt}'; tSrc];
                    T(j,kTp) = tSrc;
                    
                    % Normalize alignment score to max score
                    info.nscore = info.tscore / numel(tp);
                    
                    tv.alignInfo{k} = info;
                    tv.tempTrialNum(k) = tv.trialNum(bestTrial);
                    tt.morphFrom{k} = tSrc;
                end
                
                % Determine template times
                itvl = diff(T, 1, 2);
                if targetType == "median"
                    % Use median event intervals
                    medItvl = median(itvl, 1, 'omitnan');
                    medItvl(4) = 1; % force the mental holding interval (4th) to the designed mean of 1 sec
                    tTp = cumsum([tt.cueOn(bestTrial) medItvl])';
                elseif targetType == "best"
                    % Use the best aligned trial
                    tTp = cumsum([tt.cueOn(bestTrial) itvl(bestRep,:)])';
                end
                
                % Find and save matching times from the template
                for j = 1 : numel(trialInd)
                    k = trialInd(j);
                    info = tv.alignInfo{k};
                    if isempty(info)
                        continue
                    end
                    kTp = [1:nPre nPre+info.kk(1,:)];
                    tt.morphTo{k} = tTp(kTp);
                end
            end
        end
        
        function s = ComputeMelXC(se, ops)
            % Compute cross-correlation between the Mel spectrogram power of stim and prod
            % 
            %   s = ComputeMelXC(se, ops)
            % 
            
            % Operate on a hard copy
            se = se.Duplicate;
            
            % Remove orphan epochs
            tv = se.GetTable('taskValue');
            tv.nEp = arrayfun(@(x) sum(x==tv.trialNum), tv.trialNum);
            se.RemoveEpochs(tv.nEp < 2);
            
            % Get sampling times
            [tt, tv] = se.GetTable('taskTime', 'taskValue');
            isStim = tv.phase == "stim";
            sEdges = arrayfun(@(a,b) a : ops.rsBinSize : b, tt.stimOn(isStim), tt.stimOff(isStim), 'Uni', false);
            
            % Resample stim spectrogram
            ms = se.ResampleTimeSeries('mel', sEdges, isStim, {'speaker1'});
            ms = cell2mat(ms.speaker1);
            ms = sum(ms, 2); % use total power
            
            % Loop through time shifts
            dt = ops.rsShifts(:)';
            xr = NaN(size(dt));
            pval = xr;
            for t = 1 : numel(dt)
                % Resample prod spectrogram
                pEdges = cellfun(@(x) x+dt(t), sEdges, 'Uni', false);
                mp = se.ResampleTimeSeries('mel', pEdges, ~isStim, {'mic'});
                mp = cell2mat(mp.mic);
                mp = sum(mp, 2);
                
                % Compute cross-correlation
                [xr(t), pval(t)] = corr(ms, mp);
            end
            
            % Find lags
            [pkCoef, I] = max(xr);
            pkPval = pval(I);
            pkLag = dt(I);
            
            % Output
            s.xrCoef = xr;
            s.xrPval = pval;
            s.xrLag = dt;
            s.xrPkCoef = pkCoef;
            s.xrPkPval = pkPval;
            s.xrPkLag = pkLag;
        end
        
        function PlotMelAlignment(seTb, ce)
            % Plot Mel spectrograms of all trials 
            
            stimList = unique(seTb.stimText);
            nSen = numel(stimList);
            
            rowDist = [1 2 1 2];
            nCols = nSen;
            t = tiledlayout(sum(rowDist), nSen);
            t.Padding = 'compact';
            
            for k = 1 : nSen
                % Stim alignment
                isSE = seTb.stimText == stimList(k) & seTb.phase == "stim";
                if ~any(isSE)
                    continue
                end
                seStim = seTb.se(isSE).Duplicate;
                seStim.userData = ce.userData;
                tt = seStim.GetTable('taskTime');
                t0 = tt.stimOn(1);
                seStim.AlignTime(t0);
                
                tt = seStim.GetTable('taskTime');
                tWin = [tt.stimOn(1) tt.stimOff(1)] + [-1 1]*0.2;
                
                ntArgs = MPlot.FindTileInd(rowDist, nCols, 1, k);
                ax = nexttile(ntArgs{:});
                NP.TaskBaseClass.PlotAlignedWords(ax, seStim, [], tWin, 'combPhn');
                ax.XLabel.String = [];
                ax.YLabel.String = "Trials";
                ax.Title.String = sprintf("%s (stim)", seTb.stimText(isSE));
                
                ntArgs = MPlot.FindTileInd(rowDist, nCols, 2, k);
                ax = nexttile(ntArgs{:});
                NP.TaskBaseClass.PlotMelSpectrograms(ax, seStim, [], tWin, 'VariableName', 'speaker1', 'Waveform', false);
                ax.XLabel.String = [];
                
                % Production alignment
                isSE = seTb.stimText == stimList(k) & seTb.phase == "prod";
                if ~any(isSE)
                    continue
                end
                seProd = seTb.se(isSE).Duplicate;
                seProd.userData = ce.userData;
                seProd.AlignTime(t0);
                
                ntArgs = MPlot.FindTileInd(rowDist, nCols, 3, k);
                ax = nexttile(ntArgs{:});
                NP.TaskBaseClass.PlotAlignedWords(ax, seProd, [], tWin, 'combPhn');
                ax.XLabel.String = [];
                ax.YLabel.String = "Trials";
                ax.Title.String = sprintf("%s (prod)", seTb.stimText(isSE));
                
                ntArgs = MPlot.FindTileInd(rowDist, nCols, 4, k);
                ax = nexttile(ntArgs{:});
                NP.TaskBaseClass.PlotMelSpectrograms(ax, seProd, [], tWin, 'VariableName', 'mic', 'Waveform', false);
                ax.XLabel.String = "Aligned time (s)";
            end
        end
        
        function PlotOverlay(ce, uInd, stimIdList)
            % Plot overlays of aligned stim and prod responses
            % 
            %   PlotOverlay(ce, uInd, stimIdList)
            % 
            
            if nargin < 3
                stimIdList = LMV.Param.stimIdList4;
            end
            
            [tt, tv] = ce.GetTable('taskTime', 'taskValue');
            clusTb = NP.Unit.GetClusTb(ce);
            clusTb.clusId = clusTb.clusId - NP.Unit.GetBaseClusId(clusTb.recId);
            
            for u = uInd(:)'
                if isnan(u)
                    continue
                end
                ax = nexttile;
                
                m = tv.phase == "stim" & ismember(tv.stimId, stimIdList);
                [ts, rs] = ce.GetArray('resp', m, u+1, 'DimCat', 0);
                m = tv.phase == "prod" & ismember(tv.stimId, stimIdList);
                [tp, rp] = ce.GetArray('resp', m, u+1, 'DimCat', 0);
                
                % Apply optimal time offset
                tp = cellfun(@(x) x-clusTb.xrPkLag(u), tp, 'Uni', false);
                
                rr = cell2mat([rs; rp]);
                rrMax = max(rr)*2;
                if rrMax > 0
                    rs = cellfun(@(x) x/rrMax, rs, 'Uni', false);
                    rp = cellfun(@(x) x/rrMax, rp, 'Uni', false);
                end
                y = (1 : numel(rs)) + 0.2;
                
                MPlot.PlotTraceLadder(ts, rs, y, 'Color', 'b'); hold on
                MPlot.PlotTraceLadder(tp, rp, y, 'Color', 'r');
                
                % Phonetic markers
                m = tv.phase == "stim" & ismember(tv.stimId, stimIdList);
                m = find(m);
                for i = 1 : numel(m)
                    tge = tt.stim(m(i));
                    text(ax, 0, i, tge.GetParentLabel, 'HorizontalAlignment', 'left');
                end
                
                ax.Box = 'off';
                ax.YTick = 1:numel(rs);
                ax.XLim = [-.2 max(cellfun(@(x)x(end), ts))];
                ax.YLim = [1 numel(rs)+1] - 0.2;
                ax.Title.String = sprintf("u%i, dt = %.2fs, r = %.2f", clusTb.clusId(u), clusTb.xrPkLag(u), clusTb.xrPkCoef(u));
                MPlot.Axes(ax);
            end
        end
        
        % Sequence PETH overlay
        function s = LoadSeqData(recId)
            % Load cached sequence data
            % 
            %   s = LoadSeqData(recId)
            % 
            
            tWin = [-0.4 0.4];
            
            % Load cached sequence data
            cacheDir = fullfile(LMV.Data.GetAnalysisDir, "linker", "extracted_seq");
            s = load(fullfile(cacheDir, recId+"_seqData.mat"), "se", "seeds", "seqData");
            se = s.se;
            clusTb = NP.Unit.GetClusTb(se);
            cols = ["stim", "prod", "prod"];
            seqData = s.seqData(:,cols);
            
            % Load phone ti-RFs
            targets = ["stim", "feedback", "prod"];
            trfTb = LMV.TRF.LoadModels("phone_"+targets, recId);
            [~, trfInd] = MMath.SortLike(trfTb.clusId, clusTb.clusId);
            trfTb = trfTb(trfInd,:);
            
            % Add PETHs to seqTb
            for tIdx = 1 : numel(targets)
                C = seqData.(tIdx){"word"}{:,:};
                
                % tiRF models will be used to provide the optimal time to extract response magnitude
                mdls = trfTb.("phone_"+targets(tIdx));
                mdls = LMV.RF.UpdateModelR2(mdls);
                
                for k = 1 : numel(C)
                    seqTb = C{k}.seqTb;
                    if isempty(seqTb)
                        continue
                    end
                    seqTb(seqTb.nSample<3,:) = [];
                    seqTb = LMV.PE.AddResp2SeqTb(seqTb, se, tWin, [], mdls);
                    C{k}.seqTb = seqTb;
                end
                
                seqData.(tIdx){"word"}{:,:} = C;
            end
            
            % Denest the tables
            C = seqData{:,:};
            for i = 1 : numel(C)
                C{i}.Properties.VariableNames = targets(i);
            end
            s.triTb = cat(2, C{:});
        end
        
        function PlotSeqPethOverlayOnBrush(fig, axesStruct)
            % 
            %   PlotSeqPethOverlayOnBrush(fig, axesStruct)
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
            b = find(b);
            if numel(b) > 10 % limit to 10 sequences
                b = b(1:10);
            end
            
            % Get inputs
            ud = axesStruct.Axes.UserData;
            triTb = ud.triTb;
            uIdx = ud.uIdx;
            maxR = ud.maxR;
            seqStr = ud.seqStr(b);
            selectedSeeds = ud.seed(b);
            isSelected = ismember(triTb.Properties.RowNames, selectedSeeds);
            
            % Plotting
            if isempty(fig.UserData) || ~isvalid(fig.UserData.cbFig)
                % Create figure if absent
                fig.UserData.cbFig = MPlot.Figure( ...
                    'Name', 'Selected Units', ...
                    'NumberTitle', 'off', ...
                    'Menubar', 'none', ...
                    'Toolbar', 'figure', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
            end
            f = fig.UserData.cbFig;
            figure(f);
            clf(f);
            f.Position(3) = 250 * numel(unique(selectedSeeds));
            LMV.Linker.PlotPeriEventSeqArray(triTb(isSelected,:), 'SeqStr', seqStr, 'UnitIdx', uIdx, 'MaxSpikeRate', maxR);
        end
        
        function PlotPeriEventSeqArray(triTb, varargin)
            % Plot an array of panels each contains rasters and PETHs separated by different speech sequences
            % 
            %   PlotPeriEventSeqArray(triTb)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('UnitIdx', 1, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('MaxSpikeRate', [], @(x) isscalar(x) && isnumeric(x));
            p.addParameter('Parent', [], @(x) true);
            p.parse(varargin{:});
            uIdx = p.Results.UnitIdx;
            maxR = p.Results.MaxSpikeRate;
            axObjs = p.Results.Parent;
            
            if ismember('score', triTb.Properties.VariableNames)
                triTb.score = num2cell(triTb.score);
            end
            
            C = triTb{:,1:3};
            
            if isempty(maxR)
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
            nRow = 1;
            nCol = height(C);
            
            if isempty(axObjs)
                tl = tiledlayout(nRow, nCol);
                tl.Padding = "compact";
            end
            
            for i = 1 : nCol
                if isempty(axObjs)
                    ax = nexttile;
                else
                    ax = axObjs(i);
                end
                LMV.Linker.PlotSeqPethOverlay(ax, triTb{i,:}, 'MaxSpikeRate', maxR, 'UnitIdx', uIdx, p.Unmatched);
            end
        end
        
        function PlotSeqPethOverlay(ax, C, varargin)
            % Plot PETHs separated by different speech sequences
            % 
            %   PlotSeqPethOverlay(ax, C)
            % 
            
            p = inputParser;
            p.addParameter('SeqStr', [], @(x) isstring(x) || iscellstr(x) || ischar(x) || isempty(x));
            p.addParameter('UnitIdx', 1, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('MaxSpikeRate', [], @(x) isscalar(x) && isnumeric(x));
            p.parse(varargin{:});
            seqStr = string(p.Results.SeqStr);
            uIdx = p.Results.UnitIdx;
            maxR = p.Results.MaxSpikeRate;
            
            cc = lines(numel(C));
            seqTbs = cell(1,3);
            for i = 1 : 3
                % Organize table
                s = C{i};
                tb = s.seqTb;
                if ~isempty(seqStr)
                    m = ismember(tb.seqStr, seqStr);
                    tb = tb(m,:);
                end
                if isempty(tb)
                    return
                end
                seqTbs{i} = tb;
                
                if i == 3
                    continue
                end
                
                % PETH
                t = tb.tResp{1};
                hh = cat(3, tb.hh{:});
                hh = squeeze(hh(:,uIdx,:));
                ee = cat(3, tb.ee{:});
                ee = squeeze(ee(:,uIdx,:));
                if isempty(maxR)
                    maxR = max(hh(:));
                end
                hh = hh ./ maxR;
                ee = ee ./ maxR;
                
                MPlot.PlotHistStack(t, hh, ee, 'Scaling', 0.6, 'Color', cc(i,:), 'Parent', ax);
            end
            
            % Labels
            tWin = [-.4 .4];
            for i = 1 : height(tb)
                tge = tb.seqTge{i};
                isPlot = tge > tWin(1) & tge < tWin(2);
                tge(isPlot).Plot(i-0.5+[0 .2], 'Parent', ax);
                hold(ax, 'on');
            end
            
            if ismember('xrPkVal', seqTbs{2}.stats{1}.Properties.VariableNames)
                xc = round(cellfun(@(x) x.xrPkVal(uIdx), seqTbs{2}.stats));
                ixc = round(cellfun(@(x) x.xrPkVal(uIdx), seqTbs{3}.stats));
                labelArray = string(["XC "+xc'; "iXC "+ixc']);
                labelArray(ismissing(labelArray)) = "";
                tickLabels = sprintf('%s\\newline%s\n', labelArray(:));
                ax.YTickLabel = strtrim(tickLabels);
            end
            
            % 
            if isfield(C{1}, 'trfMdl')
                mdl = C{1}.trfMdl;
                ax.Title.String = sprintf("%s, score=%.2g", mdl.resps, C{4});
            else
                ax.Title.String = sprintf("/%s/ in %i-%s", MLing.ARPA2IPA(s.seed), s.nElem, s.level);
            end
            
            ax.XLim = tWin;
            MPlot.Axes(ax);
        end
        
        % Encoding models
        function [mdls, mRF] = MatchRF(mdls, varargin)
            % 
            % 
            %   [mdls, mRF] = MatchRF(mdls)
            % 
            
            p = inputParser;
            p.addParameter('TimeWindow', [0 .3], @isnumeric);
            p.parse(varargin{:});
            tWin = sort(abs(p.Results.TimeWindow));
            
            BB = cellfun(@(x) reshape(x.Beta, numel(x.feats), numel(x.dt)), mdls, 'Uni', false);
            % r = corr(BB{1}, BB{2});
            for i = size(BB{1},2) : -1 : 1
                for j = size(BB{2},2) : -1 : 1
                    C = cov(BB{1}(:,i), BB{2}(:,j));
                    r(i,j) = C(1,2);
                end
            end
            r = smoothdata2(r, "movmean", 5);
            
            for i = 1 : numel(mdls)
                mdl = mdls{i};
                if contains(mdl.name, "prod")
                    isWin = mdl.dt > -tWin(2) & mdl.dt < tWin(1);
                else
                    isWin = mdl.dt > tWin(1) & mdl.dt < tWin(2);
                end
                if i == 1
                    r(~isWin,:) = NaN;
                elseif i == 2
                    r(:,~isWin) = NaN;
                end
            end
            
            [rMax, I] = max(r, [], "all");
            [i, j] = ind2sub(size(r), I);
            ind = [i j];
            mRF = cell(size(mdls));
            for k = 1 : numel(mdls)
                mdls{k}.r2Idx = ind(k);
                
                stRF = mdls{k};
                RF = stRF.mdls{ind(k)};
                RF.feats = stRF.feats;
                RF.resps = stRF.resps;
                mRF{k} = RF;
            end
        end
        
    end
    
end