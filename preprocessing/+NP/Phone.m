classdef Phone
    %PHONE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        % Categories of phonemes
        vowels = {'AE','AX', 'AXR','AA','AH','AO','AW','AY','EH','ER','EY','IH','IY', 'IX','OW','OY','UH','UW', 'UX'};
        high = {'IY', 'IH', 'IX', 'UH', 'UW', 'UX'};
        mid = {'EH', 'EY', 'AX', 'AXR','ER', 'AH', 'AO', 'OW', 'OY'};
        low = {'AE', 'AA', 'AY', 'AW'};
        front = {'IH', 'IY', 'IX', 'EH', 'EY', 'AE', 'AX', 'AXR','ER'};
        back = {'UH', 'UW', 'UX', 'OW', 'AO', 'OY', 'AA', 'AW', 'AY', 'AH'};
        diphthongs = {'EY', 'AY', 'OW', 'OY', 'AW'};
        
        consonants = [NP.Phone.plosives NP.Phone.fricatives NP.Phone.affricates NP.Phone.nasals NP.Phone.approximants];
        plosives = {'P', 'B', 'T', 'D', 'K', 'G'};
        fricatives = {'F', 'V', 'TH', 'DH', 'S', 'Z', 'SH', 'ZH', 'HH'};
        affricates = {'CH', 'JH'};
        nasals = {'M', 'N', 'NG'};
        approximants = {'L', 'R', 'W', 'Y'};
        
        labial = {'P', 'B', 'M', 'W', 'F', 'V'};
        coronal = [NP.Phone.dental NP.Phone.alveolar NP.Phone.postalveolar];
        dental = {'TH', 'DH'};
        alveolar = {'T', 'D', 'S', 'Z', 'N', 'R', 'L'};
        postalveolar = {'SH', 'ZH', 'CH', 'JH'};
        palatal = {'Y'};
        velar = {'K', 'G', 'NG'};
        glottal = {'HH'};
        
        rounded = {'OW','OY','AW','AY','UH','UW','AO','R','W'}; % needs to double-check
        
        voiced = [ ...
            NP.Phone.vowels ...
            {'B','D','G'} ...           % voiced plosives
            {'V','DH','Z','ZH'} ...     % voiced fricatives
            {'JH'} ...                  % voiced affricate
            NP.Phone.nasals NP.Phone.approximants ...
            ];
    end
    
    methods(Static)
        function AddTimitGroundTruth(se)
            % 
            
            % Read TIMIT data
            srcDir = fullfile(NP.Data.GetProjectRoot, "code", "tasks", "LMV", "TIMIT");
            
            tv = se.GetTable('taskValue');
            timitIDs = unique(tv.stimId);
            
            [W, T] = MLing.ReadTimitWaveform(srcDir, timitIDs);
            tg = MLing.ReadTimitFeatures(srcDir, timitIDs);
            
            % Convert tg structs to NP.TGEvent objects
            tg = NP.TGEvent(tg);
            tg = tg.SetVfield('timitId', timitIDs);
            
            % Add data to se
            [task, mel] = se.GetTable('taskTime', 'mel');
            for k = 1 : se.numEpochs
                % Find TIMIT stim for this trial
                stimIdx = tv.stimId(k) == timitIDs;
                if ~any(stimIdx)
                    continue
                end
                w = W{stimIdx};
                t = T{stimIdx};
                
                % Find time offset using spectrogram crosscorrelation between speaker and GT audio
                tNI = mel.time{k};
                sNI = mel.speaker1{k};
                [sGT, ~, tGT] = NP.Audio.ComputeMelSpectrogram(w, t);
                sGT = sGT';
                fs = 1 / mean(diff(tNI));
                
                pGT = sum(sGT,2);
                pGT = pGT - min(pGT);
                pGT = pGT / max(pGT);
                pNI = sum(sNI,2);
                pNI = pNI - min(pNI);
                pNI = pNI / max(pNI);
                
                [r, lags] = xcorr(pNI, pGT);
                [~, I] = max(r);
                dt = lags(I) / fs;
                
                % Add NP.TGEvent objects to taskTime table
                task.stimGT(k) = tg(stimIdx) + dt;
                
                % Resample GT spectrogram to fit in mel table
                fillVal = min(sGT(:));
                sGTq = interp1(tGT+dt, sGT, tNI, 'linear', fillVal);
                mel.timit{k} = sGTq;
            end
            
            se.SetTable('taskTime', task);
            se.SetTable('mel', mel);
        end
        
        function AddPhoneEvents(se, phGroups)
            % Add phoneme NP.TGEvent objects table named 'phone' to se
            % 
            %   AddPhoneEvents(se, phGroups)
            % 
            % Inputs
            %   se              MSessionsExplore object with a single epoch.
            %   phGroups        1) Phoneme names in a cell array vector.
            %                   2) A n-by-2 cell array. The first column contains group names. 
            %                   Each element in the second column is a cell array of group memebers.
            % 
            
            if isempty(phGroups)
                phGroups = { ...
                    'high',         NP.Phone.high; ...
                    'mid',          NP.Phone.mid; ...
                    'low',          NP.Phone.low; ...
                    'front',        NP.Phone.front; ...
                    'back',         NP.Phone.back; ...
                    'rounded',      NP.Phone.rounded; ...
                    'plosives',     NP.Phone.plosives; ...
                    'fricatives',   NP.Phone.fricatives; ...
                    'nasals',       NP.Phone.nasals; ...
                    'approximants', NP.Phone.approximants; ...
                    'labial',       NP.Phone.labial; ...
                    'velar',        NP.Phone.velar; ...
                    'coronal',      NP.Phone.coronal; ...
                    'glottal',      NP.Phone.glottal; ...
                    'dental',       NP.Phone.dental; ...
                    };
            end
            
            if isvector(phGroups)
                phGroups = repmat(phGroups(:), [1 2]);
            end
            
            tt = se.GetTable('taskTime');
            if ~iscell(tt.stim)
                tt.stim = num2cell(tt.stim);
            end
            if ~iscell(tt.prod)
                tt.prod = num2cell(tt.prod);
            end
            
            ph = table;
            for i = height(tt) : -1 : 1
                % Get phonemes
                sen = cat(1, tt.stim{i}, tt.prod{i});
                tge = Cut(Cut(sen));
                lb = tge.GetParentLabel();
                lb = erase(lb, digitsPattern);
                
                % Get masks
                for j = 1 : size(phGroups,1)
                    gn = phGroups{j,1};
                    gm = cellstr(phGroups{j,2});
                    isGroup = ismember(lb, gm);
                    ph.(gn){i} = tge(isGroup);
                end
            end
            se.SetTable('phone', ph, 'eventTimes', se.GetReferenceTime('taskTime'));
        end
        
        function AddPhoneSeqEvents(se, nElem, seeds)
            % Add phoneme NP.TGEvent objects table named 'phone' to se
            % 
            %   AddPhoneSeqEvents(se)
            %   AddPhoneSeqEvents(se, nElem)
            %   AddPhoneSeqEvents(se, nElem, seeds)
            % 
            % Inputs
            %   se              MSessionsExplore object with a single epoch.
            %   seeds           Phoneme names in a vector of cell array.
            % 
            
            seeds = [NP.Phone.high, NP.Phone.mid, NP.Phone.low, NP.Phone.consonants]';
            nElem = 2;
            
            tt = se.GetTable('taskTime');
            if ~iscell(tt.stim)
                tt.stim = num2cell(tt.stim);
            end
            if ~iscell(tt.prod)
                tt.prod = num2cell(tt.prod);
            end
            
            seq = table;
            for i = height(tt) : -1 : 1
                % Get phonemes
                sen = cat(1, tt.stim{i}, tt.prod{i});
                tge = Cut(Cut(sen));
                lb = tge.GetParentLabel();
                lb = erase(lb, digitsPattern);
                
                % Get masks
                for j = 1 : size(seqGroups,1)
                    gn = seqGroups{j,1};
                    gm = cellstr(seqGroups{j,2});
                    isGroup = ismember(lb, gm);
                    seq.(gn){i} = tge(isGroup);
                end
            end
            se.SetTable("phone"+nElem, seq, 'eventTimes', se.GetReferenceTime('taskTime'));
        end
        
        function tb = PoolPhonemes(se, varNames)
            % Make a table by pooling all the phonemes from the specified event variables in se
            % 
            %   tb = PoolPhonemes(se, varNames)
            % 
            
            % Get phoneme events
            sen = se.GetColumn('taskTime', varNames);
            sen = cat(1, sen{:});
            phn = Cut(Cut(sen));
            
            % Standardize phonemic labels
            arp = erase(phn.GetParentLabel, digitsPattern);
            arp = NP.Phone.Standardize(arp);
            
            % Compute mean artic for each phoneme
            artic = se.GetTable('artic');
            t = artic.time{1};
            a = cell2mat(artic{1,2:end});
            aMean = NaN(numel(phn), size(a,2));
            for i = 1 : numel(phn)
                m = phn(i).T.tmin <= t & t < phn(i).T.tmax;
                aMean(i,:) = mean(a(m,:), 1, 'omitnan');
            end
            
            % Make a table
            tb = table;
            tb.phe = phn;
            tb.arp = arp;
            tb.cat = categorical(arp);
            tb.ipa = MLing.ARPA2IPA(arp);
            tb.artic = aMean;
        end
        
        function s = ComputeLDA(phTb, catName)
            % 
            
            switch catName
                case 'vowel'
                    phList = {'AA', 'AE', 'AH', 'UH', 'UW', 'IH', 'IY'};
                    % phList = {'AA', 'AE', 'AH', 'UH', 'IH', 'IY'};
                    
                    % phDisc = ["IY", "UW"; "IY", "AE"];
                    phDisc = ["UW", "IY"; "UW", "AE"];
                    % phDisc = ["UH", "IY"; "UH", "AE"];
                    
                case 'consonant'
                    phList = {'M', 'B', 'P', 'V', 'F', 'SH', 'Z', 'S', 'N', 'T', 'D', 'NG', 'K', 'G'};
                    phDisc = ["B", "S"; "B", "K"];
                    
                otherwise
                    error("catName must be 'vowel', or 'consonant'");
            end
            
            % Select phonemes
            isSelect = ismember(phTb.arp, phList);
            phTb = phTb(isSelect, :);
            phTb = sortrows(phTb, 'arp');
            
            % Assign sample weights inverse to phoneme probability
            for i = 1 : numel(phList)
                isPh = phTb.arp == phList{i};
                phTb.weight(isPh) = 1 / sum(isPh);
            end
            
            % Fit LDA
            nVars = size(phTb.artic, 2);
            if nVars == 18
                mdl = fitcdiscr(phTb.artic(:,1:13), phTb.cat, 'Weights', phTb.weight);
            else
                mdl = fitcdiscr(phTb.artic, phTb.cat, 'Weights', phTb.weight);
            end
            
            % Find two discriminant axes
            for i = size(phDisc, 2) : -1 : 1
                m1 = mdl.ClassNames == phDisc(i,1);
                m2 = mdl.ClassNames == phDisc(i,2);
                c(i) = mdl.Coeffs(m1,m2).Const;
                b = mdl.Coeffs(m1,m2).Linear;
                B(:,i) = b / norm(b);
            end
            
            s.phList = phList;
            s.phDisc = phDisc;
            s.mdl = mdl;
            s.B = B;
            s.c = c;
        end
        
        function PlotLDA(ax, phTb, s, ax0)
            
            % Select phonemes
            isSelect = ismember(phTb.arp, s.phList);
            phTb = phTb(isSelect, :);
            phTb = sortrows(phTb, 'arp');
            
            % Project LMV artic
            nVars = size(phTb.artic, 2);
            if nVars == 18
                Z = phTb.artic(:,1:13) * s.B;
            else
                Z = phTb.artic * s.B;
            end
            
            % Compute centroids
            [G, phCent] = findgroups(phTb.cat);
            Zcent = splitapply(@(x) median(x,1,'omitnan'), Z, G);
            
            gscatter(Z(:,1), Z(:,2), phTb.arp, MPlot.Rainbow(7), [], 8);
            text(Zcent(:,1), Zcent(:,2), MLing.ARPA2IPA(phCent), ...
                'FontSize', 24, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            h = legend();
            h.String = MLing.ARPA2IPA(h.String);
            h.Visible = 'off';
            
            ipaDisc = MLing.ARPA2IPA(s.phDisc);
            ax.XLabel.String = sprintf("LDA %s-%s", ipaDisc(1,1), ipaDisc(1,2));
            ax.YLabel.String = sprintf("LDA %s-%s", ipaDisc(2,1), ipaDisc(2,2));
            ax.Box = 'off';
            
            if exist('ax0', 'var')
                ax.XLim = ax0.XLim;
                ax.YLim = ax0.YLim;
                ax.XTick = ax0.XTick;
                ax.YTick = ax0.YTick;
            end
        end
        
        function ph = Standardize(ph)
            % Standardize phoneme labels to CMU style by grouping similar phonemes
            % 
            %   ph = NP.Phone.Standardize(ph)
            % 
            % Input
            %   ph          Phoneme labels in ARPABET.
            % Output
            %   ph          Phoneme labels with CMU style.
            % 
            
            dict = { ...
                'IX', 'IH';
                'AX', 'AH';
                'AXR', 'ER';
                'UX', 'UH';
                'EN', 'N'; % AX + N
                'EM', 'M'; % AX + M
                'EL', 'L'; % AX + L
                'ENG', 'NG'; % AX + NG
                'AX-H', 'AH';
                'BCL', 'B';
                'DCL', 'D';
                'GCL', 'G';
                'HV', 'HH';
                'KCL', 'K';
                'PCL', 'P';
                'TCL', 'T';
                };
            
            dtype = class(ph);
            
            ph = cellstr(ph);
            ph = upper(ph);
            for k = 1 : numel(ph)
                L = ph{k};
                L(L>='0' & L<='9') = []; % remove the trailing digit in phonemes
                I = find(strcmp(L, dict(:,1)), 1);
                if ~isempty(I)
                    ph{k} = dict{I,2};
                end
            end
            
            if dtype == "char"
                ph = ph{1};
            elseif dtype == "string"
                ph = string(ph);
            end
        end
        
        function c = GetPhoneColor(ph)
            % 
            
            ph = upper(ph);
            cmap = lines(7);
%             cmap = brewermap(7, 'Set1');
            
            if ismember(ph, NP.Phone.vowels)
                c = cmap(1,:);
            elseif ismember(ph, NP.Phone.fricatives)
                c = cmap(2,:);
            elseif ismember(ph, NP.Phone.plosives)
                c = cmap(3,:);
            elseif ismember(ph, NP.Phone.nasals)
                c = cmap(4,:);
            else
                c = cmap(end,:);
            end
        end
        
        function [N, C] = PhonemeHistcounts(sent, timitDir, varargin)
            % 
            
            if class(sent) == "NP.TGEvent"
                tge = sent;
            else
                tg = MLing.ReadTimitFeatures(timitDir, sent);
                tge = NP.TGEvent(tg);
            end
            
            C = string([NP.Phone.vowels NP.Phone.consonants]);
            N = zeros(numel(tge), numel(C));
            for i = 1 : numel(tge)
                ph = Cut(Cut(tge(i))).GetParentLabel;
                ph = erase(ph, digitsPattern);
                ph = categorical(ph, C);
                N(i,:) = histcounts(ph, C);
            end
        end
        
        
        
    end
end

