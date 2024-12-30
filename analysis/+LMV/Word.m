classdef Word
    
    methods(Static)
        function AddWordEvents(se, lookupPath)
            % Add word events and word level feature times to 'taskTime' table
            % 
            %   AddWordEvents(se)
            %   AddWordEvents(se, lookupPath)
            % 
            
            if nargin < 2
                lookupPath = "lmv_syllabification_lookup.csv";
            end
            lookup = readtable(lookupPath);
            
            tt = se.GetTable('taskTime');
            tgeSrc = ["stim", "prod"];
            
            % Put all event objects in cell array
            tgeCell = cell(height(tt), numel(tgeSrc));
            for i = 1 : numel(tgeSrc)
                nm = tgeSrc(i);
                if iscell(tt.(nm))
                    tgeCell(:,i) = tt.(nm);
                else
                    tgeCell(:,i) = num2cell(tt.(nm));
                end
            end
            
            % Add word and syllable objects to taskTime table
            for i = 1 : height(tt)
                % Combine sentence events from all sources
                sen = cat(1, tgeCell{i,:});
                sen = sort(sen);
                sen(isnan(sen)) = [];
                if isempty(sen)
                    continue
                end
                
                % Cut sentences to words
                wrd = Cut(sen);
                wrd = wrd.AddLookupValue(lookup);
                
                % Cut words to syllabels
                wrd = wrd.AddSyllableTier(lookup);
                syll = Cut(wrd);
                
                % Create variables in taskTime table
                tt.word{i} = wrd;
                tt.wordOn{i} = double(wrd);
                tt.wordOff{i} = tt.wordOn{i} + wrd.GetDuration;
                tt.syll{i} = syll;
                tt.syllOn{i} = double(syll);
                tt.syllOff{i} = tt.syllOn{i} + syll.GetDuration;
            end
            se.SetTable('taskTime', tt);
        end
        
        function AddWordTimeseriesTable(se)
            % Add 'word' timeseries table containing word level features to se
            % 
            %   AddWordTimeseriesTable(se)
            % 
            tt = se.GetTable('taskTime');
            tsTb = table;
            tsTb.time = cellfun(@double, tt.word, 'Uni', false);
            tsTb.nSyll = cellfun(@(x) x.GetVfield('num_syllables'), tt.word, 'Uni', false);
            tsTb.blick = cellfun(@(x) x.GetVfield('blick_score'), tt.word, 'Uni', false);
            tsTb.wordFreq = cellfun(@(x) x.GetVfield('word_freq'), tt.word, 'Uni', false);
            se.SetTable('word', tsTb, 'timeSeries', se.GetReferenceTime);
        end
        
        function AddStratifiedEventTimeTable(se, nStrat)
            % Add word level feature times to the 'wordTime' table in se
            % These include "nSyll", "blick", "wordFreq"
            % 
            %   AddStratifiedEventTimeTable(se, nStrat)
            % 
            
            fieldNames = ["num_syllables", "blick_score", "word_freq"];
            colNames = ["nSyll", "blick", "wordFreq"];
            if nargin < 2
                nStrat = [2, 3, 3];
            end
            if isscalar(nStrat)
                nStrat = repelem(nStrat, numel(fieldNames));
            end
            
            tt = se.GetTable('taskTime');
            nPerEpoch = cellfun(@numel, tt.word);
            t = cellfun(@double, tt.word, 'Uni', false);
            
            if ismember('wordTime', se.tableNames)
                wtTb = se.GetTable('wordTime');
            else
                wtTb = table;
            end
            
            for n = 1 : numel(colNames)
                % Get feature values
                fn = fieldNames(n);
                val = cellfun(@(x) x.GetVfield(fn), tt.word, 'Uni', false);
                val = cat(1, val{:});
                
                % Stratify feature values
                if fn == "num_syllables"
                    bins = val;
                    bins(bins > 1) = 2; % group all multi-syllabic words
                else
                    edges = prctile(val, 0 : 100/nStrat(n) : 100)';
                    [~, ~, bins] = histcounts(val, edges);
                end
                bins = mat2cell(bins, nPerEpoch);
                
                % Pick events from each bin
                for s = 1 : nStrat(n)
                    wtTb.(colNames(n)+s) = cellfun(@(a,b) a(b==s), t, bins, 'Uni', false);
                end
            end
            
            se.SetTable('wordTime', wtTb, 'eventTimes', se.GetReferenceTime);
        end
        
        function AddUniqueLabelTime(se, tgeName, nMin)
            % Add onset times of unique words to the 'wordTime' table in se
            % 
            %   AddUniqueObjTime(se, tgeName, nMin)
            % 
            
            if ismember('wordTime', se.tableNames)
                wtTb = se.GetTable('wordTime');
            else
                wtTb = table;
            end
            
            tge = se.GetTable('taskTime').(tgeName);
            lb = cellfun(@(x) erasePunctuation(erase(x.GetParentLabel, digitsPattern)), tge, 'Uni', false);
            allLb = cat(1, lb{:});
            [N, uniLb] = histcounts(categorical(allLb));
            uniLb = string(uniLb);
            [N, I] = sort(N, 'descend');
            I = I(N >= nMin);
            N = N(N >= nMin);
            uniLb = uniLb(I);
            
            tb = table;
            tb.label = uniLb';
            tb.N = N';
            disp(tb);
            
            for j = 1 : numel(uniLb)
                nm = uniLb(j);
                for i = numel(tge) : -1 : 1
                    isLb = lb{i} == nm;
                    wtTb.(nm){i} = double(tge{i}(isLb));
                end
            end
            
            se.SetTable('wordTime', wtTb, 'eventTimes', se.GetReferenceTime);
        end
        
        function AddUniqueWordTimes(se, nMin)
            % Add onset times of unique words to the 'wordTime' table in se
            % 
            %   AddUniqueWordTimes(se, nMin)
            % 
            if ~exist('nMin', 'var')
                nMin = 10;
            end
            LMV.Word.AddUniqueLabelTime(se, 'word', nMin);
        end
        
        function AddUniqueSyllTimes(se, nMin)
            % Add onset times of unique syllables to the 'wordTime' table in se
            % 
            %   AddUniqueSyllTimes(se, nMin)
            % 
            if ~exist('nMin', 'var')
                nMin = 20;
            end
            LMV.Word.AddUniqueLabelTime(se, 'syll', nMin);
        end
        
    end
    
end