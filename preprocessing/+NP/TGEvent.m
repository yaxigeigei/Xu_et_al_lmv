classdef TGEvent < MSessionExplorer.Event
    % NP.TGEvent is a specialized MSessionExplorer.Event that supports mPraat for TextGrid operations
    % 
    %   obj = NP.TGEvent(tg)
    % 
    % Input
    %   tg          One or a vector of structs containing TextGrid data.
    % Output
    %   obj         One or a vector of NP.TGEvent objects.
    % 
    
    properties
        tier = {};
        tmin = NaN;
        tmax = NaN;
        
        tier1 = '';
        label = '';
        parentLabel = '';
    end
    
    methods
        % Object Construction
        function obj = TGEvent(tg)
            % See class documentation above for the use of the constructor
            
            % Default values for T
            s.dt = 0;
            s.tmin = NaN;
            s.tmax = NaN;
            s.T1 = NaN;
            s.T2 = NaN;
            
            % Allow for constructing object array
            if ~nargin || isempty(tg)
                obj.T = s;
                return
            end
            
            % Add values to object
            for i = numel(tg) : -1 : 1
                
                s.tmin = tg(i).tmin;
                s.tmax = tg(i).tmax;
                s.T1 = tg(i).tier{1}.T1;
                s.T2 = tg(i).tier{1}.T2;
                
                obj(i,1).t = s.tmin;
                obj(i,1).T = s;
                
                obj(i,1).tier = tg(i).tier;
                obj(i,1).tmin = tg(i).tmin;
                obj(i,1).tmax = tg(i).tmax;
                obj(i,1).tier1 = tg(i).tier{1}.name;
                obj(i,1).label = tg(i).tier{1}.Label;
                obj(i,1).parentLabel = strjoin(tg(i).tier{1}.Label, ' ');
            end
        end
        
        function obj = Trim(obj)
            % Remove silent segments at the begining and the end. tmin will become zero 
            % and all other time values will be adjusted accordingly. 
            % 
            %   obj = Trim(obj)
            % 
            
            for i = 1 : numel(obj)
                x = obj(i);
                
                % Trim each tier
                for t = 1 : tgGetNumberOfTiers(x)
                    % Find non-silent intervals
                    nItvl = tgGetNumberOfIntervals(x, t);
                    I = tgFindLabels(x, t, '');
                    I = [I{:}];
                    I = setdiff(1:nItvl, I);
                    
                    % Remove silence
                    T1 = x.tier{t}.T1(I);
                    T2 = x.tier{t}.T2(I);
                    L = x.tier{t}.Label(I);
                    
                    % Realign time
                    if t == 1
                        tOffset = T1(1);
                    end
                    T1 = T1 - tOffset;
                    T2 = T2 - tOffset;
                    
                    % Modify gloabl time independent variables
                    x.tier{t}.T1 = T1;
                    x.tier{t}.T2 = T2;
                    x.tier{t}.Label = L;
                    if t == 1
                        x.tmin = T1(1);
                        x.tmax = T2(end);
                        x.label = L;
                        x.parentLabel = strjoin(L, ' ');
                        
                        % Modify global time dependent variables
                        x = x + tOffset;
                        x = RefreshTfields(x);
                    end
                end
                
                obj(i) = x;
            end
        end
        
        function obj = RefreshTfields(obj)
            % Refresh field values in obj.T to reflect changes in tier 1 data
            for i = 1 : numel(obj)
                x = obj(i);
                T1 = x.tier{1}.T1;
                T2 = x.tier{1}.T2;
                dt = x.T.dt;
                x.T.T1 = T1 + dt;
                x.T.T2 = T2 + dt;
                x.T.tmin = T1(1) + dt;
                x.T.tmax = T2(end) + dt;
                obj(i) = x;
            end
        end
        
        function subObjs = Cut(obj)
            % Cut TGEvent object(s) into a vector of individual next tier elements
            % e.g. from words to phones, or phones to phone
            % 
            %   subObjs = Cut(obj)
            % 
            %   obj         A TGEvent object or a vector of objects to be cut up.
            %   subObjs     A vector of individual TGEvent objects of the next tier. 
            %               The 'parentLabel' property will be set by the corresponding 'label' element in obj.
            %               If obj is a vector, subObjs concatenates all the individuals.
            %               If obj is empty, the same empty object will be returned.
            
            % Allows array operation
            if numel(obj) > 1
                subObjs = arrayfun(@Cut, obj(:), 'Uni', false);
                subObjs = cat(1, subObjs{:});
                return
            end
            
            % Return the same obj if empty
            if isnan(obj)
                subObjs = obj;
                return
            end
            
            % Find segmentation times of the top tier
            topTier = tgGetTierName(obj, 1);
            T1 = obj.tier{1}.T1;
            T2 = obj.tier{1}.T2;
            
            switch topTier
                case 'words'
                    for i = numel(T1) : -1 : 1
                        tg = tgCut(obj, T1(i), T2(i));
                        tg = tgRemoveTier(tg, 1);
                        subObjs(i,1) = NP.TGEvent(tg);
                        subObjs(i,1).parentLabel = obj.label{i}; % overwrite the default label
                    end
                    
                case 'syll'
                    for i = numel(T1) : -1 : 1
                        tg = tgCut(obj, T1(i), T2(i));
                        tg = tgRemoveTier(tg, 1);
                        subObjs(i,1) = NP.TGEvent(tg);
                        subObjs(i,1).parentLabel = obj.label{i}; % overwrite the default label
                    end
                    
                case 'phones'
                    for i = numel(T1) : -1 : 1
                        tg = tgCut(obj, T1(i), T2(i));
                        tg.tier{1}.name = 'phone';
                        subObjs(i,1) = NP.TGEvent(tg);
                        subObjs(i,1).parentLabel = obj.label{i}; % overwrite the default label
                    end
                    
                otherwise
                    error('%s is not a supported tier name', topTier);
            end
            
            % Apply any time shift back
            subObjs = subObjs + obj.T.dt;
            
        end
        
        function obj = PickupTimeseries(obj, varName, t, v, varargin)
            % 
            % 
            %   obj = PickupTimeseries(obj, varName, t, v)
            %   obj = PickupTimeseries(obj, varName, t, v, tPad)
            %   obj = PickupTimeseries(..., 'NumResample', 0)
            % 
            % Inputs
            %   varName
            %   t
            %   v
            %   tPad
            %
            
            p = inputParser;
            p.addOptional('tPad', 0, @isnumeric);
            p.addParameter('NumResample', 0, @(x) isnumeric(x) && isscalar(x));
            p.parse(varargin{:});
            tPad = p.Results.tPad;
            N = p.Results.NumResample;
            
            if size(tPad,2) == 1
                tPad = repmat(tPad, [1 2]);
            end
            if size(tPad,1) == 1
                tPad = repmat(tPad, [numel(obj) 1]);
            end
            
            tCell = cell(size(obj));
            vCell = cell(size(obj));
            for i = 1 : numel(obj)
                t1 = obj(i).T.tmin - tPad(i,1);
                t2 = obj(i).T.tmax + tPad(i,2);
                if N
                    tCell{i} = linspace(t1, t2, N);
                    vCell{i} = interp1(t, v, t_);
                else
                    m = t >= t1 & t < t2;
                    tCell{i} = t(m);
                    vCell{i} = v(m,:);
                end
            end
            obj = obj.SetTfield(varName, tCell);
            obj = obj.SetVfield(varName, vCell);
        end
        
        function obj = AddSyllableTier(obj, lookup)
            % Add syllable tier to object(s) based on a syllabification lookup table
            % 
            %   obj = AddSyllableTier(obj, lookup)
            % 
            % Inputs
            %   obj         NP.TGEvent object(s) with top tier being 'words' or 'phones'. If 'syll' tier already 
            %               exists in an object, the 'syll' tier will be removed before adding the new one.
            %   lookup      The lookup table with the following columns:
            %               1) word: text of words
            %               2) syllabification: markups
            % Output
            %   obj         The output object with a new 'syll' tier before the 'phones' tier.
            %               A single interval with all the phones will be used as a placeholder "syllable" if:
            %               1) The word cannot be found in the lookup table.
            %               2) The syllabification for the word is 'NONE'.
            %               3) The total number of phonemes from the lookup table does not match that from the object.
            % 
            
            % Call function recursively for array input
            if numel(obj) > 1
                obj = arrayfun(@(x) x.AddSyllableTier(lookup), obj);
                return
            end
            
            % Remove existing syll tier if any
            sylTierName = "syll";
            nTiers = tgGetNumberOfTiers(obj);
            for i = 1 : nTiers
                tn = tgGetTierName(obj, i);
                if tn == sylTierName
                    obj = tgRemoveTier(obj, i);
                end
            end
            
            % Cut words into single-words
            tn = tgGetTierName(obj, 1);
            if tn == "words"
                wrds = Cut(obj);
            elseif tn == "phones"
                wrds = obj;
            else
                error("The top tier of input object must be 'words' or 'phones', but was '%s'", tn);
            end
            
            % Process each word
            sylPattern = cell(size(wrds));
            for w = 1 : numel(wrds)
                % Get word and phone labels
                wrdLabel = wrds(w).GetParentLabel;
                phLables = string(wrds(w).tier{1}.Label);
                
                k = find(lookup.word == wrdLabel, 1);
                if isempty(k)
                    warning("The word '%s' cannot be found in the lookup table.", wrdLabel);
                    sylPattern{w} = strjoin(phLables, '');
                    continue
                end
                
                % Separate syllabel markups
                sylTkn = lookup.syllabification{k};
                if sylTkn == "NONE"
                    warning("The syllabification for word '%s' is 'NONE'.", wrdLabel);
                    sylPattern{w} = strjoin(phLables, '');
                    continue
                end
                sylTkn = regexp(sylTkn, '\{([^}]*)\}', 'tokens');
                
                % Count the number of phonemes in each syllable
                syl = cell(size(sylTkn));
                ph = cell(size(sylTkn));
                for i = 1 : numel(sylTkn)
                    phTkn = regexp(sylTkn{i}, '^o: (.*?)\s*, n: (.*?) \[.*\], c: (.*?)\s*$', 'tokens');
                    phTkn = phTkn{1}{1};
                    
                    %  phTkn{2} = [phTkn{2} '(?:\d)?'];
                    
                    phTkn(phTkn=="empty") = [];
                    
                    phTkn = strjoin(phTkn, ' ');
                    phTkn = strsplit(phTkn, ' ');
                    
                    ph{i} = phTkn;
                    syl{i} = cat(2, phTkn{:});
                end
                nPhPerSy = cellfun(@numel, ph);
                
                % Group phones by syllables
                nPh = numel(phLables);
                if sum(nPhPerSy) ~= nPh
                    warning("For the word '%s', the total number of phonemes from the lookup table (%i) does not match that from the object (%i). ", ...
                        wrdLabel, sum(nPhPerSy), nPh);
                    sylPattern{w} = strjoin(phLables, '');
                    continue
                end
                phLables = mat2cell(phLables', nPhPerSy);
                sylPattern{w} = cellfun(@(x) strjoin(x, ''), phLables)';
            end
            
            % Construct mPraat syllable seperation pattern, e.g. 'NOW1-BAA2-DIY2-LAY1KS-SNEY1KS'
            sylPattern = cat(2, sylPattern{:});
            sylPattern = strjoin(sylPattern, '-');
            
            % Add syllable tier before phones
            iTier = tgI(obj, 'phones');
            obj = tgDuplicateTierMergeSegments(obj, 'phones', iTier, sylTierName, sylPattern, '-');
            
            % Update obj.label with syll tier label
            obj.label = obj.tier{1}.Label;
            
            % Update obj.T if syll is now the top tier
            % TBD
        end
        
        function obj = AddLookupValue(obj, lookup, keyCol, valCols, fieldNames)
            % Add event value(s) in obj.V from a lookup table
            % 
            %   obj = AddLookupValue(obj, lookup)
            %   obj = AddLookupValue(obj, lookup, keyCol)
            %   obj = AddLookupValue(obj, lookup, keyCol, valCols)
            %   obj = AddLookupValue(obj, lookup, keyCol, valCols, fieldNames)
            % 
            % Inputs
            %   obj             NP.TGEvent object(s).
            %   lookup          The lookup table.
            %   keyCol          The name of the column in the lookup table that is used as keys.
            %   valCols         The name(s) of the column(s) in the lookup table with the values to add.
            %   fieldNames      Renaming the fields in obj.V from valCols to fieldNames.
            % Output
            %   obj             NP.TGEvent(s) objects with added or updated V fields.
            % 
            
            if ~exist('keyCol', 'var')
                keyCol = lookup.Properties.VariableNames{1};
            end
            if ~exist('valCols', 'var')
                valCols = setdiff(lookup.Properties.VariableNames, keyCol);
            end
            if ~exist('fieldNames', 'var')
                fieldNames = valCols;
            end
            
            vTb = table;
            vTb.key = obj(:).GetParentLabel;
            for i = 1 : height(vTb)
                k = find(vTb.key(i) == lookup.(keyCol), 1);
                if isempty(k)
                    warning("'%s' cannot be found in the '%s' column of the lookup table.", vTb.key(i), keyCol);
                    continue
                end
                for j = 1 : numel(valCols)
                    val = lookup.(valCols{j})(k,:);
                    vTb.(fieldNames{j})(i,:) = val;
                end
            end
            
            if width(vTb) == 1
                warning("Did not find any matching key in the '%s' column of the lookup table.", keyCol);
                return
            end
            
            for i = 1 : numel(fieldNames)
                obj = obj.SetVfield(fieldNames{i}, vTb.(fieldNames{i}));
            end
        end
        
        % Object Modifiers
        function obj = Morph(obj, gi)
            for i = 1 : numel(obj)
                x = obj(i);
                
                % Transform global times
                tOld = x.t;
                tNew = gi(tOld);
                x.t = tNew;
                x.T = structfun(@(x) gi(x), x.T, 'Uni', false);
                
                % Relative times
                x.tmin = x.T.tmin - tNew;
                x.tmax = x.T.tmax - tNew;
                
                tr = x.tier;
                for j = 1 : numel(tr)
                    tr{j}.T1 = gi(tr{j}.T1 + tOld) - tNew;
                    tr{j}.T2 = gi(tr{j}.T2 + tOld) - tNew;
                end
                x.tier = tr;
                
                obj(i) = x;
            end
        end
        
        function objs = Resample(objs, nPtInterp)
            
            for i = 1 : numel(objs)
                obj = objs(i);
                
                % Find lick window
                tHSV = obj.T.tHSV;
                tADC = obj.T.tADC;
                tWin = tHSV([1 end]);
                ti = linspace(tWin(1), tWin(2), nPtInterp)';
                
                % Interpolate HSV data
                if numel(tHSV) > 1
                    obj.length = interp1(tHSV, obj.length, ti, 'linear');
                    obj.velocity = interp1(tHSV, obj.velocity, ti, 'linear');
                    obj.angle = interp1(tHSV, obj.angle, ti, 'linear');
                else
                    obj.length = NaN(size(ti));
                    obj.velocity = NaN(size(ti));
                    obj.angle = NaN(size(ti));
                end
                
                % Interpolate ADC data
                if numel(tADC) > 1
                    obj.forceV = interp1(tADC, obj.forceV, ti, 'linear', 0);
                    obj.forceH = interp1(tADC, obj.forceH, ti, 'linear', 0);
                else
                    obj.forceV = NaN(size(ti));
                    obj.forceH = NaN(size(ti));
                end
                obj.force = sqrt(obj.forceH.^2 + obj.forceV.^2);
                
                % Normalize time
                obj.T.ti = ti;
                if ti(1) == ti(end)
                    gi = griddedInterpolant(ti(1)+[-1 1]', [-1 1]');
                else
                    gi = griddedInterpolant(ti([1 end]), [-1 1]');
                end
                obj = obj.Morph(gi);
                
                objs(i) = obj;
            end
        end
        
        % Getters
        function ss = GetParentLabel(objs)
            % Get parentLabel
            ss = arrayfun(@(x) string(x.parentLabel), objs);
        end
        
        function s = GetAllParentLabel(objs, delimiter)
            % Return a string that concatenates all parentLabels in the input object(s)
            % 
            %   s = GetAllParentLabel(objs)
            %   s = GetAllParentLabel(objs, delimiter)
            % 
            % Input
            %   delimiter       The delimiter to use when concatenating labels. The default is " ".
            % Output
            %   s               A string of concatenated parent labels.
            % 
            if nargin < 2
                delimiter = " ";
            end
            s = objs.GetParentLabel();
            if all(s == "")
                s = "";
            else
                s = join(s, delimiter);
                s = strtrim(s);
            end
        end
        
        function dur = GetDuration(objs)
            % Return the total duration of the object(s)
            % 
            %   dur = GetDuration(objs)
            % 
            dur = arrayfun(@(x) x.T.tmax - x.T.tmin, objs);
        end
        
        function tf = IsVoiced(objs)
            % Check if phone object(s) are voiced or not
            %
            %   tf = IsVoiced(objs)
            % 
            tf = false(size(objs));
            for i = 1 : numel(objs)
                obj = objs(i);
                tierName = tgGetTierName(obj, 1);
                if tierName ~= "phone"
                    error("Object must be at the 'phone' tier, but was '%s'.", tierName);
                end
                L = upper(obj.parentLabel);
                L(L>='0' & L<='9') = []; % remove the trailing digit in phonemes
                tf(i) = any(strcmp(L, NP.Phone.voiced));
            end
        end
        
        function tf = IsVowel(objs)
            % Check if phone object(s) are voiced or not
            %
            %   tf = IsVowel(objs)
            % 
            tf = false(size(objs));
            for i = 1 : numel(objs)
                obj = objs(i);
                tierName = tgGetTierName(obj, 1);
                if tierName ~= "phone"
                    error("Object must be at the 'phone' tier, but was '%s'.", tierName);
                end
                L = upper(obj.parentLabel);
                L(L>='0' & L<='9') = []; % remove the trailing digit in phonemes
                tf(i) = any(strcmp(L, NP.Phone.vowels));
            end
        end
        
        % Utilities
        function mask = MaskTimestamps(objs, t, tPad)
            % Get a mask for a vector of timestamps where samples spanned by utterances are set to true
            % 
            %   mask = MaskTimestamps(objs, t)
            %   mask = MaskTimestamps(objs, t, tPad)
            % 
            % Inputs
            %   objs        A TGEvent object or a vector of objects.
            %   t           A numeric vector of timestamps.
            %   tPad        The amount of time to pad at the onset and offset of each object. Positive 
            %               padding extends the boundaries, negative padding shortens them. Default tPad
            %               is zero and no padding is performed.
            %               1) If tPad is a scalar, the same amount is padded to both onsets and offsets
            %                  for all objects.
            %               2) If tPad is a 1-by-2 vector, the first element is used to pad onsets and the 
            %                  second element is used to pad offsets for all objects.
            %               3) If tPad is an n-by-1 or n-by-2 vector where n is the number of objects, 
            %                  each object will use the value(s) in the corresponding row of tPad.
            % Output
            %   mask        A logical vector with the same size as t.
            % 
            
            if nargin < 3
                tPad = 0;
            end
            if size(tPad,2) == 1
                tPad = repmat(tPad, [1 2]);
            end
            if size(tPad,1) == 1
                tPad = repmat(tPad, [numel(objs) 1]);
            end
            
            mask = false(size(t));
            for i = 1 : numel(objs)
                t1 = objs(i).T.tmin - tPad(i,1);
                t2 = objs(i).T.tmax + tPad(i,2);
                
                ind = t >= t1 & t < t2;
                mask(ind) = true;
            end
        end
        
        % Visualization
        function Plot(objs, varargin)
            % Plot text labels as a function of time
            %
            %   Plot(objs)
            %   Plot(objs, yy)
            %   Plot(..., 'Style', 'line')
            %   Plot(..., 'Color', [0 0 0])
            %   Plot(..., 'FontSize', 6)
            %   Plot(..., <lineArgs>)
            % 
            % Inputs
            %   objs            NP.TGEvent object(s)
            %   yy              1) A two-element vector for the Y range of the delimiter bars. Text labels 
            %                      will be positioned in the middle of this range.
            %                   2) A scalar for the Y position of text labels. To make delimiters visible, 
            %                      one can specify the plot's marker type, e.g. 'Marker', '*'.
            %   'Style'         1) 'line' (default): 
            %                   2) 'patch': 
            %                   3) Empty string '': no boundary will be plotted.
            %   'Color'         
            %   'FontSize'      Font size of text labels. Use zero to hide text labels.
            %   <lineArgs>      Other argument-value pairs for the MATLAB 'plot' function.
            % 
            
            p = inputParser();
            p.KeepUnmatched = true;
            p.addOptional('yy', [], @isnumeric);
            p.addParameter('Style', 'line', @(x) ischar(x) || isstring(x));
            p.addParameter('Color', [], @(x) isnumeric(x) || ischar(x));
            p.addParameter('FontSize', 10, @(x) true);
            p.addParameter('Parent', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
            yy = p.Results.yy;
            style = p.Results.Style;
            cc = p.Results.Color;
            fontSz = p.Results.FontSize;
            ax = p.Results.Parent;
            lineArgs = p.Unmatched;
            
            if isempty(ax)
                ax = gca;
            end
            
            if isempty(yy)
                h = gca;
                yy = h.YLim;
            end
            yy = yy(:);
            
            for i = 1 : numel(objs)
                % Get text label
                s = objs(i).parentLabel;
                if objs(i).tier1 == "phone"
                    s(s>='0' & s<='9') = []; % remove the trailing digit (stress)
                end
                
                % Get interval window
                win = [objs(i).T.tmin, objs(i).T.tmax];
                
                % Plot
                hold(ax, 'on');
                if style == "line"
                    if isempty(cc)
                        lineArgs.Color = [0 0 0];
                    else
                        lineArgs.Color = cc;
                    end
                    plot(ax, [win; win], [yy yy], lineArgs);
                elseif style == "patch"
                    if isempty(cc)
                        bkColor = lines(numel(objs));
                        if objs(i).tier1 == "phone"
                            bkColor(i,:) = NP.Phone.GetPhoneColor(s);
                        end
                    else
                        bkColor = cc;
                    end
                    MPlot.Blocks(win, yy', bkColor(i,:), 'Parent', ax);
                end
                
                % Label
                if fontSz > 0
                    if objs(i).tier1 == "phone"
                        s = MLing.ARPA2IPA(s);
                    end
                    if isempty(cc)
                        txtColor = [0 0 0];
                    else
                        txtColor = cc;
                    end
                    text(ax, mean(win), mean(yy), s, 'Color', txtColor, 'FontSize', fontSz, ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                end
            end
        end
        
        function PlotTiers(objs, varargin)
            % Plot a hierarchy of speech labels for a given trial
            % 
            %   PlotTiers(objs)
            %   PlotTiers(objs, y0)
            %   PlotTiers(objs, y0, dy)
            %   PlotTiers(objs, y0, dy, h)
            %   PlotTiers(objs, ..., 'TierInd', [])
            %   PlotTiers(objs, ..., 'Color', {'k', 'b', 'm', 'r'})
            % 
            % Inputs
            %   y0              The y-coordinate of the base. Default is 0.5.
            %   dy              The y distance between adjacent tiers. Default is 1.
            %   h               The height of separator bars. Use 0 for no separators. Default is the same as dy.
            %   'TierInd'       The indices of tiers to plot. Default is empty [], plotting all tiers.
            %   'Color'         A cell array of color vectors or letters.
            % 
            % See also NP.TGEvent.Plot
            
            p = inputParser();
            p.KeepUnmatched = true;
            p.addOptional('y0', 0.5, @(x) isnumeric(x) && isscalar(x));
            p.addOptional('dy', 1, @(x) isnumeric(x) && isscalar(x));
            p.addOptional('h', [], @(x) isnumeric(x) && isscalar(x));
            p.addParameter('TierInd', [], @(x) isnumeric(x));
            p.addParameter('Color', {'k', 'b', 'm', 'r'}, @(x) isnumeric(x) || iscellstr(x));
            p.parse(varargin{:});
            y0 = p.Results.y0;
            dy = p.Results.dy;
            h = p.Results.h;
            indT = p.Results.TierInd;
            cc = p.Results.Color;
            
            % Determine tiers to plot
            nTier = tgGetNumberOfTiers(objs(1)) + 1;
            if isempty(indT)
                indT = 1 : nTier;
            end
            indT = MMath.Bound(indT, [1 nTier]);
            indT = unique(indT);
            
            % Determine the relative separator range
            if isempty(h)
                h = dy;
            end
            if h ~= 0
                hRange = [0 h];
            else
                hRange = dy/2;
            end
            
            % Plot tiers
            k = 0;
            for i = 1 : nTier
                if ismember(i, indT)
                    k = k + 1;
                    objs.Plot(hRange+y0+(k-1)*dy, 'Color', cc{k}, p.Unmatched);
                end
                if i < nTier
                    objs = Cut(objs);
                end
            end
            if isfield(p.Unmatched, 'Parent')
                ax = p.Unmatched.Parent;
            else
                ax = gca;
            end
            ax.YDir = 'reverse';
        end
        
    end
end

