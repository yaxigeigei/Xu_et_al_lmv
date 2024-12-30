classdef Morpher < handle
    % Morpher objects morph time and timeseries data by interpolating key times
    %
    %   obj = Morpher(tFrom, tTo)
    %   obj = Morpher(tFrom, tTo, tPad)
    %
    % Inputs
    %   tFrom       Original key times to be morphed. This can be a numeric vector for 
    %               constructing a single morpher object, or an n-element cell array of such 
    %               vectors for constructing a vector of morpher objects in batch.
    %   tTo         Template key times to morph to. The dimensions must match with tFrom.
    %   tPad        The amount of time padded to the beginning and the end of tFrom and tTo.
    %               Must be an n-element vector of positive real number. If tFrom is a vector, 
    %               tPad would be a scalar.
    % Output
    %   obj         A morpher object or a vector of objects.
    % 
    
    properties
        tFrom;      % Original key times to be morphed
        tTo;        % Template key times to morph to
        giObj;      % griddedInterpolant object, or []
    end
    
    methods
        function this = Morpher(tFrom, tTo, tPad)
            % Morpher construct an instance of this class
            
            % Handle user inputs
            if nargin == 0
                % Allow for constructing array
                return;
            end
            if ~iscell(tFrom)
                tFrom = {tFrom};
            end
            if ~iscell(tTo)
                tTo = {tTo};
            end
            if ~exist('tPad', 'var')
                tPad = 1;
            end
            if numel(tPad) == 1
                tPad = repmat(tPad, size(tFrom));
            end
            
            for k = numel(tFrom) : -1 : 1
                % Cache values
                fm = tFrom{k}(:);
                to = tTo{k}(:);
                pd = tPad(k);
                
                assert(numel(fm)==numel(to), ...
                    "The number of elements in tForm (%i) and tTo (%i) for the #%i pair are not the same.\n", ...
                    numel(fm), numel(to), k);
                
                if ~isempty(fm) && ~isempty(to)
                    % Add padding
                    fm = cat(1, fm(1)-pd, fm, fm(end)+pd);
                    to = cat(1, to(1)-pd, to, to(end)+pd);
                    
                    % Make interpolant
                    gi = griddedInterpolant(fm, to, 'linear', 'linear');
                else
                    gi = [];
                end
                
                % Save values
                this(k,1).tFrom = fm;
                this(k,1).tTo = to;
                this(k,1).giObj = gi;
            end
        end
        
        function val = IsMorph(this)
            val = ~arrayfun(@(x) isempty(x.giObj), this);
        end
        
        function [t, v] = Morph(this, t, v)
            % Apply mapping to event times or time series
            
            assert(numel(this) == 1, 'This method can only be called by one object at a time');
            
            % Return inputs if no morphing is needed
            if ~this.IsMorph || isempty(t)
                return
            end
            
            % Map sample times
            if isnumeric(t)
                tOld = t;
                t = this.giObj(t);
            else
                t = t.Morph(this.giObj); % require a Morph interface that accepts griddedInterpolant object
                return
            end
            
            % Resmaple at original timepoints after morphing the timeseries
            if nargin > 2
                dtype = class(v); % remember original data type
                gi = griddedInterpolant(t, double(v), 'linear', 'nearest');
                v = gi(tOld);
                v = cast(v, dtype); % restore original data type
                t = tOld;
            end
        end
        
        function etTb = MorphEventTimes(this, etTb)
            % Apply mapping to an event time table
            
            assert(numel(this) == height(etTb), ...
                'There are %d Morpher objects but %d trials in the table', ...
                numel(this), height(etTb));
            
            for i = 1 : width(etTb)
                etCol = etTb.(i);
                for k = 1 : numel(this)
                    if ~iscell(etCol)
                        etCol(k) = this(k).Morph(etCol(k));
                    else
                        etCol{k} = this(k).Morph(etCol{k});
                    end
                end
                etTb.(i) = etCol;
            end
        end
        
        function tsTb = MorphTimeSeries(this, tsTb)
            % Apply mapping to an time series table
            
            assert(numel(this) == height(tsTb), ...
                'There are %d Morpher objects but %d trials in the table', ...
                numel(this), height(tsTb));
            
            for k = 1 : numel(this)
                tsTb.time{k} = this(k).Morph(tsTb.time{k});
%                 for i = 2 : width(tsTb)
%                     [~, tsTb.(i){k}] = this(k).Morph(tsTb.time{k}, tsTb.(i){k});
%                 end
            end
        end
        
        function MorphSE(this, se)
            % Morph times in each eventTimes or timeSeries table
            for i = 1 : numel(se.tableNames)
                if se.isEventValuesTable(i)
                    continue;
                end
                tb = se.GetTable(se.tableNames{i});
                if se.isEventTimesTable(i)
                    tb = this.MorphEventTimes(tb);
                else
                    tb = this.MorphTimeSeries(tb);
                end
                se.SetTable(se.tableNames{i}, tb);
            end
        end
    end
    
end
