classdef CodingExplorer < MSessionExplorer
    % Manage and operate with data for coding analysis
    % 
    %   ce = CodingExplorer()
    %   ce = CodingExplorer(se)
    % 
    
    properties
        
    end
    
    properties(Dependent)
        numFeat             % the number of variables in the 'feat' timeseries table
        numResp             % the number of variables in the 'resp' timeseries table
        featNames           % the names of variables in the 'feat' timeseries table
        respNames           % the names of variables in the 'resp' timeseries table
        clusTb              % ce.userData.ksMeta.clusTb
    end
    methods
        function val = get.numFeat(this)
            if ismember('feat', this.tableNames)
                val = width(this.GetTable('feat')) - 1;
            else
                val = 0;
            end
        end
        function val = get.numResp(this)
            if ismember('resp', this.tableNames)
                val = width(this.GetTable('resp')) - 1;
            else
                val = 0;
            end
        end
        function val = get.featNames(this)
            if this.numFeat
                tb = this.GetTable('feat');
                val = tb.Properties.VariableNames(2:end);
            else
                val = {};
            end
        end
        function val = get.respNames(this)
            if this.numResp
                tb = this.GetTable('resp');
                val = tb.Properties.VariableNames(2:end);
            else
                val = {};
            end
        end
        function val = get.clusTb(this)
            if isfield(this.userData, 'ksMeta')
                val = NP.Unit.GetClusTb(this.userData);
            else
                val = [];
            end
        end
        function set.clusTb(this, val)
            this.userData.ksMeta.clusTb = val;
        end
    end
    
    methods
        function this = CodingExplorer(se)
            %PETHEXPLORER Construct an instance of this class
            
            % Parent constructor
            this = this@MSessionExplorer();
            
            % Allows construction of array
            if nargin == 0
                return
            end
            
            % Copy content from se to ce
            types = se.supportedTableTypes;
            typeMask = [se.isEventTimesTable; se.isEventValuesTable; se.isTimesSeriesTable];
            for i = 1 : numel(se.tableNames)
                name = se.tableNames{i};
                type = types{typeMask(:,i)};
                this.SetTable(name, se.GetTable(name), type, se.GetReferenceTime(name));
            end
            this.userData = se.userData;
            this.isVerbose = se.isVerbose;
        end
        
        function ceCat = CatUnits(varargin)
            % Combine multiple NP.CodingExplorer objects into one by horizontally concatenating units
            % 
            %   ceCat = CatUnits(ce1, ce2, ce3, ...)
            % 
            % Inputs
            %   ce1, ce2, ce3, ...      Arbitrary number of NP.PethExplorer objects or object arrays
            % Output
            %   ceCat                   The combined NP.PethExplorer object
            % 
            
            % Collect objects
            for i = 1 : numel(varargin)
                varargin{i} = varargin{i}(:);
            end
            ceArray = cat(1, varargin{:});
            
            % Initialize concatenated pe
            ceCat = ceArray(1).Duplicate;
            
            % Concatenate and set each table
            for i = 1 : numel(ceCat.tableNames)
                % Will not concatenate eventValues table
                tn = ceCat.tableNames{i};
                m = strcmp(tn, ceCat.tableNames);
                if ~any(m)
                    fprintf("Skip '%s' table since it is not present in the first ce object.\n", tn);
                    continue
                end
                
                % Get all tables
                tbs = arrayfun(@(x) x.GetTable(tn), ceArray, 'Uni', false);
                if ceCat.isTimesSeriesTable(i)
                    tbs(2:end) = cellfun(@(x) x(:,2:end), tbs(2:end), 'Uni', false);
                end
                
                % Concatenate if all column variables are unique
                varNames = cellfun(@(x) x.Properties.VariableNames, tbs, 'Uni', false);
                varNames = cat(2, varNames{:});
                if numel(unique(varNames)) < numel(varNames)
                    fprintf("Skip '%s' table due to the presence of conflicting variable names.\n", tn);
                    continue
                end
                tbCat = cat(2, tbs{:});
                
                % Set table to cat pe
                ceCat.SetTable(tn, tbCat, ceCat.tot.tableType{i}, ceCat.tot.referenceTime{i});
            end
            
            % Concatenate clusTb
            clusTbs = arrayfun(@(x) NP.Unit.GetClusTb(x), ceArray, 'Uni', false);
            ceCat.userData.ksMeta.clusTb = cat(1, clusTbs{:});
        end
        
        function [T, V, varNames] = GetArray(this, tbIn, varargin)
            % Get PETH data in numeric array(s)
            % 
            %   [T, V, varNames] = GetArray(tbIn)
            %   [T, V, varNames] = GetArray(tbIn, rowInd)
            %   [T, V, varNames] = GetArray(tbIn, rowInd, colInd)
            %   [T, V, varNames] = GetArray(..., 'TimeShifts', 0)
            %   [T, V, varNames] = GetArray(..., 'TimeEach', false)
            %   [T, V, varNames] = GetArray(..., 'DimCat', [])
            %   [T, V, varNames] = GetArray(..., 'DimAverage', [])
            %   [T, V, varNames] = GetArray(..., 'DimCombine', [])
            %   [T, V, varNames] = GetArray(..., 'Normalization', 'none')
            % 
            % Inputs
            %   tbIn                A table in timeSeries format or the name of a timeSeries table in the object.
            %   rowInd              Integer or logical indices of rows to operate on and return. The default 
            %                       value is empty [] indicating all rows. 
            %   colInd              Integer or logical indices of columns (including the time column) to operate 
            %                       on and return. It can also be a cell array of column names of the input table. 
            %                       The default is empty [] indicating all columns.
            %   'TimeShifts'        A vector of time shifts (in sec) applied to the data. Each shift value 
            %                       results in one set of output variables. For example, the three time shifts 
            %                       [-1 0 1] applied to 10 variables will output 3x10=30 variables.
            %   'TimeEach'          Whether or not to duplicate timestamps for each timeseries and output T 
            %                       with the same size as V.
            %   'DimCat'            The dimension along which to concatenate epoch arrays. Default is 1, along 
            %                       the time dimension. Use zero for no concatenation.
            %   The following parameters apply only when 'DimCat' is non-zero
            %   'DimAverage'        The dimensions to collapse by averaging.
            %   'DimCombine'        The dimensions to reshape into a vector and put in the first dimension.
            %   'Normalization'     See the 'normType' parameter in the MMath.Normalize function.
            % Outputs
            %   T                   time-by-1 or time-by-var numeric array of timestamps, or an epoch-by-1 cell 
            %                       array of the numeric arrays if not denested.
            %   V                   time-by-var numeric array of variables, or an epoch-by-1 cell array of the 
            %                       numeric arrays if not denested.
            %   varNames            Variable names that match the columns of V.
            % 
            % See also circshift, MMath.Normalize, MMath.SqueezeDims, MMath.CombineDims
            
            p = inputParser;
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x) || iscellstr(x) || isstring(x) || ischar(x));
            p.addParameter('TimeShifts', 0, @isnumeric);
            p.addParameter('TimeEach', false, @islogical);
            p.addParameter('DimCat', 1, @isscalar);
            p.addParameter('Normalization', 'none', @(x) ischar(x) || isstring(x));
            p.addParameter('DimAverage', [], @isnumeric);
            p.addParameter('DimCombine', [], @isnumeric);
            p.parse(varargin{:});
            rowInd = p.Results.rowInd;
            colInd = p.Results.colInd;
            normType = p.Results.Normalization;
            isTimeEach = p.Results.TimeEach;
            dimCat = p.Results.DimCat;
            dimAvg = p.Results.DimAverage;
            dimCombine = p.Results.DimCombine;
            tShifts = p.Results.TimeShifts(:)';
            
            % Get table
            if ~istable(tbIn)
                assert(this.IValidateTableName(tbIn, true) == 3, '''%s'' is not a timeSeries table', tbIn);
                tbIn = this.GetTable(tbIn);
            end
            
            % Validate and standardize row and column indices
            [rowInd, colInd] = this.IValidateTableIndexing(tbIn, rowInd, colInd);
            if ~ismember(1, colInd)
                colInd = [1 colInd]; % make sure the time column is always included
            end
            
            % Select rows and columns
            tbIn = tbIn(rowInd, colInd);
            
            T = tbIn.time;
            V = cell(size(T));
            for i = 1 : height(tbIn)
                % Concatenate variables
                tbRow = tbIn{i,2:end};
                tbRow = cellfun(@double, tbRow, 'Uni', false);
                v = cell2mat(tbRow);
                
                % Expand variables with time shifts
                t = T{i};
                indShifts = round(tShifts / diff(t(1:2)));
                
                vShift = cell(size(indShifts));
                for n = 1 : numel(indShifts)
                    ns = indShifts(n);
                    vs = circshift(v, ns, 1);
                    if ns > 0
                        vs(1:ns,:) = 0;
                    elseif ns < 0
                        vs(end+ns+1:end,:) = 0;
                    end
                    vShift{n} = vs;
                end
                v = cat(2, vShift{:});
                
                % Expand timestamps
                if isTimeEach
                    t = repmat(t, [1 size(v,2)]);
                end
                
                T{i} = t;
                V{i} = v;
            end
            
            % Feature name expansion
            varNames = tbIn.Properties.VariableNames;
            for i = 1 : width(tbIn)
                nSubCol = size(tbIn.(i){1}, 2);
                if nSubCol == 1
                    varNames{i} = string(varNames{i});
                else
                    varNames{i} = varNames{i} + "_" + string(1:nSubCol);
                end
            end
            varNames = cat(2 , varNames{2:end}); % exclude time
            
            % Denest arrays
            if ~dimCat
                return % all following operations require numeric arrays
            end
            T = cat(dimCat, T{:});
            V = cat(dimCat, V{:});
            
            % Averaging and combining dimensions
            if dimAvg
                T = mean(T, dimAvg);
                V = mean(V, dimAvg, 'omitnan');
%                 T = MMath.SqueezeDims(T, dimAvg);
%                 V = MMath.SqueezeDims(V, dimAvg);
            end
            if dimCombine
                T = MMath.CombineDims(T, dimCombine);
                V = MMath.CombineDims(V, dimCombine);
            end
            
            % Normalization
            V = MMath.Normalize(V, normType, NP.Param.normAddMax);
        end
        
        function RemoveUnits(this, uInd)
            % Remove specified units across all data structures
            % 
            %   RemoveUnits(uInd)
            % 
            % Input
            %   uInd        Integer or logical indices of the units to be removed.
            % 
            
            % Get unit IDs for later use
            u2rm = "u" + this.clusTb.clusId(uInd);
            
            % Remove rows in clusTb
            this.clusTb(uInd,:) = [];
            
            % Remove these unit variables across tables
            for i = 1 : numel(this.tableNames)
                tn = this.tableNames{i};
                tb = this.GetTable(tn);
                m = ismember(tb.Properties.VariableNames, u2rm);
                tb(:,m) = [];
                this.SetTable(tn, tb);
            end
        end
        
    end
end
