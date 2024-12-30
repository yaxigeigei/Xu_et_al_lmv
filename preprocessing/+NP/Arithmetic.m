classdef Arithmetic < NP.TaskBaseClass
    
    methods(Static)
        function [uId, uIdUni] = GetSelectedClusId(recId)
            % Return handpicked cluster IDs for a given recording
            % 
            %   [uId, uIdUni] = GetSelectedClusId(recId)
            % 
            switch recId
                % case 'NP77_B2', uId = [770210453, 770210440, 770210436, 770210430, 770210431, 770210418, 770210539, 770210370, 770210530, 770210254, 770210253, 770210233, 770210520, 770210213, 770210294, 770210315, 770210439]; %auto
                case 'NP77_B2', uId = [770200115, 770200606, 770200294, 770200233, 770200599, 770200213, 770200594, 770200119, 770200047, 770210189, 770210051];
                case 'NP89_B1', uId = [890100241, 890100157, 890100360, 890100063, 890100276 890100379, 890100204, 890100177, 890100178, 890100362, 890100033, 890110089]; 
                % case 'NP89_B1', uId = [890100367, 890100157, 890100252, 890100362, 890100101, 890100089, 890100090, 890110100]; %semsr-selective units, for reference
                otherwise, uId = [];
            end
            uIdUni = NP.Unit.GetBaseClusId(recId) + uId;
        end
        
        function numberMap = createNumberMap()
                % Initialize the map
                numberMap = containers.Map();
            
                % TODO: MFA is going to split these at the hyphens into e.g.
                % "twenty" "three" – we can just work with the ParentLabel?
    
                % Manually add each number and its text representation
                numbers = {'zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', ...
                           'ten', 'eleven', 'twelve', 'thirteen', 'fourteen', 'fifteen', 'sixteen', 'seventeen', 'eighteen', 'nineteen', ...
                           'twenty', 'twenty-one', 'twenty-two', 'twenty-three', 'twenty-four', 'twenty-five', 'twenty-six', 'twenty-seven', 'twenty-eight', 'twenty-nine', ...
                           'thirty', 'thirty-one', 'thirty-two', 'thirty-three', 'thirty-four', 'thirty-five', 'thirty-six', 'thirty-seven', 'thirty-eight', 'thirty-nine', ...
                           'forty', 'forty-one', 'forty-two', 'forty-three', 'forty-four', 'forty-five', 'forty-six', 'forty-seven', 'forty-eight', 'forty-nine', ...
                           'fifty', 'fifty-one', 'fifty-two', 'fifty-three', 'fifty-four', 'fifty-five', 'fifty-six', 'fifty-seven', 'fifty-eight', 'fifty-nine', ...
                           'sixty', 'sixty-one', 'sixty-two', 'sixty-three', 'sixty-four', 'sixty-five', 'sixty-six', 'sixty-seven', 'sixty-eight', 'sixty-nine', ...
                           'seventy', 'seventy-one', 'seventy-two', 'seventy-three', 'seventy-four', 'seventy-five', 'seventy-six', 'seventy-seven', 'seventy-eight', 'seventy-nine', ...
                           'eighty', 'eighty-one'};
            
                for i = 0:81
                    numberMap(numbers{i+1}) = i;
                end
            
                % Return the map
        end

        function tv = AddTaskValueTable(se)
            % Make and add taskValue table to se
            
            % Continue with existing table or create new
            if ismember("taskValue", se.tableNames)
                tv = se.GetTable("taskValue");
            else
                tv = table;
            end
            tt = se.GetTable("taskTime");
            % TODO: Need to be careful that we're actually doing this additively in the
            % situation where tv already exists and has data from another
            % task.

            
            % Check required columns in taskTime table
            requiredCols = {'stim', 'prod'};
            assert(all(ismember(requiredCols, tt.Properties.VariableNames)), ...
                "Arithmetic task requires the '%s' and '%s' columns in the taskTime table.", requiredCols{:});
            
            numTrials = 0;
            
            % Mapping of text to numbers
            numberMap = NP.Arithmetic.createNumberMap();
            
            for i = 1:height(tt)
                % TODO: this assumes only one prod event per stim. Need to
                % account for when this is not true and take the last prod
                % event.
                estim = tt.stim(i);
                eprod = tt.prod{i};
                
                %Rare cases, NP89 label error:
                if(length(eprod) > 1)
                    eprod = eprod(end);
                end
                
                taskName = lower(estim.GetVfield('task'));
                if ~strcmp(taskName, "arithmetic")
                    continue
                end
                
                % Get stim and prod texts
                dtype = class(estim);
                if dtype == "NP.TGEvent"
                    stimText = estim.GetParentLabel;
                    prodText = eprod.GetParentLabel;
                    if(length(prodText)>1)
                        %We annotate clarifications, but the convention is
                        %that the last word transcribed is the final
                        %answer.
                        prodText = prodText(end);
                    end
                elseif dtype == "MSessionExplorer.Event"
                    stimText = eStim(1).GetVfield('type');
                    prodText = string(NaN);
                else
                    error("Unrecognized event data type '%s'", dtype);
                end
                
                numTrials = numTrials + 1;
                tv.taskName(i) = taskName;
                tv.stimText(i) = stimText;
                tv.prodText(i) = prodText;

                % Check for validity of the stimulus text
                % Note this is the "product" operation, i.e. logical AND,
                % not "production"
                validPrompt = prod([(length(estim.label) == 4), strcmp(estim.label{4}, "equals"), (~isempty(regexp(estim.label{2}, "plus|minus|times",'once')))]);
                validResponse = false;

                % For concatenation we should always use a hyphen when
                % annotating.

                % Determine stimulus type
                switch estim.label{2}
                    case "plus"
                        if(ismember(estim.label{1}, numberMap.keys) && ismember(estim.label{3}, numberMap.keys))
                            arithmetic_operation = "addition";
                            if(~isempty(eprod.label))
                                validResponse = ismember(eprod.label{end}, numberMap.keys);
                            else
                                validResponse = false;
                            end
                        elseif((length(estim.label{1}) == 1) && (length(estim.label{3}) == 1))
                            arithmetic_operation = "concatenation";
                            if(~isempty(prodText))
                                validResponse = length(prodText) == 3;
                                    % response_parts{end} should have the form
                                    % 'X-Y'
                            else
                                validResponse = false;
                            end
                        else
                            validPrompt = false;
                        end
                    case "minus"
                        if(ismember(estim.label{1}, numberMap.keys) && ismember(estim.label{3}, numberMap.keys))
                            arithmetic_operation = "subtraction";
                            if(~isempty(eprod.label))
                                validResponse = ismember(eprod.label{end}, numberMap.keys);
                            else
                                validResponse = false;
                            end
                        else
                            validPrompt = false;
                        end
                    case "times"
                        if(ismember(estim.label{1}, numberMap.keys) && ismember(estim.label{3}, numberMap.keys))
                            arithmetic_operation = "multiplication";
                            eprod.GetParentLabel
                            validResponse = ismember(eprod.label{end}, numberMap.keys);
                        else
                            validPrompt = false;
                        end
                    otherwise
                        validPrompt = false;
                        validResponse = false;
                end

                tv.trialNum(i) = numTrials;
                tv.recId(i) = string(se.userData.expInfo.recId);
                tv.arithmetic_validPrompt(i) = validPrompt;
                tv.arithmetic_validResponse(i) = validResponse;


                if(strcmp(arithmetic_operation, "concatenation"))
                    if(validPrompt)
                        tv.arithmetic_firstOperand(i) = estim.label{1};
                        tv.arithmetic_secondOperand(i) = estim.label{3};
                        tv.arithmetic_operation(i) = arithmetic_operation;
                        tv.arithmetic_trueResult(i) = string([estim.label{1} ' ' estim.label{3}]);
                        if(validResponse)
                            tv.arithmetic_providedResult(i) = prodText;
                            tv.arithmetic_isCorrect(i) = strcmp(tv.arithmetic_providedResult(i), tv.arithmetic_trueResult(i));
                        else
                            tv.arithmetic_providedResult(i) = NaN;
                            tv.arithmetic_isCorrect(i) = false;
                        end
                    end
                    
                else
                    if(validPrompt)
                        tv.arithmetic_firstOperand(i) = numberMap(estim.label{1});
                        tv.arithmetic_secondOperand(i) = numberMap(estim.label{3});
                        tv.arithmetic_operation(i) = arithmetic_operation;
                        switch arithmetic_operation
                                case "addition"
                                    tv.arithmetic_trueResult(i) = tv.arithmetic_firstOperand(i) + tv.arithmetic_secondOperand(i);
                                case "subtraction"
                                    tv.arithmetic_trueResult(i) = tv.arithmetic_firstOperand(i) - tv.arithmetic_secondOperand(i);                            
                                case "multiplication"
                                    tv.arithmetic_trueResult(i) = tv.arithmetic_firstOperand(i) * tv.arithmetic_secondOperand(i);
                        end
                        if(validResponse)
                            tv.arithmetic_providedResult(i) = numberMap(eprod.label{end}); 
                            tv.arithmetic_isCorrect(i) = (tv.arithmetic_providedResult(i) == tv.arithmetic_trueResult(i));
                        else
                            tv.arithmetic_providedResult(i) = NaN;
                            tv.arithmetic_isCorrect(i) = false;                            
                        end
                    end
                end
            end

            % tv.trialNum = height(tt);
            
            if ~numTrials
                fprintf("This se does not have any Arithmetic trial\n");
                return
            end
            
            % Make unique numeric stim IDs
            m = startsWith(tv.taskName, "arithmetic");
            [~, ~, tv.stimNumId(m)] = unique(tv.stimText(m), 'stable');
            
            % Add taskValue table to se
            se.SetTable('taskValue', tv, 'eventValues');
        end

        function [tt, tv] = FindKeyTimes(tt, tv)
            % Find the onset and offset times of all relevant words in the
            % prompt and response
            %
            %   [tt, tv] = FindKeywordTime(tt, tv)
            % 

            relevant_fields = {'arithmetic_firstOperandOn', 'arithmetic_firstOperandOff', 'arithmetic_operationOn', 'arithmetic_operationOff', 'arithmetic_secondOperandOn', 'arithmetic_secondOperandOff', 'arithmetic_equalsOn', 'arithmetic_equalsOff', 'arithmetic_responseOn', 'arithmetic_responseOff'};
            for i = 1:length(relevant_fields)
                tt.(relevant_fields{i}) = NaN(height(tt), 1);
            end
            
            for i = 1 : height(tt)
                estim = tt.stim(i);
                eprod = tt.prod{i};

                if(length(eprod) > 1)
                    eprod = eprod(end);
                    warning('Unexpected prod event length');
                end

                % Add timing of keywords
                % TODO: eventually, could do this at phonetic uniqueness
                % point.
                if(tv.arithmetic_validPrompt(i))
                    tt.arithmetic_firstOperandOn(i) = estim.T.T1(1);
                    tt.arithmetic_firstOperandOff(i) = estim.T.T2(1);
                    tt.arithmetic_secondOperandOn(i) = estim.T.T1(3);
                    tt.arithmetic_secondOperandOff(i) = estim.T.T2(3);
                    tt.arithmetic_operationOn(i) = estim.T.T1(2);
                    tt.arithmetic_operationOff(i) = estim.T.T2(2);
                    tt.arithmetic_equalsOn(i) = estim.T.T1(4);
                    tt.arithmetic_equalsOff(i) = estim.T.T2(4);
                    if(tv.arithmetic_validResponse(i))
                        % TODO: add some nuance for the multiplication and
                        % concatenation cases here
                        tt.arithmetic_responseOn(i) = eprod.T.T1(end);
                        tt.arithmetic_responseOff(i) = eprod.T.T2(end);
                        if(tt.arithmetic_responseOn(i) < tt.arithmetic_equalsOff(i))
                            tt.arithmetic_responseOn(i) = tt.arithmetic_equalsOff(i) + .1;
                            tt.arithmetic_responseOff(i) = tt.arithmetic_equalsOff(i) + .15;
                        end
                    else
                        % Need to do something with these for cases where
                        % there is no response at all.

                        % Just deleting these rows from both tables for now
                        tt.arithmetic_responseOn(i) = tt.arithmetic_equalsOff(i) + .1;
                        tt.arithmetic_responseOff(i) = tt.arithmetic_equalsOff(i) + .15;
                    end
                end
            end
        end
        
        function seSort = SortTrials(se)
            % Sort trials by categories
            % 
            %   seSort = NP.Semsr.SortTrials(se)
            % 
            
            tv = se.GetTable('taskValue');
            
            
            % Set question categories
            operationList = ["addition", "subtraction", "concatenation", "multiplication"];

            tv.operation_with_secondOperand = arrayfun(@(a,b,c) strjoin([string(a), string(b)], ' '), tv.arithmetic_operation, tv.arithmetic_secondOperand);
            operation_with_secondOperand_list = unique(tv.operation_with_secondOperand);

            % tv.operation_with_trueResult = arrayfun(@(a,b,c) strjoin([string(char(a)), string(char(b))], ' '), tv.arithmetic_operation, tv.arithmetic_trueResult);
            % operation_with_trueResult_list = unique(tv.operation_with_trueResult);
            % tv.operation_with_trueResult_categorical = categorical(tv.operation_with_trueResult, operation_with_trueResult_list, 'Ordinal', true);

            tv.operation_with_secondOperand_categorical = categorical(tv.operation_with_secondOperand, operation_with_secondOperand_list, 'Ordinal', true);
            tv.operation_categorical = categorical(tv.arithmetic_operation, operationList, 'Ordinal', true);
            
            secondOperand_list = unique(tv.arithmetic_secondOperand);
            tv.secondOperand_categorical = categorical(tv.arithmetic_secondOperand, secondOperand_list, 'Ordinal', true);

            tv.arithmetic_trueResult(isnan(tv.arithmetic_trueResult)) = -1;
            trueResult_list = unique(tv.arithmetic_trueResult);
            tv.trueResult_categorical = categorical(tv.arithmetic_trueResult, trueResult_list, 'Ordinal', true);
            
            firstOperand_list = unique(tv.arithmetic_firstOperand);
            tv.firstOperand_categorical = categorical(tv.arithmetic_firstOperand, firstOperand_list, 'Ordinal', true);

            % Sort se
            % [tv, I] = sortrows(tv, {'arithmetic_secondOperand'});
            [tv, I] = sortrows(tv, {'operation_categorical'});
            seSort = se.Duplicate;
            seSort.SortEpochs(I);
            seSort.SetTable('taskValue', tv);
        end
        
        function PlotOverview(se, varargin)
            % Plot session heatmap
            
            % Determine units to plot
            ud = se.userData;
            ut = NP.Unit.GetClusTb(se);
            if numel(varargin) == 1
                % By unit ID
                uId = varargin{1};
                uInd = NP.Unit.ClusId2Ind(uId, se);
                uInd = sort(uInd);
                p = 0;
                upp = 12;
            elseif numel(varargin) == 2
                % By page index and # of units per page
                [p, upp] = varargin{:};
                I = (1:upp)+upp*(p-1);
                I(I > height(ut)) = []; % remove out of range indices
                uInd = I;
            end
            uInd = uInd(:)';
            uId = NP.Unit.Ind2ClusId(uInd, se);
            
            % Configure plots
            nr = 3; % number of rows
            ppr = ceil(upp/nr) + 1; % plot per row = #unit / 3 rows + one label plot
            panels1 = repmat({'rate_map'}, [nr ppr]);
            panels1args = repmat({{}}, [nr ppr]);
            
            panels2 = repmat({'event_boundary'}, size(panels1));
            % panels2args = repmat({{{'arithmetic_operationOff', 'arithmetic_secondOperandOff', 'arithmetic_equalsOff', 'arithmetic_responseOn'}, @(x) brewermap(x, 'Dark2')}}, size(panels1));
            % panels2args = repmat({{{'arithmetic_operationOff', 'arithmetic_secondOperandOff', 'arithmetic_responseOn'}, @(x) brewermap(x, 'Dark2')}}, size(panels1));
            % panels2args = repmat({{{'arithmetic_firstOperandOn','arithmetic_firstOperandOff','arithmetic_operationOn','arithmetic_operationOff', 'arithmetic_secondOperandOn','arithmetic_secondOperandOff', 'arithmetic_equalsOn','arithmetic_equalsOff','arithmetic_responseOn'}, @(x) brewermap(x, 'Dark2')}}, size(panels1));
            panels2args = repmat({{{'arithmetic_firstOperandOn','arithmetic_firstOperandOff','arithmetic_operationOff','arithmetic_secondOperandOff','arithmetic_equalsOff','arithmetic_responseOn'}, @(x) brewermap(x, 'Dark2')}}, size(panels1));
            
            panels3 = repmat({'block_boundary'}, size(panels1));
            % panels3args = repmat({{{'arithmetic_secondOperand'},'Label', false}}, size(panels1));
            panels3args = repmat({{{'arithmetic_operation'},'Label', false}}, size(panels1));
            
            k = 1;
            for i = 1 : size(panels1,1)
                % first column for labels
                panels1{i,1} = ''; % not plotting spikes in the first column
                panels2{i,1} = ''; % not plotting event_boundary in the first column
                panels3args{i,1}{3} = true; % show labels in the fir
                
                % following columns for rasters
                for j = 2 : size(panels1,2)
                    if k > numel(uInd)
                        panels1{i,j} = '';
                        panels2{i,j} = '';
                        panels3{i,j} = '';
                        continue
                    end
                    panels1args{i,j} = {uInd(k)};
                    k = k + 1;
                end
            end
            
            tWin = [-2 4.5];
            
            % Plotting
            f = MPlot.Figure(38050301); clf
            f.WindowState = 'maximized';
            NP.Arithmetic.PlotSE(se, panels1, 'PanelArgs', panels1args, 'TimeWindow', tWin);
            NP.Arithmetic.PlotSE(se, panels2, 'PanelArgs', panels2args, 'TimeWindow', tWin);
            NP.Arithmetic.PlotSE(se, panels3, 'PanelArgs', panels3args, 'TimeWindow', tWin);
            
            % Save figure
            recId = string(NP.SE.GetID(ud));
            region = NP.SE.GetRegion(ud(1));
            figDir = fullfile(NP.Data.GetAnalysisRoot, "se_"+region, recId);
            if ~exist(figDir, 'dir')
                mkdir(figDir);
            end
            % figName = strjoin([recId, "semsr", "page"+p, "u"+uId], "_");
            figName = strjoin([recId, "arithmetic", "page"+p], "_"); %filename is often too long with all the unitIds included
            if ~exist(figDir, 'dir')
                mkdir(figDir);
            end
            exportgraphics(f, fullfile(figDir, figName+".png"));
        end

    end
end

