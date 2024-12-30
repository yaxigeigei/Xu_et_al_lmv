classdef Decode
    methods(Static)
        % Classification
        function xrBoot = ClassifyXT(sData, varargin)
            % Compute cross-time sentence classification accuracy
            % 
            %   xrBoot = ClassifyXT(sData)
            %   xrBoot = ClassifyXT(sData, ..., 'NBoot', 100)
            %   xrBoot = ClassifyXT(sData, ..., 'Shuffle', false)
            %   xrBoot = ClassifyXT(sData, ..., 'ErasePrint', true)
            % 
            
            p = inputParser;
            p.addParameter("NBoot", 100, @(x) isnumeric(x) && isscalar(x));
            p.addParameter("Shuffle", false, @(x) islogical(x) && isscalar(x));
            p.addParameter("ErasePrint", ~isunix, @(x) islogical(x) && isscalar(x));
            p.parse(varargin{:});
            nBoot = p.Results.NBoot;
            isShuffle = p.Results.Shuffle;
            isErase = p.Results.ErasePrint;
            
            % Get dimensions
            s = sData(1);
            nSen = size(s.indTest, 1);
            nTime = size(s.X, 3);
            
            % Make trial labels
            yTest = (1:nSen)';
            yTrain = repelem(yTest, size(s.indTrain,1)/nSen);
            
            % Bootstrap
            xrBoot = zeros(nTime, nTime, nBoot);
            for n = 1 : nBoot
                % Construct training and testing arrays
                xTrain = arrayfun(@(x) x.X(x.indTrain(:,n), x.clusTb.isUnit, :), sData, 'Uni', false);
                xTest = arrayfun(@(x) x.X(x.indTest(:,n), x.clusTb.isUnit, :), sData, 'Uni', false);
                xTrain = cat(2, xTrain{:});
                xTest = cat(2, xTest{:});
                
                % Optionally shuffle training labels
                if isShuffle
                    yTrain = randsample(yTrain, numel(yTrain));
                end
                
                % Compute cross-time classification accuracies
                str = sprintf('Iteration %i/%i\n', n, nBoot);
                fprintf(str);
                parfor a = 1 : nTime
                    for b = 1 : nTime
                        mdl = fitcecoc(xTrain(:,:,a), yTrain, 'Coding', 'onevsall');
                        xrBoot(a,b,n) = 1 - loss(mdl, xTest(:,:,b), yTest);
                    end
                end
                if isErase && n < nBoot
                    fprintf(repmat('\b', [1 numel(str)]));
                end
            end
        end
        
        function sData = ComputeXTStats(sData)
            % 
            
            for i = 1 : numel(sData)
                s = sData{i};
                [s.mxr, s.mxrCI] = iBootError(s.xrBoot);
                [s.mxrNull, s.mxrNullCI] = iBootError(s.xrNullBoot);
                sData{i} = s;
            end
            
            function [mAll, mCI, mDiagAll, mDiagCI] = iBootError(xr)
                % Overall mean
                mAll = mean(xr, 3);
                mDiagAll = LMV.Decode.IGetDiag(mAll);
                
                % Bootstrap mean
                nBoot = 250;
                nSp = 8;
                nIter = size(xr,3);
                for k = nBoot : -1 : 1
                    ind = randsample(nIter, nSp);
                    mBoot(:,:,k) = mean(xr(:,:,ind), 3);
                end
                mCI = prctile(mBoot, [2.5 97.5], 3);
                mDiagCI = LMV.Decode.IGetDiag(mCI);
            end
        end
        
        function PlotAccuracyHeatmap(s, ce)
            % 
            % 
            %   PlotAccuracyHeatmap(s, ce)
            % 
            
            ax = gca;
            
            % Plot heatmap
            M = s.mxr;
            % M(M < s.mxrNullCI(:,:,2)) = NaN;
            t = ce.GetTable("resp").time{1}(5:10:end);
            imagesc(t, t, M); hold on
            
            % Plot task phase markers
            tt = ce.GetTable("taskTime");
            wStim = [tt.stimMatchOn(1) tt.stimMatchOff(1)];
            wProd = [tt.prodMatchOn(1) tt.prodMatchOff(1)];
            plot([wStim; wStim], ax.YLim', 'Color', LMV.Param.GetTaskPhaseColors("stim"));
            plot(ax.YLim', [wStim; wStim], 'Color', LMV.Param.GetTaskPhaseColors("stim"));
            plot([wProd; wProd], ax.YLim', 'Color', LMV.Param.GetTaskPhaseColors("prod"));
            plot(ax.YLim', [wProd; wProd], 'Color', LMV.Param.GetTaskPhaseColors("prod"));
            
            b = colorbar;
            b.Label.String = "Accuracy (frac.)";
            
            axis(ax, 'square', 'equal', 'tight');
            ax.CLim(1) = 0.25;
            % ax.CLim(2) = 0.55;
            ax.YLabel.String = "Trained from (s)";
            ax.XLabel.String = "Tested on (s)";
            ax.Title.String = sprintf("%s %s", s.regionName, s.groupName);
            MPlot.Axes(ax);
        end
        
        function PlotAccuracyTimeseries(s, ce, perfType)
            % 
            
            % Find classification accuracies
            switch perfType
                case "diag"
                    m = LMV.Decode.IGetDiag(s.mxr);
                    ci = LMV.Decode.IGetDiag(s.mxrCI);
                    
                    mNull = LMV.Decode.IGetDiag(s.mxrNull);
                    ciNull = LMV.Decode.IGetDiag(s.mxrNullCI);
                    
                case "bestTest"
                    [m, I] = max(s.mxr, [], 1);
                    I = sub2ind(size(s.mxr), I, 1:numel(I));
                    
                    mxrCI = reshape(s.mxrCI, [], 2);
                    ci = mxrCI(I,:);
                    
                    mNull = s.mxrNull(I);
                    mxrNullCI = reshape(s.mxrNullCI, [], 2);
                    ciNull = mxrNullCI(I,:);
                    
                case "bestTrain"
                    [m, I] = max(s.mxr, [], 2);
                    I = sub2ind(size(s.mxr), 1:numel(I), I');
                    
                    mxrCI = reshape(s.mxrCI, [], 2);
                    ci = mxrCI(I,:);
                    
                    mNull = s.mxrNull(I);
                    mxrNullCI = reshape(s.mxrNullCI, [], 2);
                    ciNull = mxrNullCI(I,:);
                    
                otherwise
                    error("'%s' is not a valid type of performance.", perfType);
            end
            
            ax = gca;
            
            % Plot traces
            cLk = LMV.Param.GetLinkerColors(s.groupName);
            t = ce.GetTable("resp").time{1}(5:10:end);
            % plot(t([1 end]), [0.25 0.25], '--', 'Color', [0 0 0]+.7); hold on
            plot(t, m, 'Color', cLk); hold on
            MPlot.ErrorShade(t, m, ci(:,1), ci(:,2), 'IsRelative', false, 'Color', cLk);
            plot(t, mNull, 'k');
            MPlot.ErrorShade(t, mNull, ciNull(:,1), ciNull(:,2), 'IsRelative', false);
            
            % Plot task phases
            tt = ce.GetTable("taskTime");
            wStim = [tt.stimMatchOn(1) tt.stimMatchOff(1)];
            wProd = [tt.prodMatchOn(1) tt.prodMatchOff(1)];
            plot(tt.cue1([1 1]), ax.YLim, 'Color', LMV.Param.GetTaskPhaseColors("none"));
            plot([wStim; wStim], ax.YLim', 'Color', LMV.Param.GetTaskPhaseColors("stim"));
            plot(tt.cue3([1 1]), ax.YLim, 'Color', LMV.Param.GetTaskPhaseColors("none"));
            plot([wProd; wProd], ax.YLim', 'Color', LMV.Param.GetTaskPhaseColors("prod"));
            
            ax.XLim = t([1 end]);
            ax.YLim(1) = 0;
            % ax.XLabel.String = "Aligned time (s)";
            ax.YLabel.String = "Accuracy";
            ax.Title.String = sprintf("%s %s", s.regionName, s.groupName);
            MPlot.Axes(ax);
        end
        
        function a = IGetDiag(A)
            [nTime, ~, nPage] = size(A);
            I = logical(repmat(eye(nTime), [1 1 nPage]));
            a = squeeze(reshape(A(I), [nTime nPage]));
        end
        
        function mdlTb = ClassifyTrials(ce, mdlName)
            % Fit trial classification models
            % 
            %   mdlTb = ClassifyTrials(ce, mdlName)
            % 
            
            % Configure model fitting
            fprintf("\n%s\n", NP.SE.GetID(ce));
            
            phaseType = mdlName(isletter(mdlName));
            switch phaseType
                case 'speech'
                    phaseNames = {'stim', 'prod'};
                case 'delay-init'
                    phaseNames = {'delay', 'init'};
                otherwise
                    phaseNames = cellstr(phaseType);
            end
            
            
            switch mdlName(~isletter(mdlName))
                case '4'
                    stimIds = cellstr(LMV.Param.stimIdList4);
                case '14'
                    stimIds = cellstr(LMV.Param.stimIdList14);
                otherwise
                    stimIds = cellstr(LMV.Param.stimIdList14);
            end
            
            
            % Get predictors
            [~, X] = ce.GetArray('resp', 'Normalization', 'zscore');
            
            % Get class labels
            [~, Y] = ce.GetArray('feat', [], stimIds);
            assert(all(sum(Y,2)<=1), "Found more than one class label in the same observation(s)");
            Y = Y .* (1 : size(Y,2));
            y = sum(Y, 2);
            y = categorical(y, 1:numel(stimIds), stimIds);
            
            % Get task phase and trial mask
            [~, M] = ce.GetArray('feat', [], phaseNames);
            m = any(M,2) & ~isundefined(y);
            m(1:2:end) = false; % 2X downsampling
            
            
            % Search hyperparams
            fprintf("\nSearch hyperparams for trial classification.\n");
            isLocal = ispc || ismac;
            learner = templateLinear('Regularization', 'ridge', 'Solver', 'lbfgs');
            hp = {'Lambda'};
            hpOpt = struct('ShowPlot', false, 'Verbose', double(isLocal), 'Kfold', 10);
            ss = statset('UseParallel', ~isLocal);
            mdlHS = NP.Fit.BatchFitcecoc(X(m,:), y(m), ...
                'Coding', 'onevsall', ...
                'Learners', learner, ...
                'OptimizeHyperparameters', hp, ...
                'HyperparameterOptimizationOptions', hpOpt, ...
                'PredictorNames', ce.respNames, 'Options', ss);
            L = mdlHS{1}.BinaryLearners{1}.Lambda;
            
            % Fit model with optimal hyperparam
            fprintf("\nFit trial classification model with optimal lambda %.1e.\n", L);
            learner = templateLinear('Regularization', 'ridge', 'Solver', 'lbfgs', 'Lambda', L);
            [mdl, nullData] = NP.Fit.BatchFitcecoc(X(m,:), y(m), ...
                'NumShifts', 500, 'MinShift', round(0.1*sum(m)), ...
                'Coding', 'onevsall', ...
                'Learners', learner, ...
                'KFold', 10, ...
                'PredictorNames', ce.respNames, 'Options', ss);
            
            
            % Construct output table
            mdlTb = table;
            mdlTb.dt = 0;
            mdlTb.(mdlName+"HS") = mdlHS;
            mdlTb.(mdlName) = mdl;
            mdlTb.(mdlName+"Null") = nullData;
        end
        
        function mdlTb = ClassifyPhases(ce, mdlName)
            % Fit task phase classification models
            % 
            %   mdlTb = ClassifyPhases(ce, mdlName)
            % 
            
            % Configure model fitting
            fprintf("\n%s\n", NP.SE.GetID(ce));
            
            % Get predictors
            [~, X] = ce.GetArray('resp', 'Normalization', 'zscore');
            
            % Get class labels
            phaseNames = {'atten', 'stim', 'delay', 'init', 'prod', 'iti'};
            [~, Y] = ce.GetArray('feat', [], phaseNames);
            Y(sum(Y,2)>1, :) = 0; % exclude overlapping phases
            Y = Y .* (1 : size(Y,2));
            y = sum(Y, 2);
            y = categorical(y, 1:numel(phaseNames), phaseNames);
            
            % Get mask
            m = ~isundefined(y);
            m(1:2:end) = false; % 2X downsampling
            
            
            % Search hyperparams
            fprintf("\nSearch hyperparams for task phase classification.\n");
            isLocal = ispc || ismac;
            learner = templateLinear('Regularization', 'ridge', 'Solver', 'lbfgs');
            hp = {'Lambda'};
            hpOpt = struct('ShowPlot', isLocal, 'Verbose', double(isLocal), 'Kfold', 10);
            ss = statset('UseParallel', ~isLocal);
            mdlHS = NP.Fit.BatchFitcecoc(X(m,:), y(m), ...
                'Learners', learner, ...
                'OptimizeHyperparameters', hp, ...
                'HyperparameterOptimizationOptions', hpOpt, ...
                'PredictorNames', ce.respNames, 'Options', ss);
            L = mdlHS{1}.BinaryLearners{1}.Lambda;
            
            % Fit model with optimal hyperparam
            fprintf("\nFit task phase classification model with optimal lambda %.1e.\n", L);
            learner = templateLinear('Regularization', 'ridge', 'Solver', 'lbfgs', 'Lambda', L);
            [mdl, nullData] = NP.Fit.BatchFitcecoc(X(m,:), y(m), ...
                'NumShifts', 500, 'MinShift', round(0.1*sum(m)), ...
                'Learners', learner, ...
                'PredictorNames', ce.respNames, 'Options', ss);
            
            
            % Construct output table
            mdlTb = table;
            mdlTb.dt = 0;
            mdlTb.(mdlName+"HS") = mdlHS;
            mdlTb.(mdlName) = mdl;
            mdlTb.(mdlName+"Null") = nullData;
        end
        
        function mdlTb = LoadClassModels(mdlDir, recId, mdlName)
            % Load and organize classification models into a table
            % 
            %   mdlTb = LoadClassModels(mdlDir)
            %   mdlTb = LoadClassModels(mdlDir, recId)
            %   mdlTb = LoadClassModels(mdlDir, recId, mdlName)
            % 
            
            % Find model files
            if ~exist('recId', 'var') || isempty(recId)
                recId = '*';
            end
            if ~exist('mdlName', 'var') || isempty(mdlName)
                mdlName = '*';
            end
            
            % Load models
            mdlSearch = MBrowse.Dir2Table(fullfile(mdlDir, [recId '_mdlTb_' mdlName '.mat']));
            mdlTbs = cell(height(mdlSearch), 1);
            for i = 1 : numel(mdlTbs)
                load(fullfile(mdlSearch.folder{i}, mdlSearch.name{i}), 'mdlTb');
                
                recId = regexp(mdlSearch.name{i}, '^.+(?=_mdlTb)', 'match', 'once');
                mdlName = regexp(mdlSearch.name{i}, '(?<=mdlTb_).+(?=.mat$)', 'match', 'once');
                
                mdlTb.recId{1} = recId;
                mdlTb.mdlName{1} = mdlName;
                mdlTb.Properties.VariableNames = strrep(mdlTb.Properties.VariableNames, mdlName, 'mdl');
                mdlTbs{i} = mdlTb;
            end
            mdlTb = cat(1, mdlTbs{:});
            
            % Derive stats
            ciPrct = [0 100]+[1 -1]*5/2/height(mdlTb);
            
            for k = 1 : height(mdlTb)
                mdl = mdlTb.mdl{k};
                null = mdlTb.mdlNull{k};
                
                % Gross accuracy
                mdlTb.p(k) = 1 - mdl.kfoldLoss;
                mdlTb.pNull(k) = mean(1-null.kfoldLoss);
                mdlTb.pNullCI(k,:) = prctile(1-null.kfoldLoss, ciPrct);
                
                % Confusion matrix
                cm = confusionmat(mdl.Y, mdl.kfoldPredict);
                cmNullMean = mean(null.confMat, 3);
                cmNullCI = prctile(null.confMat, ciPrct, 3);
                cmSig = cm < cmNullCI(:,:,1) | cm > cmNullCI(:,:,2);
                
                mdlTb.cm{k} = cm;
                mdlTb.cmNullMean{k} = cmNullMean;
                mdlTb.cmNullCI{k} = cmNullCI;
                mdlTb.cmSig{k} = cmSig;
                
                % Accuracy by class
                cm = cm ./ sum(cm,2);
                pCla = diag(cm);
                cmNull = null.confMat;
                cmNull = cmNull ./ sum(cmNull,2);
                pClaNull = zeros(numel(pCla), null.nShifts);
                for n = 1 : null.nShifts
                    pClaNull(:,n) = diag(cmNull(:,:,n));
                end
                pClaNullCI = prctile(pClaNull, ciPrct, 2);
                
                mdlTb.pCla{k} = pCla;
                mdlTb.pClaNull{k} = mean(pClaNull, 2);
                mdlTb.pClaNullCI{k} = pClaNullCI;
            end
        end
        
        function PlotGrossAccuracy(mdlTbs)
            % Plot overall classification accuracy from tables of model data
            % Each table makes one subplot. x-axis is different models. y-axis is accuracy.
            % 
            %   PlotGrossAccuracy(mdlTbs)
            % 
            % Input
            %   mdlTbs      1) A table of model data. Each row is for a model.
            %               2) A cell array of such tables. Each table makes one subplot.
            %               Data in the following table columns will be used:
            %               p           Overall classification accuracy.
            %               pNullCI     Confidence interval of the null accuracy.
            %               mdl         Cell array of model objects.
            %               mdlName     Model names in a cell array of char strings.
            %               recId       Recording IDs in a cell array of char strings.
            % 
            
            if istable(mdlTbs)
                mdlTbs = {mdlTbs};
            end
            
            tiledlayout('flow');
            
            for k = 1 : numel(mdlTbs)
                
                mTb = mdlTbs{k};
                mdlNames = mTb.mdlName;
                nClass = numel(mTb.mdl{1}.ClassNames);
                
                y = mTb.p;
                ci = mTb.pNullCI;
                x = (1:numel(y))';
                
                ax = nexttile;
                MPlot.ErrorShade(x, y, ci(:,1), ci(:,2), 'IsRelative', false); hold on
                plot(x, y, 'ko-');
                plot(x, ones(size(x))/nClass, '--', 'Color', [0 0 0 .3]);
                ax.XLim = [0.5 numel(x)+0.5];
                ax.XTick = x;
                ax.XTickLabel = mdlNames;
                ax.YLim = [0 1];
                ax.YTick = 0:.2:1;
                ax.YLabel.String = "Accuracy";
                ax.Title.String = strrep(mTb.recId{1}, '_', ' ');
                ax.Title.Interpreter = 'none';
            end
        end
        
        function PlotClassAccuracy(mdlTb)
            % Plot classification accuracy by class
            % 
            %   PlotClassAccuracy(mdlTb)
            % 
            % Inputs
            %   mdlTb       A table of model data. Each row is for a model.
            %               Data in the following table columns will be used:
            %               recId       Recording IDs in a cell array of char strings.
            %               mdl         Model object stored in the column named with mdlName.
            %               pCla        Classification accuracies by class.
            %               pClaNullCI  Confidence intervals of the null accuracies.
            % 
            
            tiledlayout('flow');
            
            for k = 1 : height(mdlTb)
                y = mdlTb.pCla{k};
                ci = mdlTb.pClaNullCI{k};
                x = (1 : numel(y))';
                
                ax = nexttile;
                MPlot.ErrorShade(x, y, ci(:,1), ci(:,2), 'IsRelative', false); hold on
                plot(x, y, 'ko-');
                ax.XLim = [0.5 numel(x)+0.5];
                ax.XTick = x;
                ax.XTickLabel = mdlTb.mdl{k}.ClassNames;
                ax.YLim = [0 1];
                ax.YLabel.String = "Accuracy";
                ax.Title.String = strrep(mdlTb.recId{k}, '_', ' ');
                ax.Title.Interpreter = 'none';
            end
        end
        
        function PlotConfusionMat(mdlTb, varargin)
            % Plot confusion matrices from a table of model data
            % 
            %   PlotConfusionMat(mdlTb)
            %   PlotConfusionMat(..., 'Diff', false)
            %   PlotConfusionMat(..., 'SigOnly', false)
            %   PlotConfusionMat(..., 'ClassNames', [])
            %   PlotConfusionMat(..., 'Layout', tiledlayout('flow'))
            % 
            % Inputs
            %   mdlTb       A table of model data. Each row is for a model.
            %               Data in the following table columns will be used:
            %               recId       Recording IDs in a cell array of char strings.
            %               mdlName     Model names in a cell array of char strings.
            %               cm          2D numeric arrays of confusion matrix.
            %               cmSig       2D logical arrays where the significant elements in cm is true.
            %               mdl         Model object stored in the column named with mdlName.
            %   'Diff'      If true, plot the confusion matrix after subtracting the null (cmNullMean).
            %   'SigOnly'   If true, only plot significant elements. Other elements are black.
            %   'Layout'    TiledChartLayout object to continue plotting with.
            % 
            
            p = inputParser;
            p.addParameter('Diff', false, @islogical);
            p.addParameter('SigOnly', false, @islogical);
            p.addParameter('ClassNames', []);
            p.addParameter('Layout', [], @(x) isa(x, 'matlab.graphics.layout.TiledChartLayout'));
            p.parse(varargin{:});
            claNames = p.Results.ClassNames;
            t = p.Results.Layout;
            isDiff = p.Results.Diff;
            isSigOnly = p.Results.SigOnly;
            
            if isempty(t)
                tiledlayout('flow');
            end
            
            for k = 1 : height(mdlTb)
                cm = mdlTb.cm{k};
                cm = cm ./ sum(cm,2) * 100;
                
                if isDiff
                    cmn = mdlTb.cmNullMean{k};
                    cmn = cmn ./ sum(cmn,2) * 100;
                    cm = cm - cmn;
                end
                
                if isSigOnly
                    cm(~mdlTb.cmSig{k}) = NaN;
                end
                
                if isempty(claNames)
                    claNames = mdlTb.mdl{k}.ClassNames;
                end
                
                nexttile;
                h = heatmap(claNames, claNames, round(cm), 'CellLabelFormat', '%d');
                h.Title = [strrep(mdlTb.recId{k}, '_', ' ') ', ' mdlTb.mdlName{k}];
                h.ColorLimits = [0 40];
                h.ColorbarVisible = 'off';
            end
        end
        
        % Regression
        function mdlTb = FitSpeechDecoding(ce, mdlName, varargin)
            % Fit speech feature decoding models
            % 
            %   mdlTb = FitSpeechDecoding(ce, mdlName)
            %   mdlTb = FitSpeechDecoding(..., 'Lambda', [])
            % 
            % Inputs
            %   ce          An NP.CodingExplorer object.
            %   mdlName     Name of the model to fit.
            %   'Lambda'    If empty [], will search the parameter space by fitting models with a series of lambda values.
            %               If provided, this lambda value will be use to fit all models, and models with randomly shifted 
            %               data to bootstrap the null distribution of model weights.
            % Output
            %   mdlTb       
            % 
            
            p = inputParser;
            p.addParameter('Lambda', [], @isnumeric);
            p.parse(varargin{:});
            L = p.Results.Lambda;
            
            % Configure model fitting
            fprintf("\n%s\n\n", NP.SE.GetID(ce));
            ops = ce.userData.ops;
            feats = [{'peakEnv', 'peakRate'} ops.featVars.pitch ops.featVars.artic];
            switch mdlName
                case 'stim'
                    dt = -0.4 : 0.05 : 0;
                    phaseName = 'stim';
                case 'prod'
                    dt = -0.4 : 0.05 : 0.4;
                    phaseName = 'prod';
                case 'test'
                    dt = 0 : 0.25 : 0.5;
                    phaseName = 'prod';
                otherwise
                    error("'%s' is not a valid mdlName", mdlName);
            end
            
            % Fit models
            mdls = cell(numel(dt), numel(feats));
            for i = 1 : numel(dt)
                % Get predictors
                [tx, X] = ce.GetArray('resp', 'TimeShifts', dt(i), 'Normalization', 'zscore');
                
                % Get responses
                [ty, Y] = ce.GetArray('feat', [], feats, 'Normalization', 'zscore');
                k = find(ty == tx(1), 1);
                Y = Y(k:k+numel(tx)-1,:);
                Y(isnan(Y)) = 0;
                
                % Masking
                masks = ops.featVars.taskTimeSpan;
                [tm, M] = ce.GetArray('feat', [], masks);
                k = find(tm == tx(1), 1);
                M = M(k:k+numel(tx)-1,:);
                maskTb = array2table(logical(M), 'VariableNames', masks);
                m = maskTb.(phaseName);
                X = X(m,:);
                Y = Y(m,:);
                
                % Model fitting
                ss = statset('UseParallel', true);
                if isempty(L)
                    str = sprintf("Search lambda for %s decoding at %i/%i time offset.\n", mdlName, i, numel(dt));
                    fprintf(str);
                    lambdaList = logspace(-3, 2, 30);
                    mdls(i,:) = NP.Fit.BatchLasso(X, Y, 'Verbose', false, ...
                        'Alpha', eps, 'Standardize', false, 'CV', 10, 'Lambda', lambdaList, 'Options', ss); % {'Alpha', eps} is ridge regressions
                else
                    str = sprintf("Fit %s decoding models at %i/%i time offset.\n", mdlName, i, numel(dt));
                    fprintf(str);
                    mdls(i,:) = NP.TRF.BatchLasso(X, Y, 'Verbose', false, ...
                        'Alpha', eps, 'Standardize', false, 'Lambda', L, 'Options', ss, ...
                        'NumShifts', 500, 'MinShift', round(0.1*sum(m)));
                end
                
                if i < numel(dt) && (ispc || ismac)
                    fprintf(repmat('\b', 1, length(str)));
                end
            end
            
            % Construct output table
            mdlTb = table;
            mdlTb.dt = dt';
            mdlTb = [mdlTb cell2table(mdls, 'VariableNames', feats)];
        end
        
        function [r2, sig] = GetRSquared(mdls)
            % Extract r-squared values with consistent regularization in each variable
            % 
            %   [r2, sig] = GetRSquared(mdls)
            % 
            
            % Find the indices of minimal MSE
            [~, ind] = arrayfun(@(x) min(x.mdl.MSE), mdls);
            
            % Find the median index
            ind = ones(size(ind)) .* median(ind, "all");
            
            % Get R-squared values
            r2 = arrayfun(@(x,k) x.r2(k), mdls, ind);
            
            % Get the level of significance
            if ~isfield(mdls(1), 'null')
                sig = [];
                return
            end
            sig = arrayfun(@(x) NP.TRF.GetSigLevel(x.r2, x.null.r2), mdls);
        end
        
        function tb = LoadRegModels(mdlDir, mdlNames)
            % Load models and organize them into a table
            % 
            %   mdlTbs = LoadRegModels(mdlDir, mdlName)
            % 
            
            % Load models
            tb = table;
            mdlNames = cellstr(mdlNames);
            for n = 1 : numel(mdlNames)
                mn = mdlNames{n};
                mdlSearch = MBrowse.Dir2Table(fullfile(mdlDir, ['*_mdlTb_' mn '.mat']));
                for i = height(mdlSearch) : -1 : 1
                    load(fullfile(mdlSearch.folder{i}, mdlSearch.name{i}), 'mdlTb');
                    tb.recId(i) = string(NP.SE.GetID(mdlSearch.name{i}));
                    tb.(mn){i} = mdlTb;
                end
            end
            
%             % Derive stats
%             ciPrct = [0 100]+[1 -1]*1/2;
%             
%             for k = 1 : height(mdlsTb)
%                 mdl = mdlsTb.(mdlName){k};
%                 null = mdlsTb.([mdlName 'Null']){k};
%                 
%                 % Gross accuracy
%                 mdlsTb.p(k) = 1 - mdl.kfoldLoss;
%                 mdlsTb.pNull(k) = mean(1-null.kfoldLoss);
%                 mdlsTb.pNullCI(k,:) = prctile(1-null.kfoldLoss, ciPrct);
%                 
%                 % Confusion matrix
%                 cm = confusionmat(mdl.Y, mdl.kfoldPredict);
%                 cmNullMean = mean(null.confMat, 3);
%                 cmNullCI = prctile(null.confMat, ciPrct, 3);
%                 cmSig = cm < cmNullCI(:,:,1) | cm > cmNullCI(:,:,2);
%                 
%                 mdlsTb.cm{k} = cm;
%                 mdlsTb.cmNullMean{k} = cmNullMean;
%                 mdlsTb.cmNullCI{k} = cmNullCI;
%                 mdlsTb.cmSig{k} = cmSig;
%                 
%                 % Accuracy by class
%                 cm = cm ./ sum(cm,2);
%                 pCla = diag(cm);
%                 cmNull = null.confMat;
%                 cmNull = cmNull ./ sum(cmNull,2);
%                 pClaNull = zeros(numel(pCla), null.nShifts);
%                 for n = 1 : null.nShifts
%                     pClaNull(:,n) = diag(cmNull(:,:,n));
%                 end
%                 pClaNullCI = prctile(pClaNull, ciPrct, 2);
%                 
%                 mdlsTb.pCla{k} = pCla;
%                 mdlsTb.pClaNull{k} = mean(pClaNull, 2);
%                 mdlsTb.pClaNullCI{k} = pClaNullCI;
%             end
        end
        
        function PlotOptimalHyperparam(mdls, metricName)
            % Plot optimal hyperparameter values in STRF style
            % 
            %   PlotOptimalHyperparam(mdls, metricName)
            % 
            % Inputs
            %   metricName          'LambdaMinMSE' or ...
            % 
            
            mdlNames = mdls.Properties.VariableNames(2:end);
            [nc, nr] = size(mdls{:,2:end});
            
            tl = tiledlayout(nr, nc);
            tl.Padding = 'compact';
            
            for j = 1 : nr
                for i = 1 : nc
                    %
                    mn = mdlNames{j};
                    mdlTb = mdls.(mn){i};
                    x = mdlTb.dt;
                    
                    feats = mdlTb.Properties.VariableNames(2:end);
                    y = 1 : numel(feats);
                    
                    v = arrayfun(@(x) x.mdl.(metricName), mdlTb{:,2:end});
                    v = log10(v);
                    v = v';
                    
                    x = -flip(x);
                    v = flip(v, 2);
                    
                    ax = nexttile((j-1)*nc+i);
                    
                    imagesc(ax, x, y, v);
                    colorbar(ax);
                    ax.CLim = [-2 1];
                    ax.XLim = [-1 0.5];
                    ax.XLabel.String = 'Time from feature (s)';
                    ax.YTick = y;
                    ax.YTickLabel = feats;
                    ax.YDir = 'reverse';
                    ax.Title.String = mdls.recId(i) + ', ' + mn;
                    ax.Title.Interpreter = 'none';
                    MPlot.Axes(ax);
                end
            end
        end
        
        function PlotRegPerf(mdls)
            % Plot decoding performance in STRF style
            % 
            %   PlotRegPerf(mdls)
            % 
            % Inputs
            %   
            
            mdlNames = mdls.Properties.VariableNames(2:end);
            [nc, nr] = size(mdls{:,2:end});
            
            tl = tiledlayout(nr, nc);
            tl.Padding = 'compact';
            
            for j = 1 : nr
                for i = 1 : nc
                    %
                    mn = mdlNames{j};
                    mdlTb = mdls.(mn){i};
                    x = mdlTb.dt;
                    
                    feats = mdlTb.Properties.VariableNames(2:end);
                    y = 1 : numel(feats);
                    
                    [v, sig] = LMV.Decode.GetRSquared(mdlTb{:,2:end});
%                     if ~isempty(sig)
%                         v(sig < 1) = 0;
%                     end
                    
                    v = sqrt(max(v, 0));
                    
                    v = v';
                    x = -flip(x);
                    v = flip(v, 2);
                    
                    ax = nexttile((j-1)*nc+i);
                    
                    imagesc(ax, x, y, v);
                    cb = colorbar(ax);
                    cb.Title.String = 'r';
                    ax.CLim(1) = 0.1;
                    ax.XLim = [-.4 .4];
                    ax.XLabel.String = 'Time from feature (s)';
                    ax.YTick = y;
                    ax.YTickLabel = feats;
                    ax.YDir = 'reverse';
                    ax.Title.String = mdls.recId(i) + ', ' + mn;
                    ax.Title.Interpreter = 'none';
                    MPlot.Axes(ax);
                end
            end
        end
        
    end
end

