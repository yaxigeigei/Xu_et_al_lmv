classdef UnitDecode
    
    properties(Constant)
        senMdl = "sen4_lda";
        senXphaseMdl = "sen4x_ecoc";
    end
    
    methods(Static)
        % Classification
        function [mdl, nullData] = Sentences(rTb, fTb, stimIds, phaseName)
            % Fit trial classification models
            % 
            %   mdlTb = Sentences(rTb, fTb, stimIds, phaseName)
            % 
            
            % Get predictors
            X = rTb{:,:};
            
            % Get class labels
            Y = fTb{:,stimIds};
            assert(all(sum(Y,2)<=1), "Found more than one class label in the same observation(s)");
            Y = Y .* (1 : size(Y,2));
            y = sum(Y, 2);
            y = categorical(y, 1:numel(stimIds), stimIds);
            
            % Get task phase mask
            M = fTb{:,phaseName};
            m = any(M,2) & ~isundefined(y);
            
            % Print input info
%             fprintf("\n%s\n", NP.SE.GetID(ce));
            fprintf("%i predictors\n", size(X,2));
            fprintf("%i raw samples, %i masked\n", numel(m), sum(m));
            
            % Fit models
            fprintf("\nFit sentence classification models.\n");
            
            learner = templateLinear();
            isLocal = ispc || ismac;
%             ss = statset('UseParallel', ~isLocal);
            ss = statset('UseParallel', true);
            [mdl, nullData] = LMV.UnitDecode.BatchFitcecoc(X(m,:), y(m), ...
                'Coding', 'onevsone', ...
                'Learners', learner, ...
                'Leaveout', 'on', ...
                'Options', ss);
            
            % [mdl, nullData] = LMV.UnitDecode.BatchFitcdiscr(X(m,:), y(m), 'DiscrimType', 'pseudolinear', 'Leaveout', 'on');
            
        end
        
        function [test, null] = BatchFitcdiscr(X, y, varargin)
            % Run fitcdiscr for each column in X
            % 
            %   [test, null] = BatchFitcdiscr(X, y)
            %   [test, null] = BatchFitcdiscr(X, y, 'NPerm', 0)
            %   [test, null] = BatchFitcdiscr(..., fitcdiscrArgs)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('NPerm', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Verbose', true, @islogical);
            p.parse(varargin{:});
            nPerm = p.Results.NPerm;
            isVerbose = p.Results.Verbose;
            funcArgs = [fieldnames(p.Unmatched) struct2cell(p.Unmatched)]';
            
            nMdls = size(X, 2);
            test = cell(nMdls, 1);
            null = cell(nMdls, 1);
            for i = 1 : nMdls
                if isVerbose
                    str = sprintf('Fitting models %i/%i\n', i, nMdls);
                    fprintf(str);
                end
                
                % Prepare input for this iteration
                x = X(:,i);
                
                % Remove outliers in each group
                C = unique(y);
                for k = 1 : numel(C)
                    m = y == C(k);
                    xClean = x(m);
                    xClean(isoutlier(xClean)) = NaN;
                    x(m) = xClean;
                end
                
                % Compensate class imbalance with sample weights
                w = ones(size(y));
                for n = 1 : numel(C)
                    m = y==C(n) & ~isnan(x); % use valid samples
                    w(m) = 1 / sum(m);
                end
                
                % Fit test model
                try
                    mdl = fitcdiscr(x, y, 'Weights', w, funcArgs{:});
                    s = struct;
                    if isprop(mdl, 'KFold')
                        s.loss = mdl.kfoldLoss;
                        s.confMat = confusionmat(mdl.Y, mdl.kfoldPredict);
                    else
                        s.loss = mdl.loss(x, y);
                        s.confMat = confusionmat(y, mdl.predict(x));
                    end
                    s.classNames = mdl.ClassNames;
                    test{i} = s;
                catch e
                    disp(e);
                    continue
                end
                
                % Fit null models with permutation
                if ~nPerm
                    continue
                end
                
                mLoss = zeros(1, nPerm);
                cm = cell(1, nPerm);
                
                parfor n = 1 : nPerm
                    I = randperm(numel(x));
                    mdl = fitcecoc(x(I), y, 'Weights', w, funcArgs{:});
                    if isprop(mdl, 'KFold')
                        mLoss(n) = mdl.kfoldLoss;
                        cm{n} = confusionmat(mdl.Y, mdl.kfoldPredict);
                    else
                        mLoss(n) = mdl.loss(x(I), y);
                        cm{n} = confusionmat(y, mdl.predict(x(I)));
                    end
                end
                
                s = struct;
                s.nPerm = nPerm;
                s.loss = mLoss;
                s.confMat = cat(3, cm{:});
                s.classNames = mdl.ClassNames;
                null{i} = s;
            end
        end
        
        function [test, null] = BatchFitcecoc(X, y, varargin)
            % Run fitcecoc for each column in X
            % 
            %   [test, null] = BatchFitcecoc(X, y)
            %   [test, null] = BatchFitcecoc(X, y, 'NPerm', 0)
            %   [test, null] = BatchFitcecoc(..., fitcecocArgs)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('NPerm', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Verbose', true, @islogical);
            p.parse(varargin{:});
            nPerm = p.Results.NPerm;
            isVerbose = p.Results.Verbose;
            funcArgs = [fieldnames(p.Unmatched) struct2cell(p.Unmatched)]';
            
            nMdls = size(X, 2);
            test = cell(nMdls, 1);
            null = cell(nMdls, 1);
            for i = 1 : nMdls
                if isVerbose
                    str = sprintf('Fitting models %i/%i\n', i, nMdls);
                    fprintf(str);
                end
                
                % Prepare input for this iteration
                x = X(:,i);
                
                % Remove outliers in each group
                C = unique(y);
                for k = 1 : numel(C)
                    m = y == C(k);
                    xClean = x(m);
                    xClean(isoutlier(xClean)) = NaN;
                    x(m) = xClean;
                end
                
                % Compensate class imbalance with sample weights
                w = ones(size(y));
                for n = 1 : numel(C)
                    m = y==C(n) & ~isnan(x); % use valid samples
                    w(m) = 1 / sum(m);
                end
                
                % Fit test model
                try
                    mdl = fitcecoc(x, y, 'Weights', w, funcArgs{:});
                catch e
                    disp(e);
                    continue
                end
                s = struct;
                if isprop(mdl, 'KFold')
                    s.loss = mdl.kfoldLoss;
                    s.confMat = confusionmat(mdl.Y, mdl.kfoldPredict);
%                     s.AUC = NP.UnitDecode.ComputeAucFromECOC(m, x, y);
                else
                    s.loss = mdl.loss(x, y);
                    s.confMat = confusionmat(y, mdl.predict(x));
%                     s.AUC = []; % to be implemented
                end
                s.classNames = mdl.ClassNames;
                test{i} = s;
                
                % Fit null models with permutation
                if ~nPerm
                    continue
                end
                
                mLoss = zeros(1, nPerm);
                cm = cell(1, nPerm);
                
                parfor n = 1 : nPerm
                    I = randperm(numel(x));
                    mdl = fitcecoc(x(I), y, 'Weights', w, funcArgs{:});
                    if isprop(mdl, 'KFold')
                        mLoss(n) = mdl.kfoldLoss;
                        cm{n} = confusionmat(mdl.Y, mdl.kfoldPredict);
                    else
                        mLoss(n) = mdl.loss(x(I), y);
                        cm{n} = confusionmat(y, mdl.predict(x(I)));
                    end
                end
                
                s = struct;
                s.nPerm = nPerm;
                s.loss = mLoss;
                s.confMat = cat(3, cm{:});
                s.classNames = mdl.ClassNames;
                null{i} = s;
            end
        end
        
        function AUC = ComputeAucFromECOC(mdl, X, Y)
            % Compute AUC of ROC of each binary classifier in a ClassificationPartitionedECOC model
            % 
            %   AUC = ComputeAucFromECOC(mdl, X, Y)
            % 
            
            [numFold, numBinary] = size(mdl.BinaryY);
            
            m = ~isnan(X);
            X = X(m);
            Y = Y(m);
            
            AUC = NaN(numFold, numBinary);
            for i = 1 : numFold
                for j = 1 : numBinary
                    % Get the binary linear classifier
                    bl = mdl.Trained{i}.BinaryLearners{j};
                    if isempty(bl)
%                         disp("ComputeAucFromECOC: empty binary leaner");
                        continue
                    end
                    
                    % Find responses and labels for the partition
                    trainInd = mdl.Partition.training(i);
                    x = X(trainInd);
                    y = Y(trainInd);
                    
                    % Find responses and labels for this binary classifier
                    classList = mdl.ClassNames(mdl.CodingMatrix(:,j) ~= 0);
                    isCla = ismember(y, classList);
                    posCla = mdl.ClassNames(mdl.CodingMatrix(:,j) == 1);
                    x = x(isCla);
                    y = y(isCla);
                    minObs = 3;
                    if sum(y==posCla) < minObs || sum(y~=posCla) < minObs % requires at least 3 trials for each class
                        continue
                    end
                    
%                     % Check if this binary learner is involved in CV
%                     testInd = mdl.Partition.test(i);
%                     yTest = Y(testInd);
%                     if ~any(ismember(yTest, y))
%                         continue
%                     end
                    
                    % Compute AUC
                    [~, scores] = bl.predict(x);
                    scores = scores(:,1);
%                     if numel(unique(y(~isnan(scores)))) < 2
% %                         disp("ComputeAucFromECOC: only one class available for the binary learner");
%                         continue
%                     end
                    [~, ~, ~, auc] = perfcurve(y, scores(:,1), posCla);
                    AUC(i,j) = auc;
                end
            end
        end
        
        
        % Old
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
        
    end
end

