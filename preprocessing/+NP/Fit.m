classdef Fit
    methods(Static)
        % Model fitting
        function mdlData = BatchLasso(X, Y, varargin)
            % Run lasso for each response in Y
            % 
            %   mdlData = BatchLasso(X, Y)
            %   mdlData = BatchLasso(..., 'Lambda', [])
            %   mdlData = BatchLasso(..., 'NumShifts', 0)
            %   mdlData = BatchLasso(..., 'MinShift', 0)
            %   mdlData = BatchLasso(..., 'Verbose', true)
            %   mdlData = BatchLasso(..., lassoArgs)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('Lambda', [], @(x) isnumeric(x) || iscell(x));
            p.addParameter('NumShifts', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('MinShift', 1, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Verbose', true, @islogical);
            p.parse(varargin{:});
            L = p.Results.Lambda;
            nShifts = p.Results.NumShifts;
            minShift = p.Results.MinShift;
            isVerbose = p.Results.Verbose;
            lassoArgs = [fieldnames(p.Unmatched) struct2cell(p.Unmatched)]';
            
            % Expand lambda
            if ~iscell(L)
                L = {L};
            end
            if isscalar(L)
                L = repelem(L, size(Y,2));
            end
            
            % Iterate through responses
            nMdls = size(Y, 2);
            mdlData = cell(nMdls, 1);
            for i = 1 : numel(mdlData)
                if isVerbose
                    str = sprintf('Fitting models %i/%i\n', i, nMdls);
                    fprintf(str);
                end
                
                % Fit model
                y = Y(:,i);
                [B, m] = lasso(X, y, 'Lambda', L{i}, 'CV', 10, lassoArgs{:});
                
                mse0 = mean((y-mean(y,'omitnan')).^2, 'omitnan');
                m.MSE0 = mse0;
                
                md = struct;
                md.mdl = m;
                md.Beta = B;
                md.Bias = m.Intercept;
                md.r2 = 1 - m.MSE / mse0;
                
                % Fit null models with random time shifts
                if nShifts
                    if size(B,2) > 1
                        error("Bootstrap is not supported when using multiple Lambda values for each response.");
                    end
                    
                    B = zeros(size(B,1), nShifts);
                    C = zeros(1, nShifts);
                    r2 = zeros(1, nShifts);
                    
                    parfor n = 1 : nShifts
                        ys = circshift(y, randi([0 numel(y)]+[1 -1]*minShift));
                        [B(:,n), m] = lasso(X, ys, 'Lambda', L{i}, 'CV', 10, lassoArgs{:});
                        C(n) = m.Intercept;
                        r2(n) = 1 - m.MSE / mse0;
                    end
                    
                    % Estimate pval
                    pvalMethod = 'gumbel';
                    bP = MMath.EstimatePval(md.Beta, B', 'Method', pvalMethod)';
                    cP = MMath.EstimatePval(md.Bias, C, 'Method', pvalMethod);
                    rP = MMath.EstimatePval(md.r2, r2, 'Tail', 'right', 'Method', pvalMethod);
                    
                    nd = struct;
                    nd.nShifts = nShifts;
                    nd.r2 = r2;
                    nd.BetaPval = bP;
                    nd.BiasPval = cP;
                    nd.r2Pval = rP;
                    nd.pvalMethod = pvalMethod;
                    md.null = nd;
                end
                
                mdlData{i} = md;
            end
        end
        
        function mdlData = BatchSmoothLinear(X, Y, varargin)
            % Fit smooth linear model for each response in Y
            % 
            %   mdlData = BatchSmoothLinear(X, Y)
            %   mdlData = BatchSmoothLinear(..., 'Lambda', [])
            %   mdlData = BatchSmoothLinear(..., 'NumShifts', 0)
            %   mdlData = BatchSmoothLinear(..., 'MinShift', 0)
            %   mdlData = BatchSmoothLinear(..., 'Verbose', true)
            %   mdlData = BatchSmoothLinear(..., lassoArgs)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('Lambda', [], @(x) isnumeric(x) || iscell(x));
            p.addParameter('NumShifts', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('MinShift', 1, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Verbose', true, @islogical);
            p.parse(varargin{:});
            L = p.Results.Lambda;
            nShifts = p.Results.NumShifts;
            minShift = p.Results.MinShift;
            isVerbose = p.Results.Verbose;
            lassoArgs = [fieldnames(p.Unmatched) struct2cell(p.Unmatched)]';
            
            % Expand lambda
            if ~iscell(L)
                L = {L};
            end
            if isscalar(L)
                L = repelem(L, size(Y,2));
            end
            
            % Iterate through responses
            nMdls = size(Y, 2);
            mdlData = cell(nMdls, 1);
            for i = 1 : numel(mdlData)
                if isVerbose
                    str = sprintf('Fitting models %i/%i\n', i, nMdls);
                    fprintf(str);
                end
                
                % Fit model
                y = Y(:,i);
                [B, m] = NP.Fit.FitSmoothLinear(X, y, 'Lambda', L{i}, 'CV', 10, lassoArgs{:});
                
                mse0 = mean((y-mean(y,'omitnan')).^2, 'omitnan');
                m.MSE0 = mse0;
                
                md = struct;
                md.mdl = m;
                md.Beta = B;
                md.Bias = m.Intercept;
                md.r2 = 1 - m.MSE / mse0;
                
                % Fit null models with random time shifts
                if nShifts
                    if size(B,2) > 1
                        error("Bootstrap is not supported when using multiple Lambda values for each response.");
                    end
                    
                    B = zeros(size(B,1), nShifts);
                    C = zeros(1, nShifts);
                    r2 = zeros(1, nShifts);
                    
                    parfor n = 1 : nShifts
                        ys = circshift(y, randi([0 numel(y)]+[1 -1]*minShift));
                        [B(:,n), m] = NP.Fit.FitSmoothLinear(X, ys, 'Lambda', L{i}, 'CV', 10, lassoArgs{:});
                        C(n) = m.Intercept;
                        r2(n) = 1 - m.MSE / mse0;
                    end
                    
                    % Estimate pval
                    pvalMethod = 'gumbel';
                    bP = MMath.EstimatePval(md.Beta, B', 'Method', pvalMethod)';
                    cP = MMath.EstimatePval(md.Bias, C, 'Method', pvalMethod);
                    rP = MMath.EstimatePval(md.r2, r2, 'Tail', 'right', 'Method', pvalMethod);
                    
                    nd = struct;
                    nd.nShifts = nShifts;
                    nd.r2 = r2;
                    nd.BetaPval = bP;
                    nd.BiasPval = cP;
                    nd.r2Pval = rP;
                    nd.pvalMethod = pvalMethod;
                    md.null = nd;
                end
                
                mdlData{i} = md;
            end
        end
        
        function [B, m] = FitSmoothLinear(X, y, varargin)
            % Fit linear model with regularization for smooth kernel
            % 
            %   [B, m] = FitSmoothLinear(X, y)
            %   [B, m] = FitSmoothLinear(X, y, ..., 'Lambda', [])
            %   [B, m] = FitSmoothLinear(X, y, ..., 'CV', 0)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('Lambda', [], @isnumeric);
            p.addParameter('CV', 0, @(x) isnumeric(x) && isscalar(x));
            p.parse(varargin{:});
            L = p.Results.Lambda;
            nCV = p.Results.CV;
            
            if isempty(L)
                L = logspace(-2, 2, 30);
            end
            
            % Construct penalty matrix
            [nSp, nW] = size(X);
            D = eye(nW) - diag(ones(nW-1,1), 1);
            D(end) = 0;
            
            % Custom CV partitions
            if nCV > 0
                c0 = cvpartition(nSp, "KFold", nCV);
                testInd = test(c0, "all");
                testInd = [testInd; false(nW, nCV)];
            end
            
            nL = numel(L);
            B = NaN(nW, nL);
            lamMdls = cell(nL, 1);
            for i = 1 : nL
                Xr = [X; D*sqrt(L(i)*nSp)];
                yr = [y; zeros(nW,1)];
                
                % Fit full model to get coefficients
                [B(:,i), lamMdls{i}] = lasso(Xr, yr, 'Lambda', 0);
                
                % Compute crossvalidated MSE
                if ~nCV
                    continue
                end
                cvMSE = NaN(nCV,1);
                for j = 1 : nCV
                    m = ~testInd(:,j);
                    [b, info] = lasso(Xr(m,:), yr(m), 'Lambda', 0);
                    mse = mean((Xr(~m,:)*b+info.Intercept - yr(~m)).^2);
                    cvMSE(j) = mse;
                end
                [lamMdls{i}.MSE, lamMdls{i}.SE] = MMath.MeanStats(cvMSE);
            end
            
            m = lamMdls{1};
            m.Lambda = L(:)';
            fn = ["Intercept", "DF", "MSE"];
            for i = 1 : numel(fn)
                m.(fn(i)) = cellfun(@(x) x.(fn(i)), lamMdls);
            end
            
            if ~nCV
                return
            end
            
            m.SE = cellfun(@(x) x.SE, lamMdls);
            
            [MinMSE, m.IndexMinMSE] = min(m.MSE);
            m.LambdaMinMSE = L(m.IndexMinMSE);
            
            is1SE = mse < MinMSE + m.SE(m.IndexMinMSE);
            L(~is1SE) = -Inf;
            [~, m.Index1SE] = max(L);
            m.Lambda1SE = L(m.Index1SE);
            
        end
        
        function mdlData = BatchLsqnonneg(X, Y, varargin)
            % Run nonnegative linear regression for each response in Y
            % 
            %   mdlData = BatchLsqnonneg(X, Y)
            %   mdlData = BatchLsqnonneg(..., 'PredictorNames', "p"+1:size(X,2))
            %   mdlData = BatchLsqnonneg(..., 'NumShifts', 0)
            %   mdlData = BatchLsqnonneg(..., 'MinShift', 0)
            %   mdlData = BatchLsqnonneg(..., 'Verbose', true)
            %   mdlData = BatchLsqnonneg(..., lsqnonnegArgs)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('PredictorNames', "p"+(1:size(X,2)), @(x) iscellstr(x) || isstring(x));
            p.addParameter('NumShifts', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('MinShift', 1, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Verbose', true, @islogical);
            p.parse(varargin{:});
            pn = p.Results.PredictorNames;
            nShifts = p.Results.NumShifts;
            minShift = p.Results.MinShift;
            isVerbose = p.Results.Verbose;
            funcArgs = [fieldnames(p.Unmatched) struct2cell(p.Unmatched)]';
            
            % Iterate through responses
            nMdls = size(Y, 2);
            mdlData = cell(nMdls, 1);
            for i = 1 : numel(mdlData)
                if isVerbose
                    str = sprintf('Fitting models %i/%i\n', i, nMdls);
                    fprintf(str);
                end
                
                % Fit model
                y = Y(:,i);
                [B, SSE] = lsqnonneg(X, y, funcArgs{:});
                md = struct;
                md.mdl.PredictorNames = pn;
                md.Beta = B;
                md.Bias = 0;
                SSE0 = sum((y-mean(y,'omitnan')).^2, 'omitnan');
                md.r2 = 1 - SSE / SSE0;
                
                % Fit null models with random time shifts
                if nShifts
                    
                    B = zeros(size(B,1), nShifts);
                    C = zeros(1, nShifts);
                    r2 = zeros(1, nShifts);
                    
                    parfor n = 1 : nShifts
                        ys = circshift(y, randi([0 numel(y)]+[1 -1]*minShift));
                        [B(:,n), SSE] = lsqnonneg(X, ys, funcArgs{:});
                        r2(n) = 1 - SSE / SSE0;
                    end
                    
                    nd = struct;
                    nd.nShifts = nShifts;
                    nd.Beta = B;
                    nd.Bias = C;
                    nd.r2 = r2;
                    md.null = nd;
                end
                
                mdlData{i} = md;
            end
        end
        
        function mdlData = BatchFitrlinear(X, Y, varargin)
            % Run lasso for each response in Y
            % 
            %   mdlData = BatchFitrlinear(X, Y)
            %   mdlData = BatchFitrlinear(..., 'Verbose', true)
            %   mdlData = BatchFitrlinear(..., 'Lambda', 'auto')
            %   mdlData = BatchFitrlinear(..., 'NumShifts', 0)
            %   mdlData = BatchFitrlinear(..., 'MinShift', 0)
            %   mdlData = BatchFitrlinear(..., fitrlinearArgs)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('Verbose', true, @islogical);
            p.addParameter('Lambda', 'auto', @(x) isnumeric(x) || iscell(x));
            p.addParameter('NumShifts', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('MinShift', 0, @(x) isnumeric(x) && isscalar(x));
            p.parse(varargin{:});
            L = p.Results.Lambda;
            isVerbose = p.Results.Verbose;
            nShifts = p.Results.NumShifts;
            minShift = p.Results.MinShift;
            funcArgs = [fieldnames(p.Unmatched) struct2cell(p.Unmatched)]';
            
            % Expand lambda
            if ~iscell(L)
                L = {L};
            end
            if isscalar(L)
                L = repelem(L, size(Y,2));
            end
            
            nMdls = size(Y, 2);
            mdlData = cell(nMdls, 1);
            parfor i = 1 : nMdls
                if isVerbose
                    str = sprintf('Fitting models %i/%i\n', i, nMdls);
                    fprintf(str);
                end
                
                % Fit model
                y = Y(:,i);
                m = fitrlinear(X, y, 'Lambda', L{i}, funcArgs{:});
                md = struct;
                md.mdl = m;
                md.Beta = m.Beta;
                md.Bias = m.Bias;
                mse0 = mean((y-mean(y,'omitnan')).^2, 'omitnan');
                md.r2 = 1 - m.loss(X, y) / mse0;
                
                % Fit null models with random time shifts
                if nShifts
                    B = zeros(numel(md.Beta), nShifts);
                    C = zeros(1, nShifts);
                    r2 = zeros(1, nShifts);
                    
                    for n = 1 : nShifts
                        d = randi([0 numel(y)]+[1 -1]*minShift);
                        ys = circshift(y, d);
                        m = fitrlinear(X, ys, funcArgs{:});
                        B(:,n) = m.Beta;
                        C(n) = m.Bias;
                        r2(n) = 1 - m.loss(X, ys) / mse0;
                    end
                    
                    nd = struct;
                    nd.nShifts = nShifts;
                    nd.Beta = B;
                    nd.Bias = C;
                    nd.r2 = r2;
                    md.null = nd;
                end
                
                mdlData{i} = md;
            end
        end
        
        function [mdls, null] = BatchFitcecoc(X, Y, varargin)
            % Run fitcecoc for each column in Y
            % 
            %   [mdls, null] = BatchFitcecoc(X, Y)
            %   [mdls, null] = BatchFitcecoc(..., 'NumShifts', 0)
            %   [mdls, null] = BatchFitcecoc(..., 'MinShift', 0)
            %   [mdls, null] = BatchFitcecoc(..., fitcecocArgs)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParameter('NumShifts', 0, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('MinShift', 0, @(x) isnumeric(x) && isscalar(x));
            p.parse(varargin{:});
            nShifts = p.Results.NumShifts;
            minShift = p.Results.MinShift;
            funcArgs = [fieldnames(p.Unmatched) struct2cell(p.Unmatched)]';
            
            nMdls = size(Y, 2);
            mdls = cell(nMdls, 1);
            null = cell(nMdls, 1);
            for i = 1 : nMdls
                % Fit model
                y = Y(:,i);
                
                w = ones(size(y));
                C = unique(y);
                for n = 1 : numel(C)
                    m = y==C(n);
                    w(m) = 1 / sum(m);
                end
                
                mdls{i} = fitcecoc(X, y, 'Weights', w, funcArgs{:});
                
                % Fit null models with random time shifts
                if ~nShifts
                    continue
                end
                
                medLoss = zeros(1, nShifts);
                cm = cell(1, nShifts);
                B = cell(1, nShifts);
                C = cell(1, nShifts);
                
                parfor n = 1 : nShifts
                    d = randi([0 numel(y)]+[1 -1]*minShift);
                    ys = circshift(y, d);
                    ws = circshift(w, d);
                    m = fitcecoc(X, ys, 'Weights', ws, funcArgs{:});
                    
                    losses = m.kfoldLoss('Mode', 'individual');
                    medLoss(n) = median(losses);
                    [~, iMed] = min(abs(medLoss(n) - losses));
                    
                    bb = cellfun(@(x) x.Beta, m.Trained{iMed}.BinaryLearners, 'Uni', false);
                    B{n} = cat(2, bb{:});
                    C{n} = cellfun(@(x) x.Bias, m.Trained{iMed}.BinaryLearners);
                    
                    cm{n} = confusionmat(m.Y, m.kfoldPredict);
                end
                
                s.nShifts = nShifts;
                s.kfoldLoss = medLoss;
                s.confMat = cat(3, cm{:});
                % s.Beta = cat(3, B{:});
                % s.Bias = cat(2, C{:});
                
                null{i} = s;
            end
        end
        
        % Utilities
        function cTbs = LoadModels(mdlNames, varargin)
            % Load model files with tables of the same recording merged into one
            % 
            %   uTbs = LoadModels(mdlNames)
            %   uTbs = LoadModels(mdlNames, mdlDir)
            %   uTbs = LoadModels(mdlNames, mdlDir, recIds)
            %   uTbs = LoadModels(mdlNames, ..., 'clusTbs', []);
            % 
            % Inputs
            %   mdlNames        A list of model names each in the format "<feature_set>_<target>".
            %   mdlDir          The folder where models are saved.
            %   recIds          Include only recordings based on this list of recording IDs.
            %   'clusTbs'       A cell array of clusTb to load models into.
            % Output
            %   uTbs            A cell array of clusTb with models added to the last columns.
            % 
            
            p = inputParser;
            p.addOptional('mdlDir', [], @(x) isempty(x) || exist(x, 'dir'));
            p.addOptional('recIds', [], @(x) ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('clusTbs', [], @(x) iscell(x) || istable(x));
            p.parse(varargin{:});
            mdlDir = p.Results.mdlDir;
            recIds = string(p.Results.recIds);
            cTbs = p.Results.clusTbs;
            
            mdlNames = string(mdlNames);
            
            if isempty(mdlDir)
                switch mdlNames(1)
                    case "phase"
                        mdlDir = fullfile(NP.Data.GetAnalysisRoot, 'phase_resp', 'nlr', 'mdls');
                    case "sen4"
                        mdlDir = fullfile(NP.Data.GetAnalysisRoot, 'decode', 'su_sentence', 'sen4_mdls');
                    otherwise
                        error("User needs to provide the model directory.");
                end
            end
            
            if isempty(recIds)
                if isempty(mdlNames)
                    mSearch = MBrowse.Dir2Table(fullfile(mdlDir, "*_clusTb.mat"));
                else
                    mSearch = MBrowse.Dir2Table(fullfile(mdlDir, "*_clusTb_" + mdlNames(1) + ".mat"));
                end
                recIds = unique(cellfun(@(x) string(NP.SE.GetID(x)), mSearch.name));
            end
            recIds = string(recIds);
            
            if isempty(cTbs)
                cTbs = repmat({table}, size(recIds));
            end
            if ~iscell(cTbs)
                cTbs = {cTbs};
            end
            
            % Load model tables
            for i = 1 : numel(recIds)
                cTb = cTbs{i};
                
                if isempty(mdlNames)
                    load(fullfile(mdlDir, recIds(i)+"_clusTb.mat"), 'clusTb');
                    m = ~ismember(clusTb.Properties.VariableNames, cTb.Properties.VariableNames);
                    cTb = [cTb clusTb(:,m)];
                else
                    for k = 1 : numel(mdlNames)
                        load(fullfile(mdlDir, recIds(i)+"_clusTb_"+mdlNames(k)+".mat"), 'clusTb');
                        m = ~ismember(clusTb.Properties.VariableNames, cTb.Properties.VariableNames);
                        cTb = [cTb clusTb(:,m)];
                    end
                end
                
                cTbs{i} = cTb;
            end
        end
        
        function [lambda, err] = GetOptimalLambda(mdls, varargin)
            % Get optimal lambda for lasso
            % 
            %   [lambda, err] = GetOptimalLambda(mdls)
            %   [lambda, err] = GetOptimalLambda(..., 'ConvexOnly', true)
            % 
            % Input
            %   mdls            An array of model struct. Each struct contains a series of Lambda values. 
            %   'ConvexOnly'    Whether or not the MSE from the tested series of lambda values needs to 
            %                   be convex (i.e. minimum not at the two extremes) or not. If true, the 
            %                   returned optimal lambda of a model with non-convex MSE will be NaN. 
            % 
            
            p = inputParser;
            p.addParameter('ConvexOnly', true, @islogical);
            p.parse(varargin{:});
            isConv = p.Results.ConvexOnly;
            
            if ~iscell(mdls)
                mdls = num2cell(mdls);
            end
            
            lambda = NaN(size(mdls));
            err = NaN(size(mdls));
            for i = 1 : numel(mdls)
                if isempty(mdls{i})
                    continue
                end
                m = mdls{i}.mdl;
                if isfield(m, 'MSE')
                    % For models fitted by lasso
                    [err(i), I] = min(m.MSE);
                    if isConv && ismember(I, [1 numel(m.MSE)])
                        lambda(i) = NaN;
                    else
                        lambda(i) = m.Lambda(I);
                    end
                elseif isa(m, 'RegressionLinear')
                    % For models fitted by fitrlinear
                    lambda(i) = m.Lambda;
                end
            end
        end
        
        % Use only in a8_nlr.m
        function tb = GetPredictorTable(mdls, varargin)
            % Extract weights or other predictor values and organize them in a table
            % 
            %   tb = GetPredictorTable(mdls)
            %   tb = GetPredictorTable(mdls, func)
            % 
            % Inputs
            %   mdls        A cell array of structs. Each element contains a struct of model output.
            %   func        A handle to the function that extract specific values from a model struct.
            %               The default is @(x)x.Beta, which extracts the predictor weights.
            % Output
            %   tb          A table containing weights or other predictor values for each predictor.
            % 
            % See also MMath.EstimatePval
            
            p = inputParser;
            p.addOptional('func', @(x) x.Beta, @(x) isa(x, 'function_handle'));
            p.parse(varargin{:});
            func = p.Results.func;
            
            % Standardize input as cell array
            if ~iscell(mdls)
                mdls = num2cell(mdls);
            end
            
            % Preallocation
            k = find(~cellfun(@isempty, mdls), 1);
            if isempty(k)
                warning("All input models are empty. Function returns [].");
                tb = [];
                return
            end
            val = NaN(numel(mdls), numel(func(mdls{k})));
            
            % Extract predictor values
            for i = numel(mdls) : -1 : 1
                if isempty(mdls{i})
                    continue
                end
                val(i,:) = func(mdls{i});
            end
            
            % Make table
            mdl1 = mdls{k};
            if isempty(mdl1.mdl.PredictorNames)
                predNames = "p" + (1 : numel(mdl1.Beta));
            else
                predNames = string(mdl1.mdl.PredictorNames);
            end
            tb = array2table(val, 'VariableNames', predNames);
        end
        
        function [pTb, sTb] = GetSigTable(mdls, varargin)
            % Get p-value and the level of significance for each predictor weight in a table
            % 
            %   [pTb, sTb] = GetSigTable(mdls)
            %   [pTb, sTb] = GetSigTable(mdls, 'Tail', 'two')
            %   [pTb, sTb] = GetSigTable(mdls, 'AlphaList', [0.05 0.01 0.001])
            % 
            % Inputs
            %   mdls            A cell array of structs. Each element contains a struct of model output.
            %   'Tail'          See MMath.EstimatePval
            %   'AlphaList'     See MMath.EstimatePval
            % Outputs
            %   pTb             A table containing p-values for each predictor.
            %   sTb             A table containing significance levels for each predictor.
            % 
            % See also MMath.EstimatePval
            
            p = inputParser;
            p.addParameter('Tail', 'two');
            p.addParameter('AlphaList', [0.05 0.01 0.001], @(x) isnumeric(x));
            p.parse(varargin{:});
            tailMode = p.Results.Tail;
            alphaList = p.Results.AlphaList;
            
            % Compute p-values
            pTb = NP.Fit.GetPredictorTable(mdls, @(x) MMath.EstimatePval(x.Beta, x.null.Beta', 'Tail', tailMode));
            
            % Get significance levels
            alphaList = permute(alphaList(:), [2 3 1]);
            L = sum(pTb{:,:} < alphaList, 3);
            sTb = pTb;
            sTb{:,:} = L;
        end
        
        function tb = CatSameVariables(tbs)
            % Horizontally concatenate columns of the same variables across multiple pval tables
            % 
            %   tb = CatSameVariables(tbs)
            % 
            
            % Preallocation
            k = find(~cellfun(@isempty, tbs), 1);
            tb = tbs{k};
            for i = 1 : width(tb)
                tb.(i) = NaN(height(tb), numel(tbs));
            end
            
            % Fill in values
            for i = 1 : numel(tbs)
                if isempty(tbs{i}) % certain sentences may not be avialble in some recordings
                    continue
                end
                for k = 1 : width(tb)
                    tb.(k)(:,i) = tbs{i}.(k);
                end
            end
        end
        
        % Not in use
        function [lvl, lower, upper] = GetSigLevel(X, XBoot, alphaList, tailMode)
            % Get significance level of each predictor weight
            % 
            %   [lvl, lower, upper] = GetSigLevel(X, XBoot)
            %   [lvl, lower, upper] = GetSigLevel(X, XBoot, alphaList)
            %   [lvl, lower, upper] = GetSigLevel(X, XBoot, alphaList, tailMode)
            % 
            
            if ~exist('tailMode', 'var') || isempty(tailMode)
                tailMode = 'both';
            end
            
            if ~exist('alphaList', 'var') || isempty(alphaList)
                alphaList = [0.05 0.01 0.001];
            end
            p = alphaList(:)' * 100;
            
            switch tailMode
                case 'both'
                    lower = prctile(XBoot, p/2, 2);
                    upper = prctile(XBoot, 100-p/2, 2);
                case 'left'
                    lower = prctile(XBoot, p, 2);
                    upper = Inf(size(lower));
                case 'right'
                    upper = prctile(XBoot, 100-p, 2);
                    lower = -Inf(size(upper));
                otherwise
                    error("tailMode must be 'both', 'left', or 'right', but was '%s'.", tailMode);
            end
            
            sig = X < lower | X > upper;
            lvl = sum(sig, 2);
            lower(isinf(lower)) = 0;
            upper(isinf(upper)) = 0;
        end
        
    end
end
