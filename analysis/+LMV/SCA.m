classdef SCA
    
    methods(Static)
        function tbs = LoadResults(filePaths, nCompList)
            % Reorganize SCA output into a table with additional metrics
            % 
            %   tbs = LoadResults(s)
            %   tbs = LoadResults(s, nCompList)
            % 
            
            filePaths = string(filePaths);
            ss = arrayfun(@load, filePaths);
            
            for i = numel(ss) : -1 : 1
                s = ss(i);
                X = s.X;
                nComp = double(s.n_components(:));
                s = rmfield(s, ["X", "n_components"]);
                sComp = structfun(@(x) x, s);
                
                if exist('nCompList', 'var') && ~isempty(nCompList)
                    isComp = ismember(nComp, nCompList);
                    nComp = nComp(isComp);
                    sComp = sComp(isComp);
                end
                
                tb = table;
                tb.nComp = nComp;
                tb.X(:) = {X};
                tb = [tb struct2table(sComp, "AsArray", true)];
                if ~iscell(tb.b_u)
                    tb.b_u = mat2cell(tb.b_u, ones(height(tb),1));
                end
                if ~iscell(tb.b_v)
                    tb.b_v = mat2cell(tb.b_v, ones(height(tb),1));
                end
                
                tbs{i,1} = tb;
            end
        end
        
        function tb = EnrichResultTable(tb)
            % Extract SCA results for a given number of components
            % 
            %   tb = EnrichResultTable(tb)
            % 
            
            % Compute VE from multivariate Gaussian
            maxComp = max(tb.nComp);
            nInputDims = unique(cellfun(@(x) size(x,2), tb.X));
            nullVE = arrayfun(@(x) LMV.SCA.ComputeNullVE(x, maxComp), nInputDims, 'Uni', false);
            
            % Iterate through each model
            for i = 1 : height(tb)
                % Make activation positive
                sn = LMV.SCA.GetLatentDirections(tb.Z{i});
                tb.Z{i} = tb.Z{i}.*sn;
                tb.U{i} = tb.U{i}.*sn;
                tb.b_u{i} = tb.b_u{i}.*sn;
                tb.V{i} = tb.V{i}.*sn';
                
                % Compute metrics
                [tb.r2(i), tb.varExplained{i}, tb.reconLoss(i)] = LMV.SCA.ComputeMetrics(tb.X{i}, tb.U{i}, tb.V{i});
                
                iNull = nInputDims == size(tb.X{i},2);
                tb.relVE{i} = tb.varExplained{i} ./ nullVE{iNull}(1:tb.nComp(i));
            end
        end
        
        function s = GetLatentDirections(varargin)
            % Find whether each latent component is pointing up (1) or down (-1)
            % 
            %   s = GetLatentDirections(Z)
            %   s = GetLatentDirections(X, U, bu)
            % 
            switch numel(varargin)
                case 1
                    Z = varargin{1};
                case 3
                    [X, U, bu] = varargin{:};
                    Z = X*U+bu;
                otherwise
                    error("Wrong number of input arguments. Must be 1 or 3.");
            end
            % [~, I] = max(abs(Z));
            % s = sign(arrayfun(@(a,b) Z(a,b), I, 1:numel(I)));
            s = sign(mean(Z));
        end
        
        function [r2, varExplained, reconLoss] = ComputeMetrics(X, U, V)
            % Compute r-squared and reconstruction loss of a decomposition
            % 
            %   [r2, varExplained, reconLoss] = ComputeMetrics(X, U, V)
            % 
            
            Xhat = X*U*V;
            reconLoss = sum((Xhat - X).^2, "all");
            
            varExplained = MMath.VarExplained(X, U);
            
            n = size(X,1);
            p = size(U,2);
            r2 = 1 - (n-1)/(n-p) * sum((X-Xhat).^2)./sum((X-mean(X)).^2);
            r2 = mean(r2, "omitmissing", "Weights", var(X));
        end
        
        function veNull = ComputeNullVE(nInputDims, nComp)
            % Compute the percent variance explained by PCs of multivariate uniform Gaussian distribution
            % 
            %   veNull = ComputeNullVE(nInputDim, nPC)
            % 
            rng(61);
            nNull = 10;
            for i = nNull : -1 : 1
                X = randn(1e4, nInputDims);
                [~, ~, ~, ~, veNull(i,:)] = pca(X, NumComponents=nComp);
            end
            veNull = veNull(:, 1:nComp);
            veNull = median(veNull, 1);
        end
        
        function tb = QuantifySentenceActivations(ce, Z, CI)
            % Quantify peak, span, and relative variance of the latent dynamics for each sentence in each component
            % 
            %   tb = QuantifySentenceActivations(ce, Z)
            % 
            
            % Find the activated sentences
            % indAct = LMV.SCA.FindActivatedSentences(ce, Z);
            indAct = LMV.SCA.FindActivatedSentences2(Z, CI, ce);
            
            % Reshape latents
            [respTb, tt, tv] = ce.GetTable("resp", "taskTime", "taskValue");
            nComp = size(Z,2);
            Z = reshape(Z, [], height(respTb), nComp);
            
            t = respTb.time{1};
            isPre = t <= tt.cueOn(1);
            isLMV = t > tt.stimOn(1) & t < tt.prodMatchOff(1);
            
            tbs = cell(nComp,1);
            for i = 1 : nComp
                % Make defections pointing up
                z = Z(:,:,i);
                [~, I] = max(abs(z(:)));
                if z(I) < 0
                    z = -z;
                end
                
                % Get features
                zBase = z(isPre,:);
                th = mean(zBase) + std(zBase) * 3;
                isActT = z > th;
                
                v = sum((z-mean(zBase)).^2);
                v = v./sum(v);
                
                tb = table;
                tb.stimId = tv.stimId;
                tb.stimText = tv.stimText;
                tb.peak = max(abs(z))';
                tb.span = mean(isActT(isLMV,:) | isActT(isLMV,:), 1)';
                tb.relVar = v';
                tb.isAct(indAct{i}) = true;
                tb.comp(:) = i;
                
                tbs{i} = tb;
            end
            tb = cat(1, tbs{:});
        end
        
        function tb = QuantifyPhaseActivations(ce, Z)
            % Quantify the variance of the latent dynamics in each task phase of each component
            % 
            %   tb = QuantifyPhaseActivations(ce, Z)
            % 
            
            [respTb, tt] = ce.GetTable("resp", "taskTime");
            nComp = size(Z,2);
            Z = reshape(Z, [], height(respTb), nComp);
            
            t = respTb.time{1};
            % isPre = t <= tt.cueOn(1);
            phases = ["atten", "stim", "delay", "init", "prod"];
            
            tb = table;
            tb.comp = (1:nComp)';
            for i = 1 : nComp
                z = Z(:,:,i);
                % zBase = z(isPre,:);
                % v = (z-mean(zBase)).^2;
                v = z.^2;
                
                for j = 1 : numel(phases)
                    pn = phases(j);
                    evt = tt.(pn)(1);
                    if iscell(evt)
                        evt = evt{1};
                    end
                    m = evt.MaskTimestamps(t);
                    tb.(pn)(i) = mean(v(m,:), "all");
                end
            end
        end
        
        function [indAct, stimText] = FindActivatedSentences(ce, Z)
            % Find activated sentences in each component
            % 
            %   [indAct, stimText] = FindActivatedSentences(ce, Z)
            % 
            
            % Reshape Z to time-by-sen-by-comp
            nComp = size(Z,2);
            Z = reshape(Z, [], ce.numEpochs, nComp);
            
            % Find the activated sentences
            tv = ce.GetTable("taskValue");
            indAct = cell(nComp, 1);
            stimText = cell(nComp, 1);
            
            for k = 1 : nComp
                % Find ranges of outlier values
                z = Z(:,:,k);
                tOut = isoutlier(z', "median", ThresholdFactor=3.5)';
                
                % Ignore periods shorter than 200ms
                for i = 1 : size(tOut,2)
                    bb = MMath.Logical2Bounds(tOut(:,i));
                    for j = 1 : size(bb,1)
                        bInd = bb(j,1) : bb(j,2);
                        tOut(bInd,i) = numel(bInd) >= 20; % 200ms
                    end
                end
                
                indAct{k} = find(any(tOut,1));
                stimText{k} = tv.stimText(indAct{k});
            end
        end
        
        function [indAct, stimText] = FindActivatedSentences2(Z, CI, ce)
            % Find activated sentences in each component
            % 
            %   [indAct, stimText] = FindActivatedSentences(Z, CI, ce)
            % 
            
            % Find locations of significant values
            M = Z<CI(:,:,1) | Z>CI(:,:,2);
            nComp = size(M,2);
            M = reshape(M, [], ce.numEpochs, nComp); % reshape to time-by-sen-by-comp
            Z = reshape(Z, [], ce.numEpochs, nComp);
            
            % Get phase masks
            [tt, tv, respTb] = ce.GetTable("taskTime", "taskValue", "resp");
            stim = tt.stim(1);
            prod = tt.prod{1};
            cue3 = tt.cue3(1);
            t = respTb.time{1};
            isStim = stim.MaskTimestamps(t);
            isDelay = t > stim.T.tmax & t < cue3.T.tOn;
            isInit = t > cue3.T.tOn & t > prod.T.tmin;
            isProd = prod.MaskTimestamps(t);
            
            % Iterate thorugh each component
            indAct = cell(nComp, 1);
            stimText = cell(nComp, 1);
            for k = 1 : nComp
                % Get sig mask
                isSig = M(:,:,k);
                
                % Find ranges of outlier values
                isOut = isoutlier(Z(:,:,k), 2);
                
                % Combine both criteria
                m = isSig & isOut;
                
                % Ignore periods shorter than 200ms
                for i = 1 : size(m,2)
                    bb = MMath.Logical2Bounds(m(:,i));
                    for j = 1 : size(bb,1)
                        bInd = bb(j,1) : bb(j,2);
                        m(bInd,i) = numel(bInd) >= 20; % 200ms
                    end
                end
                
                % Check activation in all phases
                isActPhase = [any(m(isStim,:)); any(m(isDelay,:)); any(m(isInit,:)); any(m(isProd,:))];
                
                indAct{k} =  find(all(isActPhase));
                stimText{k} = tv.stimText(indAct{k});
            end
        end
        
        function PlotLatentNoCI(ce, Z, relVE, tl)
            % 
            
            % Find the activated sentences
            [indAct, stimText] = LMV.SCA.FindActivatedSentences(ce, Z);
            
            % Reshape latents
            respTb = ce.GetTable("resp");
            t = respTb.time{1};
            Z = reshape(Z, [], height(respTb), size(Z,2));
            
            if nargin < 4
                tl = tiledlayout("flow");
                tl.Padding = "compact";
            end
            
            if nargin < 3 || isempty(relVE)
                relVE = NaN(1, height(respTb));
            end
            
            for i = 1 : size(Z,3)
                ax = nexttile(tl);
                
                plot(t, Z(:,:,i)); hold(ax, 'on');
                ax.YLim = [-5 20];
                % ax.YLim(2) = ax.YLim(2) + diff(ax.YLim)*0.2;
                
                colNames = {'cue1', 'stim', 'prod', 'cue3'};
                cc = LMV.Param.GetTaskPhaseColors(colNames);
                NP.TaskBaseClass.PlotEventWindows(ax, ce, 1, t([1 end])', colNames, 'YRange', ax.YLim, 'Colors', cc, 'Alpha', 0.1, 'Text', false);
                
                yText = linspace(ax.YLim(2), ax.YLim(1), 10);
                cc = lines(ce.numEpochs);
                for j = 1 : numel(stimText{i})
                    text(-0.1, yText(j), stimText{i}(j), 'Color', cc(indAct{i}(j),:), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
                end
                
                ax.Title.String = sprintf("# %i (rVE = %.1f)", i, relVE(i));
                ax.XTick = [];
                ax.YTick = [];
                ax.XLabel.String = [];
                ax.YLabel.String = [];
                axis off
            end
        end
        
        function PlotLatent(ce, Z, CI, relVE, tl)
            % 
            % 
            %   PlotLatent(ce, Z, CI, relVE, tl)
            % 
            
            % Find the activated sentences
            [indAct, stimText] = LMV.SCA.FindActivatedSentences2(Z, CI, ce);
            
            % Reshape latents
            respTb = ce.GetTable("resp");
            t = respTb.time{1};
            Z = reshape(Z, [], height(respTb), size(Z,2));
            CI = reshape(CI, [], height(respTb), size(CI,2), 2);
            
            % Get sentence colors
            ccSen = LMV.Param.GetSentenceColors(ce.GetTable("taskValue").stimId);
            
            % Plot each component in a separate panel
            if ~exist('tl', 'var')
                tl = tiledlayout("flow");
                tl.Padding = "compact";
            end
            for i = 1 : size(Z,3)
                ax = nexttile(tl);
                
                % Latents
                z = Z(:,:,i);
                hh = plot(t, z); hold(ax, 'on');
                for n = 1 : numel(hh)
                    hh(n).Color = ccSen(n,:);
                end
                
                % Sentence-null CI
                ci = squeeze(mean(CI(:,:,i,:), 2, "omitmissing"));
                MPlot.ErrorShade(t, mean(ci,2), ci(:,1), ci(:,2), 'IsRelative', false, 'Alpha', 0.2); hold(ax, 'on');
                
                % Phases
                ax.YLim = [-5 20];
                colNames = {'cue1', 'stim', 'prod', 'cue3'};
                cc = LMV.Param.GetTaskPhaseColors(colNames);
                NP.TaskBaseClass.PlotEventWindows(ax, ce, 1, t([1 end])', colNames, 'YRange', ax.YLim, 'Colors', cc, 'Alpha', 0.1, 'Text', false);
                
                % Sentence labels
                yText = linspace(ax.YLim(2), ax.YLim(1), 10);
                for j = 1 : numel(stimText{i})
                    text(-0.1, yText(j), stimText{i}(j), 'Color', ccSen(indAct{i}(j),:), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
                end
                
                if ~exist('relVE', 'var') || isempty(relVE)
                    ax.Title.String = sprintf("# %i", i);
                else
                    ax.Title.String = sprintf("# %i (rVE = %.1f)", i, relVE(i));
                end
                ax.XTick = [];
                ax.YTick = [];
                ax.XLabel.String = [];
                ax.YLabel.String = [];
                axis off
            end
        end
        
        function A = ExtractActivations(ce, Z, phaseName)
            % Extract the magnitude of latent activation in a given phase
            % 
            %   A = ExtractActivation(ce, Z)
            %   A = ExtractActivation(ce, Z, phaseName)
            % 
            % Inputs
            %   Z           A time*sen-by-rep or time*sen-by-comp array.
            % Output
            %   A           A sen-by-rep or sen-by-comp array.
            
            % Reshape latents
            Z = squeeze(Z); % time*sen-by-rep
            respTb = ce.GetTable("resp");
            Z = reshape(Z, [], height(respTb), size(Z,2)); % time-by-sen-by-rep
            
            % Find phase mask
            t = respTb.time{1};
            if exist("phaseName", "var") && ~isempty(phaseName)
                evt = ce.GetTable("taskTime").(phaseName)(1);
                isPhase = evt.MaskTimestamps(t);
            else
                isPhase = true(size(t));
            end
            
            % Compute mean phase activation for each repeat
            for i = size(Z,3) : -1 : 1
                A(:,i) = mean(Z(isPhase,:,i), "omitmissing"); % sen-by-rep
            end
        end
        
        function [A1, A2] = MatchActivations(Z1, Z2)
            % Extract the magnitude of latent activation in a given phase
            % 
            %   [A1, A2] = MatchActivations(Z1, Z2)
            % 
            % Inputs
            %   Z1, Z2      Two time*sen-by-comp arrays of the same size.
            % Output
            %   A1, A2      Two sen-by-comp arrays.
            
            % Match components using correlation
            M = corr(Z1, Z2);
            [~, I] = max(M, [], 2);
            Z2 = Z2(:,I);
            
            % Match sentences
            nComp = size(Z1, 2);
            nSen = 14;
            Z1 = reshape(Z1, [], nSen, nComp); % time-by-sen-by-comp
            Z2 = reshape(Z2, [], nSen, nComp); % time-by-sen-by-comp
            for c = 1 : nComp
                z1 = Z1(:,:,c);
                z2 = Z2(:,:,c);
                M = corr(z1, z2);
                [~, I] = max(M, [], 2);
                
            end
            
            % Find phase mask
            t = respTb.time{1};
            evt = ce.GetTable("taskTime").(phaseName)(1);
            isPhase = evt.MaskTimestamps(t);
            
            % Compute mean phase activation for each repeat
            for i = size(Z,3) : -1 : 1
                A(:,i) = mean(Z(isPhase,:,i), "omitmissing"); % sen-by-rep
            end
        end
        
        
        
        % Sentence decoding
        function s = ClassifySentences(ce, Z)
            % 
            
            tv = ce.GetTable("taskValue");
            nTime = size(Z,1) / height(tv);
            y = repelem(tv.stimId, nTime);
            
            mdl = fitcecoc(Z, y, 'KFold', 10);
            yHat = kfoldPredict(mdl);
            isSame = strcmp(y, yHat);
            r = mean(reshape(isSame, nTime, []), 2);
            
            s = struct;
            s.r = r;
            
            nShift = 0;
            if ~nShift
                return
            end
            rNull = cell(nShift,1);
            for n = 1 : nShift
                mdl = fitcecoc(Z, y, 'KFold', 10);
                yHat = kfoldPredict(mdl);
                isSame = strcmp(y, yHat);
                rNull{n} = mean(reshape(isSame, nTime, []), 2);
            end
            rNull = cat(2, rNull{:});
            rPval = MMath.EstimatePval(r, rNull);
            [mNull, ~, ~, ciNull] = MMath.MeanStats(rNull, 0.05);
            
            s.rPval = rPval;
            s.rNullMean = mNull;
            s.rNullCI = ciNull;
        end
        
    end
end