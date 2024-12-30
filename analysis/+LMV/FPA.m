classdef FPA
    
    methods(Static)
        function s = ComputeRelativeDist(vTest, vRef, varargin)
            % Compute relative distribution and test for statistical significance
            % 
            %   s = ComputeRelativeDepth(vTest, vRef)
            % 
            % 
            
            p = inputParser;
            p.addParameter('Name', "depth", @(x) isstring(x) || ischar(x));
            p.parse(varargin{:});
            dName = p.Results.Name;
            
            % Test difference
            s = struct;
            s.vTest = vTest;
            s.vRef = vRef;
            [~, s.pval, s.ks2stat] = kstest2(vRef, vTest);
            
            % Fit distribution
            switch dName
                case 'depth'
                    x = linspace(0, 6000, 100);
                    % d = d(d>x(1) & d<x(end));
                    dArgs = {'Kernel', 'Kernel', 'normal', 'Bandwidth', 500};
                otherwise
                    error("'%s' is not a valid distribution name.", dName);
            end
            pdRef = fitdist(vRef, dArgs{:});
            pdTest = fitdist(vTest, dArgs{:});
            yRef = pdf(pdRef, x);
            yTest = pdf(pdTest, x);
            
            % Compute relative distribution
            yRef = yRef/sum(yRef);
            yTest = yTest/sum(yTest);
            yRatio = yTest./yRef;
            
            s.x = x;
            s.yRef = yRef;
            s.yTest = yTest;
            s.yRatio = yRatio;
        end
        
        function [x, P] = BootDepthDist(d, varargin)
            % Bootstrap fitting of probability distribution across depth
            % 
            %   [x, P] = BootDepthDist(d, varargin)
            % 
            
            p = inputParser;
            p.addParameter("NBoot", 200, @(x) isnumeric(x) && isscalar(x));
            p.parse(varargin{:});
            nBoot = p.Results.NBoot;
            
            x = linspace(0, 6000, 100);
            % d = d(d>x(1) & d<x(end));
            dArgs = {'Kernel', 'Kernel', 'normal', 'Bandwidth', 500};
            
            P = zeros(numel(x), nBoot);
            for i = 1 : nBoot
                % Random sampling
                k = round(numel(d)*0.5);
                dRs = randsample(d, k);
                
                % Fit distribution
                pd = fitdist(dRs, dArgs{:});
                P(:,i) = pdf(pd, x);
            end
        end
        
        function s = ComputeDist(uTb, xVar, yVar, varargin)
            % Compute conditional probability distribution between two unit attributes
            % 
            %   ComputeDist(uTb, xVar, yVar)
            %   ComputeDist(uTb, xVar, yVar, 'NBoot', 1000)
            % 
            
            p = inputParser;
            p.addParameter('NBoot', 1000, @(x) isscalar(x) && isnumeric(x));
            p.parse(varargin{:});
            nBoot = p.Results.NBoot;
            
            % Process variables
            varSets = {xVar, yVar};
            for i = numel(varSets) : -1 : 1
                vn = varSets{i};
                val = uTb{:,vn};
                
                if iscolumn(val)
                    % Convert single variable form to logical array form
                    if vn == "depth"
                        % Custom bin edges for depth
                        edges = [0:.5:5 Inf] * 1e3;
                        centers = MMath.BinEdges2Centers(edges(1:end-1));
                        ticks = interp1(centers, 1:numel(centers), (0:5)*1e3, 'linear', 'extrap');
                        tickLabels = 0:5;
                        label = 'Depth (mm)';
                    elseif iscategorical(val) || isstring(val) || iscellstr(val) || (numel(unique(val)) < 30)
                        % Convert possible categorical variable to indices, e.g. {'nmfcId', 'nlrId', 'wfId', 'recId'}
                        [val, tickLabels] = findgroups(val);
                        ticks = unique(val);
                        edges = MMath.BinCenters2Edges(ticks);
                        label = '';
                    else
                        % Treat as continuous variable
                        warning("'%s' is not a supported variable", vn);
                        [~, edges, I] = histcounts(val);
                        ticks = unique(I);
                        tickLabels = ticks;
                        label = vn;
                    end
                    [~, ~, I] = histcounts(val, edges);
                    X = false(numel(val), numel(edges)-1);
                    for j = 1 : size(X,2)
                        X(I==j,j) = true;
                    end
                else
                    % Input is logical array form
                    X = val;
                    ticks = 1 : numel(vn);
                    tickLabels = vn;
                    label = [];
                end
                
                sVar(i).X = X;
                sVar(i).ticks = ticks;
                sVar(i).tickLabels = tickLabels;
                sVar(i).label = label;
            end
            
            % Compute conditional probabilities
            X1 = sVar(1).X;
            X2 = sVar(2).X;
            P = computeCondProb(X1, X2);
            Pnull = zeros([size(P) nBoot]);
            for i = 1 : nBoot
                I = randperm(size(X1,1));
                Pnull(:,:,i) = computeCondProb(X1(I,:), X2);
            end
            
            function P = computeCondProb(X1, X2)
                P = zeros(size(X2,2), size(X1,2));
                for k = 1 : size(P,1)
                    m = X2(:,k);
                    P(k,:) = sum(X1(m,:), 1);
                end
                P = P ./ sum(P,2);
                P(isnan(P)) = 0;
            end
            
            % Test statistical significance
            pval = NaN(size(P));
            for i = 1 : size(P,1)
                disp(i);
                pval(i,:) = MMath.EstimatePval(P(i,:), permute(Pnull(i,:,:), [3 2 1]), 'Tail', 'right');
            end
            
            s.sVar = sVar;
            s.P = P;
            s.Pnull = Pnull;
            s.pval = pval;
        end
        
        function PlotDist(s, varargin)
            % Plot a heatmap of probability distribution between two unit attributes
            % 
            %   PlotDist(s)
            %   PlotDist(s, 'Style', "heatmap")
            % 
            
            p = inputParser;
            p.addParameter('Style', "heatmap", @(x) any(x == ["heatmap", "trace"]));
            p.parse(varargin{:});
            style = p.Results.Style;
            
            % Unpack input struct
            P = s.P; % - mean(s.Pnull, 3);
            pval = s.pval * numel(P);
            sVar = s.sVar;
            
            % Get coordinates of sig stars
            [x, y] = ind2sub(size(pval), find(pval<0.05));
            
            % Make plots
            ax = gca;
            if style == "heatmap"
                imagesc(ax, P);
                hold(ax, 'on');
                if ~isempty(x)
                    plot(ax, x, y, '*', 'Color', 'r');
                end
%                 cmap = colormap(ax, 'gray');
%                 colormap(ax, flip(cmap));
                colormap(ax, 'gray');
                colorbar(ax);
                ax.CLim(1) = 0;
            elseif style == "trace"
                x = 1 : size(P,2);
                y = (1 : size(P,1));
                MPlot.PlotTraceLadder(x', P'*2, flip(y'), 'Parent', ax);
                ax.XLim = [-.5 .5] + x([1 end]);
                ax.YLim = [0 1] + y([1 end]);
                sVar(2).tickLabels = flip(sVar(2).tickLabels);
            end
            ax.XLabel.String = sVar(1).label;
            ax.YLabel.String = sVar(2).label;
            ax.XTick = sVar(1).ticks;
            ax.YTick = sVar(2).ticks;
            ax.XTickLabel = sVar(1).tickLabels;
            ax.YTickLabel = sVar(2).tickLabels;
            ax.TickLabelInterpreter = 'none';
            MPlot.Axes(ax);
        end
        
        function tl = PlotDendroHeatmap(D, ut, varargin)
            % Plot heatmap with clustering dendrograms on top and/or side
            % 
            %   tl = PlotDendroHeatmap(D, ut)
            %   tl = PlotDendroHeatmap(..., 'SortBy', '')
            %   tl = PlotDendroHeatmap(..., 'AddOnFeature', '')
            % 
            
            p = inputParser;
            p.addParameter('SortBy', '', @(x) ischar(x) || isstring(x));
            p.addParameter('AddOnFeature', '', @(x) ischar(x) || isstring(x));
            p.parse(varargin{:});
            sortName = p.Results.SortBy;
            addonFeatName = p.Results.AddFeature;
            
            % Use optimal leaf order
            [ut, nmfOrder] = sortrows(ut, {'nmfcCat', 'nmfcW'}, {'ascend', 'descend'});
            
            % Get cluster IDs
            cid = ut.nmfcCat;
            
            % Determine y-axis ordering
            N = height(ut);
            if isempty(sortName)
                featOrder = 1:N;
            else
                yVal = ut.(sortName);
                tb = table;
                tb.val = yVal;
                tb.cid = cid;
                [~, featOrder] = sortrows(tb);
                yVal = yVal(featOrder);
            end
            
            % Plot similarity heatmap
            tl = tiledlayout('flow');
            ax = nexttile;
            if isvector(D)
                Dsq = squareform(D);
            else
                Dsq = D;
            end
            Dsq2 = -Dsq(nmfOrder, nmfOrder);
            Dsq2 = Dsq2(featOrder, :);
            imagesc(Dsq2);
            MPlot.Axes(ax);
            ax.CLim = prctile(Dsq2(:), [1 99]);
            ax.Visible = 'off';
            
            % Plot clusters
            ax = nexttile('north');
            cidList = unique(cid, 'stable'); % maintain leaf order
            nClus = numel(cidList);
            [cx, cy] = MPlot.GroupRibbon(cid, [0 1], lines(nClus), 'Groups', cidList);
            text(ax, cx, cy+2, string(cidList), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8);
            ax.XLim = [1 numel(cid)] + [-.5 .5];
            ax.YLim = [0 10];
            ax.Visible = 'off';
            
            % Plot ordering feature
            ax = nexttile('west');
            if ~isempty(sortName)
                LMV.FPA.PlotFractionStack(ax, yVal, cid(featOrder), cidList, lines(nClus));
            end
            
            % Plot additional feature
            ax = nexttile('east');
            switch addonFeatName
                case 'depth'
                    yVal = ut.depth(featOrder) / 1e3; % um to mm
                    plot(ax, yVal, (1:N)', 'ko');
                    MPlot.Axes(ax);
                    ax.XLim(1) = 0;
                    ax.XTick = 1 : 6;
                    ax.XGrid = 'on';
                    ax.YLim = [0 N+1];
                    ax.YTick = [];
                    ax.YDir = 'reverse';
                    ax.XLabel.String = 'Depth (mm)';
                    
                otherwise
                    
            end
        end
        
        function PlotFractionStack(ax, G, C, catC, colorC)
            % Helper function used in PlotClustCombo
            
            % Compute cluster fractions in each group
            G = categorical(G);
            C = categorical(C);
            catG = unique(G, 'stable');
            catC = categorical(catC);
            numG = numel(catG);
            numC = numel(catC);
            N = zeros(numG, numC);
            for i = 1 : numG
                isG = G==catG(i);
                N(i,:) = histcounts(C(isG), catC, 'Normalization', 'probability');
            end
            
            % Plot groups
            [~, yG] = MPlot.GroupRibbon(-0.1+[-0.1 0], G, @prism);
            
            % Plot fraction
            cumN = [zeros(numG,1) cumsum(N,2,'omitnan')];
            plot(ax, cumN, yG, 'Color', [0 0 0 .3]);
            
            b = barh(ax, yG, N, 0.15, 'stacked', 'EdgeColor', 'none');
            for i = 1 : numel(b)
                b(i).FaceColor = colorC(i,:);
            end
            
            ax.YTick = yG;
            ax.YTickLabel = string(catG);
            ax.TickLabelInterpreter = 'none';
            ax.YLim = [1 numel(G)] + [-0.5 0.5];
            ax.YDir = 'reverse';
            ax.XLim = [-0.2 1];
            ax.XLabel.String = 'Fraction';
            MPlot.Axes(ax);
        end
        
    end
    
end

