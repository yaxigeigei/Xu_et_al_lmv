classdef Movie
    
    methods(Static)
        function [vidMat, Ye] = ExtrapolateMotionTraces(Fscale, Y, t)
            % 
            
            % Initialize variables
            cc = lines(size(Y,2));
            vidMat = struct('cdata', [], 'colormap', []);
            frIdx = 1;
            
            % Prepare extrapolation
            [nTm, nTr] = size(Y);
            isVal = ~isnan(Y);
            len = sum(isVal)';
            Ye = Y;
            
            indTr = 1 : nTr;
            for i = indTr
                indTmp = setdiff(indTr, i);
                
                % Compute temporal distance (zero if overlap)
                dt = zeros(nTr, 1);
                for j = indTmp
                    bb = MMath.Logical2Bounds(isVal(:,i) | isVal(:,j));
                    if size(bb,1) == 1
                        dt(j) = 0;
                    else
                        dt(j) = bb(2,1) - bb(1,2);
                    end
                end
                
                % Sort templates first by temporal and then by trace length
                [~, order] = sortrows([dt len], {'ascend', 'descend'});
                
                % Extrapolate traces
                for j = order(:)'
                    y = NP.Motion.ExtrapolateTrace(Fscale, Ye(:,i), Y(:,j));
                    if sum(isnan(y)) == sum(isnan(Ye(:,i)))
                        continue
                    end
                    Ye(:,i) = y;
                    
                    cla;
                    plot(t, Ye); hold on
                    plot(t, Ye(:,i), 'Color', cc(i,:), 'LineWidth', 3);
                    ax = MPlot.Axes(gca);
                    ax.XLim = [250 1100];
                    ax.YLim = [500 5500];
                    ax.XLabel.String = 'Time (s)';
                    ax.YLabel.String = 'Distance from tip (um)';
                    ax.LooseInset = [0 0 0 0];
                    
                    drawnow;
                    vidMat(frIdx) = getframe(gcf);
                    frIdx  = frIdx + 1;
                    pause(0.5);
                end
            end
            
            vidMat = cat(4, vidMat.cdata);
        end
        
        function MotionRecon(seUnit, seBk, seTask, t, dur, varargin)
            % 
            
            p = inputParser();
            p.addParameter('YWin', [0 7660], @(x) isnumeric(x) && numel(x)==2);
            p.addParameter('UnitColors', lines(seUnit.numEpochs), @(x) isnumeric(x));
            p.addParameter('PropertyBinSize', .2, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('MarkerSize', 8, @(x) isnumeric(x) && isscalar(x));
            p.parse(varargin{:});
            
            tWin = t + [-dur 0];
            yWin = p.Results.YWin;
            smWin = [-1 1] * p.Results.PropertyBinSize/2;
            cc = p.Results.UnitColors;
            
            f = gcf;
            cla;
            rFun = @(x) x/diff(yWin)*diff(tWin)/f.Position(3)*f.Position(4);
            
            % Plot background spikes
            tb = seBk.SliceTimeSeries('spikes', tWin);
            tSpk = tb.time{1};
            yc = tb.y{1};
            y = yc + seUnit.userData.F2.fwF(t*ones(size(yc)), yc);
            a = tb.amp{1};
            ampRange = 8 : 100;
            for k = 1 : numel(ampRange)
                % for each amplitude bin, plot all the spikes of that size in the
                % same shade of gray
                ind = a == ampRange(k); % the amplitudes are rounded to integers
                if ~any(ind)
                    continue
                end
                plot(tSpk(ind), y(ind), '.', ...
                    'MarkerSize', p.Results.MarkerSize, ...
                    'Color', [1 1 1] * max(0, 1-ampRange(k)/20), ... % the marker color here has been carefully tuned
                    'HitTest', 'off');
            end
            
            % Plot units
            tTb = seUnit.SliceTimeSeries('original', tWin, [], 1);
            pTb = seUnit.SliceTimeSeries('resampled', t+smWin, [], 2:4);
            
            x = cellfun(@(x) median(x), pTb.spkCentX);
            x(:) = (x-35)/24 / diff(yWin) * diff(tWin);
            yc = cellfun(@(x) median(x), pTb.spkCentY);
            y = yc + seUnit.userData.F2.fwF(t*ones(size(yc)), yc);
            amp = cellfun(@(x) median(x), pTb.spkTempAmp);
            k = 4;
            
            for i = 1 : seUnit.numEpochs
                if y(i) < yWin(1) || y(i) > yWin(2)
                    continue
                end
                
                % Plot cell
                if ~isnan(y(i))
                    NP.Movie.DrawNeuron(t+x(i), y(i), rFun(amp(i))*k, amp(i)*k, [cc(i,:) .5]);
                end
                
                % Plot spike train
                if ~isempty(tTb.time{i})
                    tSpk = tTb.time{i} + x(i);
                    MPlot.PlotPointAsLine(tSpk, y(repmat(i, size(tSpk))), amp(i)*3, 'Color', cc(i,:)); hold on
                end
            end
            plot(t+x, y, 'k.', 'MarkerSize', p.Results.MarkerSize);
            
            % Find speech text
            tRef = seTask.GetReferenceTime;
            iTrial = find(tRef>=t, 1);
            sv = seTask.GetTable('speechValue');
            txt = sv.speechText{iTrial};
            
            ax = MPlot.Axes(gca);
            ax.XLim = [tWin(1) t+dur/15];
            ax.YLim = yWin;
            ax.XLabel.String = 'Recording time (s)';
            ax.YLabel.String = 'Distance from tip (\mum)';
            ax.Title.String = txt;
            MPlot.Paperize(ax);
%             axis off
%             ax.LooseInset = [0 0 0 0];
            
            % Capture frame
            
        end
        
        function h = DrawNeuron(x, y, dx, dy, c)
            % Plot a neuron
            % 
            %   h = DrawNeuron(x, y, r, c)
            %
            % Inputs:
            %   x       X-coordinate of the center
            %   y       Y-coordinate of the center
            %   r       Radius of the circle
            %   c       Color of the circle
            % Output:
            %   h       Object handle of the circle shape. By nature, it is a rectangle with rounded corners. 
            
            px = x-dx/2;
            py = y-dy/2;
            h = rectangle('Position', [px py dx dy], 'Curvature', [1 1], 'FaceColor', c, 'LineStyle', 'none');
        end
        
        function Soundtrack(seTask, t, dur)
            % 
            
            f = gcf;
            cla;
            
            % Plot sound
            tWin = t + [-dur 0];
            anaTb = seTask.SliceTimeSeries('analog', tWin, 1, {'mic'}, 'Fill', 'bleed');
            plot(anaTb.(1){1}, anaTb.(2){1}, 'k');
            
            % Find speech text
            tRef = seTask.GetReferenceTime;
            iTrial = find(tRef>=t, 1);
            sv = seTask.GetTable('speechValue');
            txt = sv.speechText{iTrial};
            
            ax = MPlot.Axes(gca);
            ax.XLim = [tWin(1) t+dur/15];
            ax.YLim = [-1 1]*3e4;
            ax.XLabel.String = 'Recording time (s)';
            ax.YLabel.String = 'AU';
            ax.Title.String = txt;
            MPlot.Paperize(ax);
%             axis off
%             ax.LooseInset = [0 0 0 0];
            
            % Capture frame
            
        end
        
    end
end

