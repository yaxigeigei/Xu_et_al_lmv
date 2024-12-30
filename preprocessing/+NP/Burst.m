classdef Burst
    
    methods(Static)
        function [btTb, brTb] = MakeBurstTables(stTb)
            % Compute burst as in 
            
            
            
        end
        
        function [isBst, info] = FindBurstSpikes(st, minNumSpk)
            % Find burst and burst-related spikes from a spike train based on 
            % Kapucu et al. 2012, Front. Comput. Neurosci.
            % 
            %   [isBst, info] = FindBurstSpikes(st)
            %   [isBst, info] = FindBurstSpikes(st, minNumSpk)
            % 
            % Input:
            %   st          A numeric vector of spike times in second.
            %   minNumSpk   Minimal number of spikes required to define a burst.
            % Outputs:
            %   isBst       A n-by-2 logical array. n is the number of input spikes. 
            %               The first column indicates if spikes are burst spikes.
            %               The second column indicates if spikes are burst OR "burst-related" spikes.
            %   info        A struct that contains the following information:
            %               Xm      The ISI value where maximal CMA is reached.
            %               a       A two-element vector where the first is the tolerance coeff for 
            %                       burst ISI, and the second for burst-related ISI.
            %               Xt      A two-element vector where the first is the threshold for 
            %                       burst ISI, and the second for burst-related ISI.
            % 

            if ~exist('minNumSpk', 'var')
                minNumSpk = 3;
            end

            % Compute histogram of inter-spike intervals
            isi = diff(st(:));
            tEdges = 0 : 1e-3 : 20;
            tBin = tEdges(1:end-1) + diff(tEdges);
            H = histcounts(isi, tEdges);
            H = H';

            % Find threshold from the peak of cumulative ISI histogram
            CH = cumsum(H);
            CMA = CH ./ (1:numel(tBin))';
            [CMAm, m] = max(CMA);
            Xm = tBin(m);

            % Add tolerace to the threshold
            sk = skewness(isi);
            if sk < 1
                a = [1 0.5]; % coeffs for burst and "burst-related" ISI threshold, respectively
            elseif sk < 4
                a = [0.7 0.5];
            elseif sk < 9
                a = [0.5 0.3];
            else
                a = [0.3 0.1];
            end

            CMAtail = CMA;
            CMAtail(1:m) = NaN; % mask out bins <= m
            
            aCMAm = CMAm .* a; % compute two thresholds at the same time with array expansion
            [~, t] = min(abs(CMAtail(:) - aCMAm));
            Xt = tBin(t);

            % Find bursting masks
            isBst = isi < Xt;
            isBst = [isBst; false(size(Xt))]; % pad the end to match the number of spike times
            for i = 1 : numel(Xt)
                % Convert logical mask to boundaries
                bb = MMath.Logical2Bounds(isBst(:,i));
                
                if i == 1
                    % Exclude bursts having less than three spikes
                    bb2rm = diff(bb,1,2)+1+1 < minNumSpk; % +1 for the total # of intervals, then +1 to the # of spikes
                elseif i == 2
                    % Remove burst-related periods that do not contain any burst spikes
                    bb2rm = false(size(bb,1), 1);
                    for j = 1 : numel(bb2rm)
                        bb2rm(j) = ~any(isBst(bb(j,1):bb(j,2), 1));
                    end
                end
                bb(bb2rm,:) = [];

                % Include last burst spike
                bb(:,2) = bb(:,2) + 1;

                % Convert boundaries back to logical mask
                isBst(:,i) = MMath.Bounds2Logical(bb, size(isBst,1));
            end

            % Return info
            info = struct;
            info.H = H;
            info.ISIskew = sk;
            info.CMA = CMA;
            info.Xm = Xm;
            info.a = a;
            info.Xt = Xt;
        end
        
        function ComputeSynchrony(isBst)
            % Compute the variance-to-mean ratio for the burst signal.
            % Burst signal is the fraction of units bursting.
            % 
            % 

            if iscell(isBst)
                isBst = cat(2, isBst{:});
            end

            sync = var(isBst,2) ./ mean(isBst,2);

        end
        
    end
end

