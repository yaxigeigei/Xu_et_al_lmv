classdef Param
    
    properties(Constant)
        % Ephys
        RP = 1.5e-3;                % refractory period in sec, within which refractory period violation occurs
        minISI = 0.5e-3;            % minimally allowed inter-spike interval in sec, within which spikes are treated as duplicates
        maxRPV = 1.5;               % maximal percent refractory period violation
        maxContam = 20;             % maximal percent contanimation
        minSNR = 3;                 % minimal waveform SNR for single-units
        minNumSpk = 300;            % minimal number of spikes for single-units
        minFracSpan = 0.5;          % minimal fraction of task where a unit is present
        fsPyrCutoff = [.4 .5];      % peak to trough time that separate putative FS and Pyr units
        
        activePrct = 99;            % percentile interspike-interval, below which the period of time is deemed active
        normAddMax = 5;             % spike rate added to the maximum
    end
    
    methods(Static)
        function ops = Enrich(numOps)
            % Default options to enriching features in se
            
            ops.description = '';
            
            % Meta
            ops.isMergeMeta = false;
            
            % Spike rate
            ops.isSpkRate = false;
            ops.spkLagInSec = 0;
            ops.spkBinSize = 2.5e-3;    % NP.Param.RP
            ops.spkKerSize = 0.015;
            ops.isSpkSpan = false;
            
            % Speech features
            ops.isFiltSpeech = false;   % Filter speech audio waveform
            ops.isMel = false;          % Compute Mel spectrogram
            ops.isPitch = false;        % Derive pitch features
            ops.isArtic = false;        % Derive articulatory features
            ops.isTimitGT = false;      % Import ground truth labels of TIMIT stim
            
            % Replicate ops
            if nargin > 0
                ops = repmat(ops, [numOps 1]);
            end
        end
        
        function ops = Resample(ops)
            % Default options for resampling data in se and outputing numeric arrays
            % e.g. used in NP.SE.GetStimArray, NP.SE.GetRespArray, NP.SE.SetStimRespArrays
            
            if ~exist('ops', 'var')
                ops = struct;
            end
            
            % Variables to resample
            ops.featVars = [];
            
            % Resampling window, resolution, and method
            ops.rsWin = [];             % resampling time window in sec, e.g. [-.5 .5]
            ops.rsBinSize = 0.01;       % resampling bin size in sec
            ops.rsShifts = 0;           % time shifts in sec, e.g. [-1 : 0.2 : 1]
            ops.rsArgs = {'Method', 'linear', 'Extrap', 'none'};
            
%             % Average and reorganize output matrices
%             ops.dimAverage = [];        % dimensions to average. 1 time, 2 variable, 3 trial
%             ops.dimCombine = [];        % dimensions to collapse. 1 time, 2 variable, 3 trial
        end
        
        function cc = GetRegionColors(regions)
            regions = cellstr(regions);
            lineColors = lines();
            for i = numel(regions) : -1 : 1
                switch regions{i}
                    case 'mPrCG'
                        cc(i,:) = lineColors(3,:); % yellow
                    case 'IFG'
                        cc(i,:) = lineColors(4,:); % purple
                    case 'STG'
                        cc(i,:) = lineColors(1,:); % blue
                    case 'vPrCG'
                        cc(i,:) = lineColors(2,:); % orange
                    otherwise
                        cc(i,:) = [0 0 0];
                end
            end
        end
        
    end
end


