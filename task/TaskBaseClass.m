classdef TaskBaseClass < handle
    %TASKBASECLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = TaskBaseClass()
            %TASKBASECLASS Construct an instance of this class
        end
    end
    
    methods(Static)
        function w = AdaptAudio(w, fromFreq, toFreq, gain)
            % 
            w = resample(w, toFreq, fromFreq); % upsample to the target frequency
            w = w * gain; % increase the volume to be comparable to the cue sounds
        end
        
        function t = MakeTriggerWave(w, Fs, dur, delay)
            delaySamples = round(delay * Fs);
            durSamples = max(round(dur * Fs), 1); % must have at least one sample
            t = zeros(size(w));
            a = min(delaySamples+1, length(w));
            b = min(a+durSamples-1, length(w));
            t(a : b) = 10;
        end
        
        function ms = Millis()
            ms = round(GetSecs * 1e3);
        end
        
        function val = SampleExp(mu, muMin, muMax)
            % Sample a value from exponential distribution with limits
            
            if mu == muMin && mu == muMax
                val = mu;
                return
            end
            
            if mu < muMin || mu > muMax
                val = mu;
                warning('mu value %g is out of the range [%g, %g]. Use %g.', mu, muMin, muMax, mu);
                return
            end
            
            val = exprnd(mu);
            
            if val > muMax || val < muMin
                val = TaskBaseClass.SampleExp(mu, muMin, muMax); % resample recursively until limits are satisfied
            end
        end
        
    end
end

