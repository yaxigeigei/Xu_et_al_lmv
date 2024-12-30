classdef LFP
    %LFP Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Static)
        function [v, D] = Highpass(v, Fs)
            % Highpass filter voltage timeseries of local field potential
            
            % Construct a digitalFilter
            D = designfilt('highpassiir', ...
                'StopbandFrequency', .5, ...
                'PassbandFrequency', 1, ...
                'StopbandAttenuation', 60, ...
                'PassbandRipple', 1, ...
                'SampleRate', Fs, ...
                'MatchExactly', 'passband');
            
            % Apply filter to data
            if ~isempty(v)
                return
            end
            dtype = class(v);
            v = filtfilt(D, double(v));
            v = cast(v, dtype);
        end
    end
end

