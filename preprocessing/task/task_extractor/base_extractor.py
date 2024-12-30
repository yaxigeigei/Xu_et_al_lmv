import os
from scipy.io.wavfile import read
from audio import *

class BaseExtractor():
    """
    Base class of a task extractor
    """
    
    def __init__(self, preproc_root, subject_block, time_ranges) -> None: 
        """
        Constructor of the base task extractor
        
            BaseExtractor(preproc_root, subject_block, time_ranges)
        
        Inputs
            preproc_root    The audio signal to extract from.
            subject_block   The sampling rate of the audio.
            time_ranges     A tuple of task time range in sec, or a list of such tuples for discontinuous ranges.
        """
        
        # Ensures time_tuples is a list of tuple(s)
        if type(time_ranges) is not list:
            time_ranges = [time_ranges]
        
        self.preproc_root = preproc_root
        self.subject_block = subject_block
        self.time_ranges = time_ranges
    
    
    def extract_sound_wins(self, audio, sample_rate, time_ranges=None, th=7.0, pad=0.05, max_gap=0.4, min_sound=0.2, plot=None, plot_ind=range(5)):
        """
        Detect the onset of offset times of sounds
        
            extract_sound_wins(audio, sample_rate, time_ranges=None, th=5.0, pad=0.05, min_sound=0.2, max_gap=0.4, plot='sound', plot_ind=range(5))
        
        Inputs
            audio           The audio signal to extract from.
            sample_rate     The sampling rate of the audio.
            time_ranges     A tuple of time range of interest in sec, or a list of such tuples for discontinuous ranges.
                            If this input is None, self.time_ranges will be used by default.
            th              The threshold of detection.
            pad             The amount of time to pad before and after the detected windows.
            max_gap         Maximum gap in sec allowed within each sound.
            min_sound       Minimal duration in sec of a valid sound.
            plot            Plot each 'sound' (default) or 'range'. Use None to disable plots.
            plot_ind        Indices of the sound or range windows to plot.
            
        Outputs
            t_wins          A list of time window vectors in sec.
        """
        
        # Use self.time_ranges by default
        if time_ranges is None:
            time_ranges = self.time_ranges
        
        # Extract sound windows
        t_wins = []
        for i, tt in enumerate(time_ranges):
            # Plot detection for example time ranges
            is_plot = plot=='range' and np.isin(i, plot_ind)
            if is_plot:
                print('Example time range: #{}'.format(i))
            
            t0, t1 = tt
            clip = get_clip(audio, sample_rate, t0, t1)
            t_win = detect_sounds(clip, sample_rate, th=th, min_sound=min_sound, max_gap=max_gap, is_plot=is_plot)
            t_win = t_win + [-pad, pad] + t0
            t_wins.append(t_win)
        
        # Plot detection for example sounds
        if plot=='sound':
            sound_wins = np.vstack(t_wins)
            for i in plot_ind:
                print('Example sound: #{}'.format(i))
                clip = get_clip(audio, sample_rate, sound_wins[i,0]-0.5, sound_wins[i,1]+0.5)
                detect_sounds(clip, sample_rate, th=th, min_sound=min_sound, max_gap=max_gap, is_plot=True)
        
        return t_wins
    