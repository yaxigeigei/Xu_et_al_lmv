import numpy as np
import pandas as pd
import sys
import os
import glob
import re
from pathlib import Path
from scipy.io.wavfile import read, write
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

from .base_extractor import BaseExtractor
from audio import *


class LMVExtractor(BaseExtractor):
    """
    Extract labels for the LMV task
    
    First find audio event times which include:
    1. Cue onset and offset times
    2. Stim onset and offset times
    2. Speech onset and offset times
    
    """
    
    def __init__(self, preproc_root, subject_block, time_ranges, cue_freq=None, cue_dur=None):
        """
        Constructor of LMVExtractor
        
            LMVExtractor(preproc_root, subject_block, time_ranges, cue_freq=[392,262,330], cue_dur=0.3)
        
        Inputs:
            preproc_root    The audio signal to extract from.
            subject_block   The sampling rate of the audio.
            time_ranges     A tuple of task time range in sec, or a list of such tuples for discontinuous ranges.
            cue_freq        A list of cue tone frequencies in Hz.
            cue_dur         The duration of the tones in sec.
        """
        
        # invoking the __init__ of the parent class
        super().__init__(preproc_root, subject_block, time_ranges)
        
        # Read event csv file
        self.log_df = self.read_log_csv()
        
        # Cue duration changed from 0.5 to 0.3 since NP44_B2 (inclusive)
        if cue_dur is None:
            match = re.search(r'NP(\d+)_', subject_block)
            if int(match.group(1)) >= 44:
                cue_dur = 0.3
            else:
                cue_dur = 0.5
                print("This recording uses old cue duration of 0.5 sec")
        self.cue_dur = cue_dur
        
        # Tone frequencies of the cues
        if cue_freq is None:
            cue_freq = [392, 262, 330] # Hz
        self.cue_freq = cue_freq
        
    
    def read_log_csv(self):
        """
        
        """
        
        # Instead of using hardcoded path, let's search for the _event.csv file
        csv_pattern = self.preproc_root / self.subject_block / 'labels' / 'lmv' / '_'.join([self.subject_block, '*.csv'])
        print('Search CSV file with pattern: ', csv_pattern)
        results = glob.glob(str(csv_pattern))
        if len(results) == 0:
            print('Found no matching task file')
            assert False
        elif len(results) > 1:
            print('Found more than one matching files', results)
        
        # Read the first matching csv file
        csv_file = Path(results[0])
        print('Reading the task event csv file', csv_file)
        log_df = pd.read_csv(csv_file)
        pd.set_option('display.width', 1000)
        print(log_df)
        
        return log_df
    
    
    def extract_cue(self, audio, sample_rate, t_inc=0.01, peak_th=3e6):
        """
        Detect the onset of offset times of tones in audio
        
            extract_cue(audio, sample_rate, tone_freq=[392,262,330], tone_dur=0.3, t_inc=0.025, peak_th=3e6)
        
        Inputs
            audio           The audio signal to extract from.
            sample_rate     The sampling rate of the audio.
            t_inc           The sliding increment of the detection window.
            peak_th         The threshold of detection. The default value should be pretty robust.
            
        Outputs
            self.tone_on    A list of vectors. Each vector stores the onset times of the tone corresponding to elements in tone_freq.
            self.tone_off   Same as tone_on but for the tone offset times.
        """
        
        tone_dur = self.cue_dur
        tone_freq = self.cue_freq
        
        # Trim audio
        t = np.arange(len(audio))/sample_rate
        is_task = np.zeros(len(audio), dtype=bool)
        for time_range in self.time_ranges:
            is_task[(t >= time_range[0]) & (t < time_range[1])] = True
        t = t[is_task]
        audio = audio[is_task]
        
        # Make detection timestamps
        detect_samples = tone_dur * sample_rate
        sample_inc = t_inc * sample_rate
        i_detect = np.arange(len(audio)-detect_samples, step=sample_inc).astype(np.int64)
        t_detect = t[i_detect]
        
        # Run detection
        tone_mag, bg_mag = detect_tone(audio, sample_rate, tone_freq, tone_dur, t_inc)
        
        # Agressive background subtraction
        tone_supra = tone_mag - bg_mag * 5
        tone_supra[tone_supra < 0] = 0
        
        # Find peak magnitudes, i.e. cue onsets
        peak_ind = []
        peak_time = []
        
        for i in range(len(tone_freq)):
            ind, _ = find_peaks(tone_supra[:,i], height=peak_th, distance=tone_dur/t_inc)
            peak_ind.append(ind)
            peak_time.append(t_detect[ind])
        
        # print('The task log contains {} trials.'.format(self.log_df.shape[0]))
        [print('Detected {} tones at {} Hz'.format(len(x[0]), x[1])) for x in zip(peak_time, tone_freq)]
        
        # Plot results
        figure(figsize=(20,4))
        plt.plot(t_detect, tone_supra)
        plt.plot(peak_time[0], tone_supra[peak_ind[0],0], 'x')
        plt.plot(peak_time[2], tone_supra[peak_ind[2],2], 'x')
        plt.plot(t_detect, np.ones_like(t_detect)*peak_th, '--', color=np.ones(3)*0.5)
        plt.xlim(t_detect[[0,-1]])
        
        # Outputs
        self.cue_on = peak_time
        self.cue_off = [pk+tone_dur for pk in peak_time]
    
        
    def extract_stim(self, audio, sample_rate, th=5.0, pad=0.05, plot_ind=range(0,21,5)):
        """
        Detect the onset of offset times of speech stimulation
        
            extract_stim(audio, sample_rate, th=5.0, plot_ind=range(5))
        
        Inputs
            audio           The audio signal to extract from.
            sample_rate     The sampling rate of the audio.
            th              The threshold of detection.
            pad             The amount of time to pad before and after the detected windows.
            plot_ind        Indices of trials to plot.
            
        Outputs
            self.stim_win   A list of vectors. Each vector stores the onset times of the tone corresponding to elements in tone_freq.
        """
        
        # Get onset and offset times of the detection windows
        t0 = self.cue_off[0]+pad
        t1 = self.cue_on[2]-pad
        
        # Limit window to a max duration
        max_win_dur = 3.0 # all stim should end in 3s
        is_out = (t1 - t0) > max_win_dur
        t1[is_out] = t0[is_out] + max_win_dur
        
        # Detect sounds within the windows
        det_win = list(zip(t0, t1))
        self.stim_win = super().extract_sound_wins(audio, sample_rate, time_ranges=det_win, th=th, pad=pad, plot='range', plot_ind=plot_ind)
    
    
    def extract_prod(self, audio, sample_rate, th=3.0, pad=0.05, max_win_dur=10.0, plot_ind=range(0,21,5)):
        """
        Detect the onset of offset times of speech production
        
            extract_prod(audio, sample_rate, th=3.0, pad=0.05, max_win_dur=10.0, plot_ind=range(5))
        
        Inputs
            audio           The audio signal to extract from.
            sample_rate     The sampling rate of the audio.
            th              The threshold of detection.
            pad             The amount of time to pad before and after the detected windows.
            max_win_dur     The maximal duration of detection for each trial.
            plot_ind        Indices of trials to plot.
            
        Outputs
            self.stim_win   A list of vectors. Each vector stores the onset times of the tone corresponding to elements in tone_freq.
        """
        
        # Get onset and offset times of the detection windows
        t0 = self.cue_off[2]+pad
        t1 = np.append(self.cue_on[0][1:], np.Inf)-pad
        
        # Limit window to a max duration
        is_out = (t1 - t0) > max_win_dur
        t1[is_out] = t0[is_out] + max_win_dur
        
        # Detect sounds within the windows
        det_win = list(zip(t0, t1))
        self.prod_win = self.extract_sound_wins(audio, sample_rate, time_ranges=det_win, th=th, pad=pad, plot='range', plot_ind=plot_ind)
    
    
    def write_timing_files(self):
        
        preproc_root = self.preproc_root
        subject_block = self.subject_block
        log_df = self.log_df
        
        # Create labels for cues
        cue_name = ['cue{}'.format(x+1) for x in range(len(self.cue_on))]
        cue_label = [[pair[0] for v in pair[1]] for pair in zip(cue_name, self.cue_on)]
        cue_df = pd.DataFrame({
            'label': np.hstack(cue_label),
            'on': np.hstack(self.cue_on),
            'off': np.hstack(self.cue_off)
        })
        cue_df.sort_values(by=['on'])
        
        # Cue timing file
        cue_timing_file = preproc_root / subject_block / 'labels' / 'lmv' / 'cues_timing.txt'
        cue_timing_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(cue_timing_file, 'w') as f:
            for index, row in cue_df.iterrows():
                s = '\t'.join(['%.6f'%row['on'], '%.6f'%row['off'], row['label']]) + '\n'
                f.write(s)
        
        # Stim timing file
        stim_timing_file = preproc_root / subject_block / 'labels' / 'lmv' / 'speech_stim_timing.txt'
        stim_timing_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(stim_timing_file, 'w') as f:
            for ww, sen in zip(self.stim_win, log_df['label']):
                for t0, t1 in ww:
                    s = '\t'.join(['%.6f'%t0, '%.6f'%t1, sen]) + '\n'
                    f.write(s)
        
        stim_timing_file = preproc_root / subject_block / 'labels' / 'lmv' / 'stim_timing.txt'
        
        with open(stim_timing_file, 'w') as f:
            for ww, sen in zip(self.stim_win, log_df['stim_id']):
                for t0, t1 in ww:
                    s = '\t'.join(['%.6f'%t0, '%.6f'%t1, sen]) + '\n'
                    f.write(s)
        
        # Speech on/off time
        speech_timing_file = preproc_root / subject_block / 'labels' / 'lmv' / 'speech_prod_timing_auto.txt'
        speech_timing_file.parent.mkdir(parents=True, exist_ok=True)

        with open(speech_timing_file, 'w') as f:
            for ww, sen in zip(self.prod_win, log_df['label']):
                for t0, t1 in ww:
                    s = '\t'.join(['%.6f'%t0, '%.6f'%t1, sen]) + '\n'
                    f.write(s)
