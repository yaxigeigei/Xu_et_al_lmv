import numpy as np
import pandas as pd
from scipy.io.wavfile import read
from audio import get_clip, detect_sounds


def find_speech_times(cue, mic_wav_file, buffer=0.1, labels=None, 
                      fallback_label='speech', threshold=1.0): 
    
    """
    Use speech cues to get out the active times from the microphones wav file.
    
    Inputs:
    cue, a dataframe of Onset time, Offset time, and Type (which label this is)
    mic_wav_file: the microphone wav file. 
    fallback_label: If we dont have the label of what was said, and dont infer it with otter.ai, 
                    then we can just put in the fallback label. 
    
    
    Outputs: a dataframe of the times when the person was speaking. 
    
    If labels are available, we can add the labels. 
    If not, they will just be the fallback_label
    """
    
    
    # Make the rep_win array from the onset and offset times. 
    rep_win = np.zeros((len(cue), 2))
    rep_win[:, 0] = cue['Offset time'].values + buffer
    rep_win[:-1, 1] = cue['Onset time'].values[1:] 
    rep_win[-1, 1] = rep_win[-1, 0] + 4 # last value
    
    ### Detect onsets using many's functions
    speech_win = []
    sample_rate, mic = read(mic_wav_file)
    for t0, t1 in rep_win: 
        clip = get_clip(mic, sample_rate, t0, t1)
        t_win = detect_sounds(clip, sample_rate, th=threshold)
        t_win += t0
        speech_win.append(t_win)
        
    speech_win = np.array(speech_win)
    print(speech_win.shape, t_win)
 
    
    ons, offs = [], []
    for ww in speech_win: 
        for t0, t1 in ww: 
            ons.append(t0)
            offs.append(t1)
            
    
    if labels is None: 
        labels = [fallback_label]*len(ons)
    print(len(ons), len(labels))
    
    return pd.DataFrame({
        'Onset time':ons,
        'Offset time':offs,
        'Type':labels
    })
        
    ### put everything into the pandas dataframe. 