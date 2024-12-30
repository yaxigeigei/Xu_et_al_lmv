import numpy as np
import pandas as pd
import cv2 as cv
import librosa
from scipy.signal import argrelextrema
import pandas as pd
import matplotlib.pyplot as plt
"""
Template matching algo. 

give a wav file, and a template. 

Look for that template throughout the file.

Return the corresponding times
"""
def template_matching(base_wav_file, target_sound, start_time, end_time):
    
    print('target sound', target_sound)
    print('base file', base_wav_file)
    
    # Get the template matching algo that i have for LMV
    # Try to make it more robust across sounds. 
    # need to find the camera sound from the wavfile recorded during the speakers

    # Step 1: load the tone
    target, sr = librosa.load(target_sound)

    
    mfccs = librosa.feature.melspectrogram(target, sr=sr, n_fft=512, hop_length=512//4)
    mfccs = np.clip(mfccs, 0, np.percentile(mfccs, 98))

    trk, sr = librosa.load(base_wav_file)
    trk = trk[int(start_time*sr):int(end_time*sr)]
    mffc_trk = librosa.feature.melspectrogram(trk, sr=sr, n_fft=512, hop_length=512//4)
    mffc_trk = np.clip(mffc_trk, 0, np.max(mfccs))
    res = cv.matchTemplate(mffc_trk, mfccs, 0)
    res = res[0]

    og_to_new = trk.shape[0]/mffc_trk.shape[-1]
    effective_sr= sr/og_to_new
    res = res - np.min(res)
    res /= np.percentile(res, 98)
    amins = argrelextrema(res, np.less, order=int(effective_sr))
    amins = amins[0]
    
    plt.plot(res)
    plt.show()
    amins = [a for a in amins if res[a] < .1]
    print('Found %d instances of the target sound' %len(amins))
    print('If you are not happy, check the cue_utils.py file and adjust the threshold.')

    cue_times = [a/effective_sr + start_time for a  in amins]
    df= pd.DataFrame({ 
        'Onset time': cue_times, 
        'Offset time': [t + 0.52 for t in cue_times], 
        'Type': ['cue']*len(cue_times)})
    return df


def many_LMV_extraction(base_wav_file):
    # TODO: implement this. 
    pass 
    #return LMV_times

