import os
import numpy as np
import pandas as pd
import librosa
import crepe
import scipy
import matplotlib.pyplot as plt

def linear_interpolation(x, connect_beginning=False, max_gap=200):
    f_idx = np.isfinite(x)
    loc_idx = np.where(f_idx)[0]  # index of non nan
    gap_sizes = np.diff(loc_idx)   # difference in contiguity
    gap_loc = loc_idx[:-1][(gap_sizes > 1)]   # where are the gaps
    gap_loc_size = gap_sizes[(gap_sizes > 1)]  # how big are the gaps
    # make a copy of the pitch
    pitch_x = np.copy(x)
    # add the initial line
    if connect_beginning:
        first_idx = loc_idx[0]
        pitch_x[0:first_idx] = np.linspace(0, x[first_idx], first_idx)
    # loop through gaps
    for start, gap_size in zip(gap_loc, gap_loc_size):
        # only connect if it is small enough
        if gap_size >= max_gap:
            continue
        # span
        stop = start+gap_size
        # linear line
        val_start = x[start]
        val_stop = x[stop]
        line = np.linspace(val_start, val_stop, gap_size)
        # put it all together
        pitch_x[start:stop] = line
    return pitch_x
# make a function that extracts pitch from the waveforms
def crepe_pitch(audio, sr=1000., step_size=10, confidence_threshold=0.8):
    # run crepity crepe
    time, frequency, confidence, activation = crepe.predict(audio, sr, model_capacity='full',
                                                            viterbi=True, step_size=step_size)
    # mask based on confidence
    idx = np.logical_or(confidence < confidence_threshold, np.isfinite(confidence)==False)
    # undefined regions are zero.... not sure how i feel about this
    frequency_na = np.copy(frequency)
    frequency_na[idx] = np.nan
    return frequency_na[:-1]


def make_pitch(audio_dir, pitch_dir):
    
    """
    Inputs: audio_dir - where the wav files live
    they should be in the format wav_t0_t1.wav, where t0 is the start time
    """
    import numpy as np
    wavfiles = os.listdir(audio_dir)
    print('n wav', len(wavfiles))
    for w in wavfiles: 
    #     print(w)
        if not 'wav' in w: 
            continue
        start_time = w.split('_')[1]
        start_time = float(start_time)
        start_time = np.round(start_time, 2)
        y, sr = librosa.load(os.path.join(audio_dir, w), sr=1000)
#         print('y.shape', y.shape)
        res = crepe_pitch(y)
        x = res  # crepe output

        # smooth
#         print('x.shape', x.shape)
        x_finite_idx = np.roll(np.isfinite(x), 5)[1:]
        pitch_x = linear_interpolation(x, max_gap=500)
        print(pitch_x.shape)
       
#         print('pitch_x', pitch_x.shape)
#         print('pitch_x', pitch_x)

        if pitch_x.shape[0] > 50:
        
            y = np.convolve(pitch_x, np.ones(50), mode='same')
        else: 
            y = np.convolve(pitch_x, pitch_x.shape[0], mode='same')
        if np.sum(np.isnan(y)) == y.shape[0]:
            print('using smaller window')
            y = np.convolve(pitch_x, np.ones(5), mode='same')
#             plt.plot(y)
    
#         plt.plot(y, label='y')
#         plt.show()
        
      
        try:
            y_norm = y/np.nanmax(np.abs(y))
            y_norm[:np.where(np.isfinite(y_norm))[0][0]] = y_norm[np.isfinite(y_norm)][0]
            y_norm[~np.isfinite(y_norm)] = y_norm[np.isfinite(y_norm)][-1]
            no_minmax=False
        except Exception: 
            print('Exception occured!!!')
            y_norm = np.zeros_like(y) + 1e-15
            no_minmax = True
            

        # init
        pitchMin = np.zeros(len(y_norm))
        pitchMax = np.zeros(len(y_norm))

        # just the dumb max and min
        if not no_minmax: 
            pitchMax[np.nanargmax(x)] = 1
            pitchMin[np.nanargmin(x)] = 1

        # discrete delta upwards
        diff_up = np.diff(y_norm)
        diff_up_ = np.convolve(diff_up, np.ones(50), mode='same')
        if np.sum(np.isnan(diff_up_)) == diff_up_.shape[0] or (diff_up.shape[0] != 50 and diff_up_.shape[0] == 50) :
            print('diffdiff')
            diff_up_ = np.convolve(diff_up, np.ones(5), mode='same')

        diff_up = diff_up_
#         print(diff_up 'pre', diff_up)
        diff_up[diff_up < 0] = 0
#         print('diff up FINAL', diff_up)
        if np.sum(diff_up) != 0:
            diff_up = diff_up / np.nanmax(diff_up)
        
        peak_idx, props = scipy.signal.find_peaks(diff_up, height=0.1)
        pitchUp = np.zeros(len(diff_up))
        if len(peak_idx) > 0:
            idx = np.concatenate([(np.diff(peak_idx) > 200)*1, [1]]) == 1
            peak_idx = peak_idx[idx]
            pitchUp[peak_idx] = props['peak_heights'][idx]

        # discrete delta downwards
        diff_down = np.diff(y_norm)
#         diff_down = np.convolve(diff_down, np.ones(50), mode='same')
        diff_down_ = np.convolve(diff_down, np.ones(50), mode='same')
        if np.sum(np.isnan(diff_down_)) == diff_down_.shape[0] or (diff_down.shape[0] != 50 and diff_down_.shape[0] == 50) :
            print('diffdiff')
            diff_down_ = np.convolve(diff_down, np.ones(5), mode='same')
        diff_down = diff_down_
        
       
        diff_down[diff_down > 0] = 0
        if np.sum(diff_down) !=0:
            diff_down = np.abs(diff_down / np.nanmax(np.abs(diff_down)))
        peak_idx, props = scipy.signal.find_peaks(diff_down, height=0.1)
        pitchDown = np.zeros(len(diff_down))
        if len(peak_idx) > 0:
            idx = np.concatenate([(np.diff(peak_idx) > 200)*1, [1]]) == 1
            peak_idx = peak_idx[idx]
            pitchDown[peak_idx] = props['peak_heights'][idx]
    #     print(res.shape, pitchMin.shape, pitchMax.shape, pitchUp.shape, pitchDown.shape)
        try: 
            feats_we_want = np.vstack([res[:-1], pitchMin[:-1], pitchMax[:-1], pitchUp, pitchDown])
        except Exception:
            try:
                feats_we_want = np.vstack([res, pitchMin, pitchMax, pitchUp, pitchDown])
            except Exception:
                print('dimensions error', res.shape, pitchMin.shape, 
                      pitchMax.shape, pitchUp.shape, pitchDown.shape)
                assert False
    #     print(feats_we_want.shape)
#
        this_filename = w.replace('wav_', 'pitch_')
        this_filename = w.replace('.wav', '.npy')
        np.save(os.path.join(pitch_dir, this_filename) ,feats_we_want.T)
