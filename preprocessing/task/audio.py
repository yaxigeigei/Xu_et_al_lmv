import numpy as np
import pandas as pd
import glob
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

import scipy
from scipy.io.wavfile import read, write
from scipy.signal import hilbert, spectrogram, medfilt
import librosa
import librosa.display
from librosa.feature import melspectrogram

from readSGLX import readMeta, SampRate, makeMemMapRaw


def save_ni_audio(preproc_root, subject_block, mic_chan=1, speaker_chan=2, trig_chan=3, pdiode_chan=4):
    """
    Read speaker and mic audio from NIDQ file and save them as wav files
    """
    
    # Instead of using hardcoded path, let's search for the nidq.bin file
    ni_bin_pattern = preproc_root / subject_block / 'sglx' / '*' / '_'.join([subject_block, '*.nidq.bin'])
    results = glob.glob(str(ni_bin_pattern))
    
    if len(results) == 0:
        print('no matching nidq.bin file was found')
    elif len(results) > 1:
        print('more than one matching nidq.bin files were found', results)
    
    ni_bin_file = Path(results[0])
    print('extracting audio from', ni_bin_file)
    
    
    # Load metadata
    meta = readMeta(ni_bin_file)
    
    # Map data array to memory
    ni_bin = makeMemMapRaw(ni_bin_file, meta)
    
    # Load microphone audio
    mic = ni_bin[mic_chan,:]
    # mic = mic / np.amax(mic)
    
    # Load speaker audio
    speaker = ni_bin[speaker_chan,:]

    # Load speaker trigger signal
    trig = ni_bin[trig_chan,:]

    # Load speaker audio
    pdiode = ni_bin[pdiode_chan,:]
    
    # Save audio as WAV files
    mic_wav_file = preproc_root / subject_block / 'audio_files' / (subject_block + '_mic.wav')
    mic_wav_file.parent.mkdir(parents=True, exist_ok=True)
    
    speaker_wav_file = preproc_root /  subject_block / 'audio_files' /(subject_block + '_speaker.wav')
    speaker_wav_file.parent.mkdir(parents=True, exist_ok=True)

    pdiode_wav_file = preproc_root /  subject_block / 'audio_files' /(subject_block + '_pdiode.wav')
    pdiode_wav_file.parent.mkdir(parents=True, exist_ok=True)
    
    sample_rate = int(SampRate(meta))
    write(mic_wav_file, sample_rate, mic)
    write(speaker_wav_file, sample_rate, speaker)
    write(pdiode_wav_file, sample_rate, pdiode)
    
    print('saved mic wav', mic_wav_file)
    print('saved speaker wav', speaker_wav_file)
    print('saved pdiode wav', pdiode_wav_file)
    
    # return mic_wav_file, speaker_wav_file, mic, speaker, sample_rate


def chunk_audio(labels, audio_file, clip_dir, 
               label_offset=0.0, onset_pad=0, offset_pad=0, 
               audio_offset=0.0, 
               plot_interval=20):
    """
    Given a full audio wav file, slice it up using the times defined.
    
    Inputs
        labels dataframe with onset time and offset time of the clip of interest to be
        extracted from the audio. 
        wavfile fp (wavfiles)
    
    Outputs
        Clips of wav files will be saved to your specified directory
    """
    
    # Create the folder
    clip_dir.mkdir(parents=True, exist_ok=True)
    
    # Load the wavfile
    fs, w = read(audio_file)
    if len(w.shape) > 1: 
        w = w[:,0]
    
    # Cut the audio file into clips
    onsets = labels['Onset time'].values
    offsets = labels['Offset time'].values
    
    for k, (on, off) in enumerate(zip(onsets, offsets)):
        # Get the clip
        on += label_offset - onset_pad
        off += label_offset + offset_pad
        clip = w[int(on*fs):int(off*fs)]
        
        # Save clip
        clip_path = clip_dir / ('wav_%.3f_%.3f.wav'%(on+audio_offset, off+audio_offset))
        write(clip_path, fs, clip)
        
        # Plot clip for inspection
        if plot_interval is not None and k%plot_interval == 0:
            # Slice an extended clip for display
            n_pad = int(0.2*fs); # leave 0.2 sec before and after
            clip2plot = w[int(on*fs)-n_pad:int(off*fs)+n_pad]
            
            # 
            plt.plot(clip2plot)
            plt.title('Clip #%d' %k + ' ' + labels.iloc[k]['Type'])
            plt.axvline(n_pad)
            plt.axvline(clip2plot.shape[0]-n_pad)
            plt.show()
        
    print('wav files chunked and stored in ', clip_dir, k+1, 'files processed')


def get_clip(sig, fs, t0, t1):
    """
    Return a part of the signal based on start and end time in seconds
    """
    i0 = np.int64(np.floor(t0 * fs))
    i1 = np.int64(np.floor(t1 * fs))
    return sig[i0:i1]


def detect_tone(sig, fs, tone_freq, tone_dur, t_inc):
    """
    TBW
    """
    
    # Convert time in sec to number of samples
    detect_samples = tone_dur * fs
    sample_inc = t_inc * fs
    
    # Calculate the number of time steps to run
    i_detect = np.arange(len(sig)-detect_samples, step=sample_inc).astype(np.int64)
    n_steps = len(i_detect);
    tone_mag = np.zeros((n_steps,len(tone_freq)))
    bg_mag = np.zeros((n_steps,1))
    
    for i in range(n_steps):
        # a = int(i*sample_inc)
        # b = int(i*sample_inc+detect_samples)
        a = i_detect[i]
        b = np.int64(a+detect_samples)
        tone_mag[i,:], bg_mag[i] = compute_freq_magnitude(sig[a:b], fs, tone_freq)
    
    return tone_mag, bg_mag


def compute_freq_magnitude(sig, fs, freq):
    """
    TBW
    """
    
    # Run FFT on the signal snippet
    y = np.fft.fft(sig)
    mag = np.abs(y)
    
    # Find the bin closest to the target frequencies
    all_freq = np.fft.fftfreq(sig.size, d=1.0/fs)
    freq = np.array(freq)
    dist = all_freq.reshape((len(all_freq),1)) - freq.reshape((1,len(freq)))
    i_freq = np.argmin(np.abs(dist), axis=0)
    
    # Include some neighborghing bins to form frequency ROIs
    roi_freq = i_freq.reshape((1,len(i_freq))) + np.arange(-1,2).reshape((3,1))
    roi_freq = roi_freq.flatten()
    
    # Get a frequency limit
    max_freq = 500 # Hz
    i_cut = np.argmin(np.abs(all_freq - max_freq))
    
    # Compute background magnitude from outside frequency ROIs
    bg_mag = np.delete(mag[:i_cut], np.s_[roi_freq])
    bg_mag_sd = np.std(bg_mag)
    
#     print('Bin frequency: {} Hz'.format(f[i_target]))
#     print('Bin magnitude: {}'.format(mag[i_target]))
#     print('Back ground magnitude (std): {}'.format(bg_mag_sd))
    
#     plt.plot(f[:50], mag[:50])
#     plt.plot(f[:50], np.ones_like(f[:50])*bg_mag_sd*3)
    
    return mag[i_freq], bg_mag_sd


def compute_mel_spec_saliency(sig, fs):
    """
    TBW
    """
    
    # Compute, scale, and filter spectrum
    spec = melspectrogram(y=np.float64(sig), sr=fs)
    t = np.linspace(0, len(sig)/fs, spec.shape[1])
    low_freq_mel = 0
    high_freq_mel = 2595 * np.log10(1 + (fs/2) / 700)                     # Convert Hz to Mel
    mel_points = np.linspace(low_freq_mel, high_freq_mel, spec.shape[0])  # Equally spaced in Mel scale
    freqs = 700 * (10**(mel_points/2595) - 1)
    
    spec = librosa.power_to_db(spec)
    spec = medfilt(spec, [11,1]) # this reduces some background beeping
    
    # Extract specific frequencies
    y = np.mean(spec, axis=0)
    
    # Baseline subtraction
    y = y - np.percentile(y, 20)
    
    return y, freqs, t, spec


def compute_linear_spec_saliency(sig, fs):
    """
    TBW
    """
    
    # Compute, scale, and filter spectrum
    freqs, t, spec = spectrogram(sig, fs)
    spec = 10 * np.log10(spec)
    spec = medfilt(spec, [11,1]) # this reduces some background beeping
    
    # Extract specific frequencies
    y = np.mean(spec[10:110,:], axis=0)
    
    # Baseline subtraction
    y = y - np.percentile(y, 20)
    
    return y, freqs, t, spec


def filt_signal(sig, fs, ker_sec=0.01):
    """
    TBW
    """
    
    # Apply running average
    N = np.int64(np.rint(ker_sec * fs))
    ker = np.ones(N) / N
    sig = np.convolve(sig, ker, mode='same')
    return sig


def find_windows(x):
    """
    TBW
    """
    
    # Find windows of continuous non-zero and continuous zero values
    d = np.diff(x)
    onset = np.where(d == 1)[0]
    offset = np.where(d == -1)[0]
    on_win = np.stack((onset, offset), axis=1)
    off_win = np.stack((offset[:-1], onset[1:]), axis=1)
    return on_win, off_win


def detect_sounds(sig, fs, th, max_gap=0.5, min_sound=0.5, is_plot=False):
    """
    Detect the onset of offset times of sounds

        detect_sounds(sig, fs, th, max_gap=0.5, min_sound=0.5, is_plot=False)

    Inputs
        sig             The audio signal to detect sound from.
        fs              The sampling rate of the audio signal.
        th              The threshold of detection.
        max_gap         Maximum gap in sec allowed within each sound.
        min_sound       Minimal duration in sec of a valid sound.
        is_plot         Whether or not (default) to plot the detection.

    Outputs
        on_win          An (n,2) numpy.array of sound time windows in sec. n is the number of sounds. 
                        2 is for the onset and offset.
    """
    
    # Compute spectrum "saliency"
    y, freqs, t, spec = compute_linear_spec_saliency(sig, fs)
    fs = 1 / np.diff(t[:2])
    y = filt_signal(y, fs, ker_sec=0.05)
    
    # Threshold envelop
    y_supra = y - th
    is_supra = y_supra > 0
    is_supra = np.hstack((0, is_supra, 0))
    
    # Find sound on and off windows
    on_win, off_win = find_windows(is_supra)
    
    # Fill small gaps
    max_gap = max_gap * fs
    for a, b in off_win:
        if b-a < max_gap:
            is_supra[a:b+1] = 1
    
    # Find sound on and off windows
    on_win, off_win = find_windows(is_supra)
    
    # Remove short sound on windows
    min_sound = min_sound * fs
    for a, b in on_win:
        # print((a,b))
        if b-a < min_sound:
            is_supra[a:b+1] = 0
    
    # Find sound on and off windows
    on_win, off_win = find_windows(is_supra)
    
    # Convert sample to second
    on_win = on_win / float(fs)
    
    # for inspection only
    if is_plot:
        fig, ax = plt.subplots(2, 1, figsize=(20,6))
        
        ax[0].plot(t, y, linewidth=0.5)
        ax[0].plot(t, is_supra[:-2]*th)
        ax[0].set_ylabel('AU')
        ax[0].set_xlabel('Time (sec)')
        ax[0].set_xlim(t[[1,-1]])
        
        low_freq_mel = 0
        high_freq_mel = 2595 * np.log10(1 + (fs/2) / 700) # Convert Hz to Mel
        mel_points = np.linspace(low_freq_mel, high_freq_mel, spec.shape[0]) # Equally spaced in Mel scale
        ax[1].pcolormesh(t, mel_points, spec)
        ax[1].set_yticks(mel_points[::20].flatten())
        ax[1].set_yticklabels(np.int64(freqs[::20]).astype(str))
        
        # ax[1].pcolormesh(t, freqs, spec)
        ax[1].set_ylabel('Frequency (Hz)')
        ax[1].set_xlabel('Time (sec)')
    
    return on_win

