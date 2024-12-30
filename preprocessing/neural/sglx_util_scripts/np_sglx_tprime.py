# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 19:28:09 2021

@author: Many
"""

import sys
import subprocess
import os
from pathlib import Path
from tkinter import Tk
from tkinter import filedialog
import re
import numpy as np


# %% Paths and names

if sys.platform.startswith('linux'):
    # Determine the path of the TPrime bash
    pkg_dir = Path('/userdata/dxu/pkg/NP')
    tPrime_path = pkg_dir / 'TPrime-linux/runit.sh'
    
    # Specify the catgt ap.meta path
    ap_meta_path = '/userdata/dxu/project_np/preproc/sglx/NP32/catgt_NP32_B2_g0/NP32_B2_g0_imec0/NP32_B2_g0_tcat.imec0.ap.meta'
    
    # Specify the Kilosort folder
    ks_dir = Path('/userdata/dxu/project_np/preproc/kilosort/NP30_B12_mc')
    
elif sys.platform.startswith('win'):
    # Determine the path of the TPrime executable
    pkg_dir = Path('C:/Users/many/Lab/pkg/NP')
    tPrime_path = pkg_dir / 'TPrime-win/TPrime.exe'
    
    # Specify the catgt ap.bin path by browsing
    root = Tk()                         # create the Tkinter widget
    root.withdraw()                     # hide the Tkinter root window
    root.attributes("-topmost", True)   # Windows specific; forces the window to appear in front
    ap_meta_path = filedialog.askopenfilename(title='Select the ap.meta file from CatGT output', \
                                              filetypes=[('*.ap.meta', '*.meta')], \
                                              initialdir=Path('C:/Users/many/Lab/preproc/sglx'))
    root.destroy()
    
    # Specify the Kilosort folder
    ks_dir = Path('C:/Users/many/Lab/preproc/kilosort/NP28_B1_mc')
    
else:
    print('Unknown system, cannot run TPrime')

# Derive other names and paths
ap_meta_path = Path(ap_meta_path)
if ap_meta_path.is_file():
    # Extract data directory
    prb_dir = ap_meta_path.parents[0]
    run_dir = ap_meta_path.parents[1]
    subj_dir = ap_meta_path.parents[2]
    
    # Extract run name from file name
    ap_meta_name = ap_meta_path.name
    suffix_start = re.search(r'\.imec[0-9]+\.ap\.meta$', ap_meta_name).start()
    run_name = ap_meta_name[0:suffix_start]
    
else:
    print("No ap.meta was found")


# %% TPrime

from DemoReadSGLXData.readSGLX import readMeta, SampRate#, makeMemMapRaw, GainCorrectIM, GainCorrectNI

# Path of extracted NI sync times
tostream_path = run_dir / (run_name + '.nidq.XA_0_500.txt')

# Path of extracted probe sync times
fromstream_path = prb_dir / (run_name + '.imec0.ap.SY_384_6_500.txt')

# Get the AP sampling rate
meta = readMeta(ap_meta_path)
sampling_rate = SampRate(meta)

# Convert spike times from sample index to time in sec
st_ind_path = ks_dir / 'spike_times.npy'
st_ind = np.load(st_ind_path)
st_sec = st_ind / sampling_rate
st_sec_path = ks_dir / 'spike_times_in_sec.npy'
np.save(st_sec_path, st_sec)

# Path of adjusted spike times in sec
st_sec_adj_path = ks_dir / 'spike_times_in_sec_adj.npy'

# Assemble command
params = list()
params.append('-syncperiod=1.0')
params.append('-tostream=' + str(tostream_path))
params.append('-fromstream=1,' + str(fromstream_path))
params.append('-events=1,' + str(st_sec_path) + ',' + str(st_sec_adj_path))
params = ' '.join(params)

if sys.platform.startswith('linux'):
    catGT_cmd = ['bash', str(tPrime_path), params]
else:
    tPrime_cmd = ' '.join([str(tPrime_path), params])

print('TPrime command:')
print(tPrime_cmd)


# %% Run tPrime

subprocess.run(tPrime_cmd)


# %% Load adjusted spike times

st_sec_adj = np.load(st_sec_adj_path)


# %% Compare original and the adjusted spike times

import matplotlib.pyplot as plt

st_sec_diff = st_sec_adj - np.squeeze(st_sec)

fig, ax = plt.subplots()
ax.plot(st_sec_adj/60, st_sec_diff*1000)
ax.set(xlabel='Recording time (min)', ylabel='Time shift (ms)')

