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
    # Determine the path of the CatGT bash
    pkg_dir = Path('/userdata/dxu/pkg/NP')
    catGT_path = pkg_dir / 'CatGT-linux/runit.sh'
    
    # Source SpikeGLX data directory
    src_dir = Path('/datastore_spirit/human/Neuropixels')
    out_dir = Path('/userdata/dxu/project_np/preproc/sglx')
    
    # Specify the original ap.bin path
    ap_bin_path = '/datastore_spirit/human/Neuropixels/NP32/NP32_B3_g0/NP32_B3_g0_imec0/NP32_B3_g0_t0.imec0.ap.bin'
    
elif sys.platform.startswith('win'):
    # Determine the path of the CatGT executable
    pkg_dir = Path('C:/Users/many/Lab/pkg/NP')
    catGT_path = pkg_dir / 'CatGT-win/CatGT.exe'
    
    # Source SpikeGLX data directory
    src_dir = Path('C:/Users/many/Lab/preproc/sglx')
    out_dir = src_dir
    
    # Specify the original ap.bin path by browsing
    root = Tk()                         # create the Tkinter widget
    root.withdraw()                     # hide the Tkinter root window
    root.attributes("-topmost", True)   # Windows specific; forces the window to appear in front
    ap_bin_path = filedialog.askopenfilename(title='Select an ap.bin file', 
                                             filetypes=[('*.ap.bin', '*.bin')], 
                                             initialdir=src_dir)
    root.destroy()
    
else:
    print('Unknown system, cannot run CatGT')

# Derive other names and paths
ap_bin_path = Path(ap_bin_path)
if ap_bin_path.is_file():
    # Extract data directory
    prb_dir = ap_bin_path.parents[0]
    run_dir = ap_bin_path.parents[1]
    subj_src_dir = ap_bin_path.parents[2]
    subj_out_dir = Path(str(subj_src_dir).replace(str(src_dir), str(out_dir)))
    
    # Extract run name from file name
    ap_bin_name = ap_bin_path.name
    suffix_start = re.search(r'_g[0-9]+_t[0-9]+\.imec[0-9]+\.ap\.bin$', ap_bin_name).start()
    run_base_name = ap_bin_name[0:suffix_start]
    
else:
    print("No ap.bin was found")


# %% Assemble CatGT command

# Assemble command
params = list()
params.append('-dir=' + str(subj_src_dir))  # subject directory of source data, e.g. '/datastore_spirit/human/Neuropixels/NP11'
params.append('-run=' + run_base_name)      # run name, e.g. 'NP11_B4'
params.append('-g=' + "0,0")                # gate
params.append('-t=' + "0,0")                # trigger
params.append('-zerofillmax=0')             # zero for not leaving any gap when concatenating
params.append('-prb=' + "0")                # probe
params.append('-ap -lf -ni')                # stream
# params.append('-aphipass=300')              # ap highpass filtering
params.append('-gblcar')                    # ap common average referencing
params.append('-gfix=0.40,0.10,0.02')       # ap artifact removal: ||amp(mV)||, ||slope(mV/sample)||, ||noise(mV)||
params.append('-SY=0,384,6,500')            # extract sync edges from imec SY (probe,word,bit,millisec)
params.append('-XA=0,1,3,500')              # extract sync edges from nidq XA (word,thresh1(V),thresh2(V),millisec)
params.append('-prb_fld -out_prb_fld')      # a folder per probe for both input and output
params.append('-dest=' + str(subj_out_dir)) # subject directory of output data, e.g. '/userdata/dxu/project_np/preproc/sglx/NP11'
params = ' '.join(params)

if sys.platform.startswith('linux'):
    catGT_cmd = ['bash', str(catGT_path), params]
else:
    catGT_cmd = ' '.join([str(catGT_path), params])

print('CatGT command:')
print(catGT_cmd)


# %% Run CatGT

if not os.path.exists(subj_out_dir):
    os.mkdir(subj_out_dir)

subprocess.run(catGT_cmd)

