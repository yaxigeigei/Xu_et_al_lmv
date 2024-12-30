import numpy as np
import pandas as pd
import glob
from pathlib import Path
import os
import shutil
import shlex
import subprocess


def extract_intensity(wav_dir, output_dir):
    """
    Launch a matlab job to execute find_peakRate function that extracts envelope, peak envelope, 
    and peak rate of envelope change using the method in Oganian & Chang, 2019
    
    Inputs
        wav_dir: directory of source wav files.
        output_dir: where you want the output files stored.
    
    """
    
    # Create output folder
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate job script from template
    this_dir = Path(__file__).parent;
    with open(this_dir / 'intensity_script_template.m') as f:
        script_lines = f.readlines()
    
    final_lines = []
    for m in script_lines:
        if '{wav_dir}' in m:
            m = m.replace('{wav_dir}', str(wav_dir))
        elif '{output_dir}' in m:
            m = m.replace('{output_dir}', str(output_dir))
        final_lines.append(m)
    
    script_file = output_dir / 'intensity_script.m'
    with open(script_file, 'w+') as f:
        f.writelines(final_lines)
    
    # Submit the job to the server
    import subprocess
    import shlex
    
    job_name = 'intensity_extration'
    log_file = output_dir / 'run_log.txt'
    
    cmd = "submit_job -q pia-batch.q -c 8 -m 20" 
    cmd += " -o" + ' ' + str(log_file) + ' ' + '-n ' + job_name
    cmd += f" -x /data_store2/MATLAB/R2019a/bin/matlab {str(script_file)}"
    print('command!', cmd)
    cmd = shlex.split(cmd)
    subprocess.run(cmd, stderr=subprocess.STDOUT)
    print('submitted intensity_extration job')
    

def setup_intensity_script(corpus_dir, pr_output_dir):
    """
    Launch a matlab job that extract intensity features including envelop, peak envelop, peak rate.
    
    Inputs
        corpus_dir: the directory used for the audio files from phone extraction (numbering is useful)
        pr_output_dir: where you want the peakrate files stored.
    
    """
    
    # get number of trials
    maxint = -1
    for d in os.listdir(corpus_dir + '/speaker1/'):
        if '.wav' in d:
            maxint = max(maxint, int(d.split('_')[1].split('.')[0]))
    print('ntrials', maxint)
    
    # Create output folder
    if not os.path.isdir(pr_output_dir):
        os.mkdir(pr_output_dir)
    
    wavfile_dir = corpus_dir + '/speaker1/'
    print(wavfile_dir)
    
    # Generate job script from template.
    curdir = os.path.dirname(__file__)
    with open(os.path.join(curdir, 'matlab_wrapper_template.m')) as f:
        matlab_lines = f.readlines()
    
    final_lines = []
    for m in matlab_lines: 
        if '{wavfile_dir}' in m: 
            m =m.replace('{wavfile_dir}', wavfile_dir)
        elif '{output_dir}' in m: 
            if not pr_output_dir[-1] == '/':
                pr_output_dir += '/'
            m = m.replace('{output_dir}', pr_output_dir)
        elif '{nstims}' in m: 
            m = m.replace('{nstims}', str(maxint))
        final_lines.append(m)
        
    with open(os.path.join(curdir, 'matlab_wrapper.m'), 'w') as f: 
        f.writelines(final_lines)
    
    
    # Submit the job to the server
    import subprocess
    import shlex

    name = 'peakrate_extract'
    filename = f'{curdir}/run.txt'
    string = "submit_job -q pia-batch.q -c 8 -m 20" 
    string += " -o" + ' ' + filename + ' ' + '-n ' + name
    string += f" -x /data_store2/MATLAB/R2018b/bin/matlab {curdir}/matlab_wrapper.m"
    print('command!', string)
    cmd = shlex.split(string)
    subprocess.run(cmd, stderr=subprocess.STDOUT)
    print('submitted peakrate job')
        
    return filename
