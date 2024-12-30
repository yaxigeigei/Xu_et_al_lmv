import numpy as np
import pandas as pd
import glob
from pathlib import Path
import shutil

def extract_pitch(wav_dir, phone_dir, output_dir):
    """
    Launch a matlab job that extract pitch features using the STRAIGHT algorithm from VoiceSauce
    
    Inputs
        wav_dir: directory of source wav files.
        phone_dir: the directory used for the audio files from phone extraction (numbering is useful)
        output_dir: where you want the pitch files stored.
    
    """
    
    # Create output folder
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate job script from template
    this_dir = Path(__file__).parent;
    with open(this_dir / 'pitch_script_template.m') as f:
        script_lines = f.readlines()
    
    final_lines = []
    for m in script_lines:
        if '{wav_dir}' in m:
            m = m.replace('{wav_dir}', str(wav_dir))
        elif '{phone_dir}' in m:
            m = m.replace('{phone_dir}', str(phone_dir))
        elif '{output_dir}' in m:
            m = m.replace('{output_dir}', str(output_dir))
        final_lines.append(m)
    
    script_file = output_dir / 'pitch_script.m'
    with open(script_file, 'w+') as f:
        f.writelines(final_lines)
    
    # Submit the job to the server
    import subprocess
    import shlex
    
    job_name = 'STRAIGHT_pitch_extration'
    log_file = output_dir / 'run_log.txt'
    
    cmd = "submit_job -q pia-batch.q -c 8 -m 20" 
    cmd += " -o" + ' ' + str(log_file) + ' ' + '-n ' + job_name
    cmd += f" -x /data_store2/MATLAB/R2019a/bin/matlab {str(script_file)}"
    print('command!', cmd)
    cmd = shlex.split(cmd)
    subprocess.run(cmd, stderr=subprocess.STDOUT)
    print('submitted STRAIGHT_pitch_extration job')
    