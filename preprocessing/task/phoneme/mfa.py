import numpy as np
import pandas as pd
import glob
from pathlib import Path
import os
import shutil
import shlex
import subprocess

from .prep_alignment_data import add_speaker_id, prep_wavs_and_labels, make_pronounciation_dict
from .edit_alignment_template import make_alignment_files, make_spanish_alignment_files

def run_mfa(chunk_wav_dir, speech_df, phone_dir, language='english'):
    
    # Label each speech with speaker identity
    speech_df = add_speaker_id(speech_df)
    
    # Prepare the .wav clips and .lab labels
    input_dir = phone_dir
    input_dir.mkdir(parents=True, exist_ok=True)
    phone_ind_to_iloc, ind_to_time = prep_wavs_and_labels(chunk_wav_dir, input_dir, speech_df)
    
    # Customize bash script
    this_dir = Path(__file__).parent
    output_dir = phone_dir / 'results'
    # aligner_env = '/home/smetzger/.conda/envs/aligner'
    aligner_env = '/userdata/qgreicius/conda_envs/mfa'
    
    if language == 'english':
        lexicon_file = this_dir / 'librispeech_lex.txt'
        make_alignment_files(phone_dir, lexicon_file, output_dir, language, aligner_env)
        script_name = 'run_aligner.sh'
        
    elif language == 'spanish':
        # lexicon_file = this_dir / 'spanish_dict.txt'
        lexicon_file = this_dir / 'spanish_lex.txt'
        make_spanish_alignment_files(phone_dir, lexicon_file, output_dir, language, aligner_env)
        script_name = 'run_aligner_spanish.sh'
    
    # Run MFA
    script_file = this_dir / script_name
    subprocess.run(['bash', str(script_file)], stderr=subprocess.STDOUT)
    
    
def pool_speakers(phone_dir):
    # Make copy of speaker files in a common folder
    
    base_dir = phone_dir / 'results'
    all_dir = base_dir / 'all'
    all_dir.mkdir(parents=True, exist_ok=True)
    for speaker in base_dir.glob('speaker*'):
        spk_dir = base_dir / speaker
        for f in spk_dir.glob('*'):
            f = f.name
            shutil.copy(spk_dir/f, all_dir/f)
    
    base_dir = phone_dir
    all_dir = base_dir / 'all'
    all_dir.mkdir(parents=True, exist_ok=True)
    for speaker in base_dir.glob('speaker*'):
        spk_dir = base_dir / speaker
        for f in spk_dir.glob('*'):
            f = f.name
            shutil.copy(spk_dir/f, all_dir/f)
            