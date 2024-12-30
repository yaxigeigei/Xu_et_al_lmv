import subprocess
import shlex

def invert(chunked_wav_dir, artic_dir):

    # Create a folder to save MFCC data
    mfcc_dir = artic_dir / 'mfcc'
    mfcc_dir.mkdir(parents=True, exist_ok=True)
    
    # Construct command
    log_file = artic_dir / 'invert.txt'
    string = "submit_job -q mind-gpu" 
    string += " -m 78 -g 1 -o" + ' ' + str(log_file) + ' ' + '-n ' + 'invert'
    string += ' -x /userdata/dxu/anaconda3/envs/artic_inversion/bin/python'
    string += ' /userdata/dxu/project_np/code/third_party_python/articulatory_inversion/Predictions_arti/predictions_arti.py'
    string += ' --model_name F01_indep_Haskins_loss_90_filter_fix_bn_False_0_setting2'
    string += ' --wav_folder ' + str(chunked_wav_dir)
    string += ' --mfcc_folder ' + str(mfcc_dir)
    string += ' --output_folder ' + str(artic_dir)
    
    cmd = shlex.split(string)
    subprocess.run(cmd, stderr=subprocess.STDOUT)
    print('running inversion!')
    print('returning artic output dir and prog. file')
