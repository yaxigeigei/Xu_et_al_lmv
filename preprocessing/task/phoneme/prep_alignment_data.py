import os
import shutil
from shutil import copyfile
import string
import pandas as pd

def clean_str(txt ,lower=True): 
    if lower:
        txt= txt.lower()
    else: 
        txt = txt.upper()
    # txt = txt.translate(str.maketrans('', '', string.punctuation))
    return txt


def wavsort(wav_files):
    wavfiles_sort = sorted(wav_files, key=lambda x: float(x.split('_')[1]))
    return wavfiles_sort


# add the speaker id to each column... 
def add_speaker_id(speech_df):
    id_dictionary = {} # Will contain ids for different people (E.g. if the source is many, EC, sean)
    speaker_ids = []
    speaker_id_ct = 2
    for source, task, sent in zip(speech_df['Source'], speech_df['Task'], speech_df['Type']):
        if source == 'prod':
            speaker_ids.append(1)
        elif source == 'stim': 
            if not sent in id_dictionary:
                print('adding speaker id for source, sent', source + ',' , sent)
                
                id_dictionary[sent] = speaker_id_ct
                speaker_id_ct += 1
            speaker_ids.append(id_dictionary[sent])
        else: 
            if not source in id_dictionary:
                print('adding speaker id for source', source)
                id_dictionary[source] = speaker_id_ct
                speaker_id_ct+=1
            speaker_ids.append(id_dictionary[source])
                
    speech_df['Speaker id'] = speaker_ids
    return speech_df


def prep_wavs_and_labels(wav_fp, phn_fp, labels, audio_offset=0, format_wav=True, return_ind_to_time=False):
    
    assert audio_offset== 0 
    # Find all WAV files in wav_fp folder
    wav_files = os.listdir(wav_fp)
    wav_files = [w for w in wav_files if (not '.ipynb' in w and w.endswith('.wav'))]
    labels["Ons labels"] = ['%.3f' %(o-audio_offset) for o in labels['Onset time']]
    # Move the chunked wav files into the directory
    if not os.path.isdir(phn_fp): 
        os.mkdir(phn_fp)
    ind_to_iloc = {}
    ind_to_time = {}
    for k, w in enumerate(wavsort(wav_files)):  #### This assumes there is a wav file for each stimulus you're looking
        
        # 
        timestr = w.replace('wav_', '').replace('.wav', '')
        
        spkid=labels.iloc[k]['Speaker id'] # Get the speaker id
    
        if not os.path.isdir(os.path.join(phn_fp, f'speaker{spkid}')):
            os.mkdir(os.path.join(phn_fp, f'speaker{spkid}'))
            
        shutil.copy(os.path.join(wav_fp, w), os.path.join(phn_fp, f'speaker{spkid}', f'stim_{timestr}.wav'))
        ons_time = w.split('_')[1]
        try:
            # print('ons time', ons_time)
            sdf = labels.loc[labels['Ons labels'] == ons_time]
            ind_to_iloc[k] = sdf.index
            ind_to_time[k] = (float(sdf['Ons labels'].values[0]), sdf['Offset time'].values[0])
            sdf= sdf['Type'].values[0]
        except Exception:
            print('Oh no!!!', w, labels.iloc[:2], 'wav file did not have a start time associated in the labels')
            print(type(ons_time), type(labels['Ons labels']))
            assert False
        label_text = sdf
        with open(os.path.join(phn_fp, f'speaker{spkid}', f'stim_{timestr}.lab'), 'w') as f: 
            txt = clean_str(label_text)
            f.writelines(txt)
            
    print('prepped wav and txt files for mfa')
    
    if return_ind_to_time:
        print('returning stim indices to timepoints')
        return ind_to_time
    
    return ind_to_iloc, ind_to_time

            
def make_pronounciation_dict(lines): 
    
    # CVs arent in the pronounciation dictionary. 
    cvs = ['BAA  B AA1\n', 'GAA  G AA1\n', 'DAA D AA1\n', 
      'BOO B UW1\n', 'GOO G UW1\n', 'DOO D UW1\n', 'BEE B IY1\n', 
      'DEE D IY1\n', 'GEE G IY1\n']
    
    
    pronounce_dict = {}
    with open('./phoneme_extract/librispeech_lex.txt', 'r') as f: 
        d = f.readlines()
    partd = d
    for d in partd + cvs: 
        if not ('\t') in d: 
            word = (d.split(' '))[0]
            pronounciation = (' '.join(d.split(' ')[1:]))
            pronounce_dict[word] = pronounciation
#             print(pronounciation)
        elif ('\t') in d: 
            word = d.split('\t')[0]
            pronounciation = d.split('\t')[1]
            pronounce_dict[word] = pronounciation
            
    words = []
    for l in lines: 
        new = clean_str(l, lower=False)
#         print(new)
        for n in new.split(' '): 
            if not n in pronounce_dict: 
                for k in range(1):
                    print('Missing pronounciation: add in pronounciation for', n)
                    assert False # break the code until 
            else: 
                words.append(n)
    words = sorted(list(set(words)))

    
    # Get the text files
    final_lines = []
    for w in words: 
        final_str = w
        final_str += '  '
        final_str += pronounce_dict[w]
        final_lines.append(final_str)
    with open('/userdata/smetzger/prod_feat_extract/phoneme_extract/custom_lexicon.txt', 'w') as f: 
        f.writelines(final_lines)
        
    print('pronounciation dict made')
