import numpy as np
import pandas as pd
import os
"""
A collection of functions that take in the textgrid, and output phoneme and word labels. 

"""

def parse_textgrid(textgrid_lines, stim_start):
    """
    This function takes in your textgrid lines and returns parsed df with the times of each
    word onset OVERALL in the teask
    """
#     print(textgrid_lines)
    new_tup = False
    xmin, xmax, text = [], [], []
    phonelabel = []
    for t in textgrid_lines: 
#         print(t)
        white = (len(t) - len(t.lstrip(' ')))
        if white == 8: 
            if 'name' in t: 
                if 'words' in t: 
                    wordflag=True
                elif 'phones' in t: 
                    wordflag = False
        if white == 12: 
            t= t.replace(' ', '')
            t = t.replace('\n', '')
            if t.startswith('xmin'):
                xmin.append(float(t.split('=')[-1]) + stim_start)
            elif t.startswith('xmax'):
                xmax.append(float(t.split('=')[-1]) + stim_start)
            elif t.startswith('text'): 
                text.append(t.split('=')[-1][1:-1])
                phonelabel.append(not(wordflag))
    df = pd.DataFrame({'start':xmin, 
                     'end':xmax,
                     'label':text, 
                      'phone':phonelabel})
    return df

def load_timing(data, start_time):
    """
    Input: a bunch of indices corresponding to stimulus indices
    Output: A dataframe of individual components of the stimulus, and their onsets and offset
    """
    # print('ind to iloc', ind_to_iloc)
    times = parse_textgrid(data, start_time)
    res_df = times
    return res_df


def extract_phonemes(phn_output_dir,
                     sr, 
                     arr_len, offset=0,  
                     n_phones=41, join=True): 
    
    assert offset == 0
    first= True
    phones = np.zeros((arr_len, n_phones))
    if join: 
        fp = os.path.join(phn_output_dir, 'speaker1')
    else: 
        fp = phn_output_dir
    for file_ind, file in enumerate(os.listdir(fp)): 
        if not ('.TextGrid') in file: 
            continue
        with open(os.path.join(fp, file)) as f: 
            start_time = float(file.split('_')[1])
            lines = f.readlines()
            t = load_timing(lines, start_time)
            if first: 
                final_df= t
                first = False
            else: 
                final_df= pd.concat((final_df, t))
    return final_df

def phone_df_to_arr(phone_df, sr, neural_arrshape, ks_offset): 
    assert ks_offset == 0
    
    phone_feats = phone_df.loc[phone_df['phone'] == True]
    label_encoding = set(phone_feats['label'].values)
    label_enc = {v:k for k, v in enumerate(label_encoding)}
    
    # initialize the array
    res = np.zeros((neural_arrshape, len(label_enc)))
    print('shape', res.shape)
    
    for s, e, l in zip(phone_feats['start'].values, phone_feats['end'].values, phone_feats['label'].values):
        res[int(s*sr):int(e*sr), label_enc[l]] = 1
        
    return res, label_enc

def phone_df_to_word_level_feats(phone_df, sr, neural_arrshape, ks_offset): 
    assert ks_offset == 0
    
    #### This should be called word feats but I'm too lazy. 
    phone_feats = phone_df.loc[phone_df['phone'] == False]
    res = np.zeros((neural_arrshape, 1))
    w_ons = [None]*neural_arrshape
    for s, e, l in zip(phone_feats['start'].values, phone_feats['end'].values, phone_feats['label'].values):
        res[int(s*sr)] = 1
        w_ons[int(s*sr)] = l
        
    for k, w in enumerate(w_ons): 
        if w == '': 
            w_ons[k] = None
            
    ### Repeat but for offsets. 
    res_off = np.zeros((neural_arrshape, 1))
    w_offs = [None]*neural_arrshape
    for s, e, l in zip(phone_feats['start'].values, phone_feats['end'].values, phone_feats['label'].values):
        res_off[int(e*sr)] = 1
        w_offs[int(e*sr)] = l
        
    for k, w in enumerate(w_offs): 
        if w == '': 
            w_offs[k] = None
    return res, w_ons, res_off, w_offs


def get_phoneme_clusters(phone_enc_dict): 
    import matplotlib.pyplot as plt
    phone_set= list(phone_enc_dict.keys())

    print(phone_set)
    vocalics = ['AO1', 'OW0', 'IH0', 'IY1','AH2', 'AA0', 'EY1','UW1','AO2','AH0','AA2','AY2','IY0','EH1','EY2','AE1',
             'AA1','IH2', 'AE0','AY1', 'AH1', 'OW1','UW2', 'IH1','UW0','AE2', 'R', 'HH', 'ER0', 'Y',
               'ER1', 'UH1', 'L', 'OY1', 'EH2', '']
    coronals = ['S','Z','JH', 'D', 'SH', 'DH', 'TH', 'T', 'CH']
    labials = ['B', 'P', 'V', 'F', 'M', 'W']
    dorsals = ['G', 'K', 'NG', 'N']

    new_inds = []
    for v in vocalics: 
        new_inds.append(phone_set.index(v))
    for c in coronals: 
        new_inds.append(phone_set.index(c))
    for l in labials: 
        new_inds.append(phone_set.index(l))
    for d in dorsals: 
        new_inds.append(phone_set.index(d))
        
    return vocalics+coronals + labials + dorsals, new_inds


def get_sent_onset_offset(sr, shape, offset, labels): 
    assert offset == 0
    newlabs, newstarts, newends= [], [], []
    res = np.zeros((shape, 3))
    for s, e, t in zip(labels['Onset time'].values, labels['Offset time'].values, labels['Type'].values):
        sdf = words_df.loc[words_df['start'] > s-.2]
        sdf = sdf.loc[sdf['end'] < e + .2]
        try:
            s_ = np.min(sdf['start'].values)
            e_ = np.max(sdf['end'].values)
            s = s_
            e = e_
            newlabs.append(t)
            newstarts.append(s)
            newends.append(e)
            
        except Exception:
            print(sdf, s, e)
      
        res[int(s*sr), 0] = 1
        res[int(e*sr), 1] = 1
        res[int(s*sr):int(e*sr), -1] =1
        
        # TODO: get the label at the in
        
    new_df = pd.DataFrame({
        'Onset time':newstarts,
        'Offset time':newends,
        'Type':newlabs
    })
#     labels.to_csv(args['labels_fp'])
#     new_df.to_csv('/')
    new_fp = args['labels_fp'].split('.csv')[0] + '_fa.csv'
    new_df.to_csv(new_fp %(args['subject'], args['subject'], args['block']))
    print('saved at', new_fp %(args['subject'], args['subject'], args['block']))
    return res