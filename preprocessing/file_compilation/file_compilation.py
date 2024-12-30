import pandas as pd
import numpy as np


def get_sent_onset_offset(sr, shape, offset, labels, labels_fp, words_df): 
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
    new_fp = labels_fp.split('.csv')[0] + '_fa.csv'
    new_df.to_csv(new_fp)
    print('saved at', new_fp)
    return res


def postprocess_phonemes(phonemes_and_words_df, 
                                  phone_features, 
                                phone_enc_dict, 
                                  word_onsets, 
                                  onset_labels, 
                                  word_offsets, 
                                  offset_labels, 
                         neural_array,
                        labels, labels_fp): 
    """
    This function will basically just take the phonemes and then 
    """
    
    import nltk

    from nltk.corpus import cmudict

    nltk.download('cmudict')
    d = cmudict.dict()
    exceptdict = {
        'charlies':2, 
        'accomodated':5,
        'barbaras':3,
        'bruh':1,
        'youve':1,
        'youre':1,
        'dont':1,
        None:0,
        'thats':1,
        'didnt':2,
        'somehwere':2,
        'understandingly':5,
        'nobodys':3,
        'weve':1, 
        'thatll':2, 
        'oclock':2, 
        'ive':1,
        'theyve':1, 
        'theyre':1, 
        'theyll':1,
        'theres':1, 
        'youll': 1
    }
    for c in ['b', 'd', 'g']:
        for v in ['aa', 'oo', 'ee']:
            exceptdict[c+v] = 1

    def nsyl(word):
        try:
            return [len(list(y for y in x if y[-1].isdigit())) for x in d[word.lower()]][0]
        except Exception: 
            if word == 'activ':
                return 2
            else: 
                try:
                    return exceptdict[word]
                except Exception:
                    print('word', word)
                    return 2 # This is bad...need to switch to estimator. 

    sylctr = {}
    maxs = 0
    for w in set(onset_labels):
        ns = nsyl(w)
        sylctr[w] = ns
        maxs = max(ns, maxs)

    complexity_feature = np.zeros((len(onset_labels), maxs))

    for k, o in enumerate(onset_labels): 
        if not (o is None):
            complexity_feature[k, sylctr[o]-1] = 1

    complexity_df = pd.DataFrame(data=complexity_feature, columns= ['complex_%d' %k for k in range(complexity_feature.shape[-1])])

    # Can we derive a stress feature? 
    phonemes_and_words_df.head()

    # Now we want to get another thing: 

    words_df = phonemes_and_words_df.loc[~phonemes_and_words_df['phone']]

    words_df = words_df.loc[words_df['label'] != '']
    ustress_df = phonemes_and_words_df.loc[phonemes_and_words_df['label'].str.endswith('0')]
    stress_df = phonemes_and_words_df.loc[phonemes_and_words_df['label'].str.endswith('1')]
    secondary_df = phonemes_and_words_df.loc[phonemes_and_words_df['label'].str.endswith('2')]

    len(stress_df), len(ustress_df), len(secondary_df)

    stress_feats = np.zeros((len(neural_array), 3))



    for k, df in enumerate([ustress_df, stress_df, secondary_df]):
        starts = df['start'].values
        for s in starts:
            ind = int(s*100)
            stress_feats[ind, k] = 1

    stress_df = pd.DataFrame(data=stress_feats, columns= ['stress_%d' %k for k in range(stress_feats.shape[-1])])

    sent_feats = get_sent_onset_offset(100, neural_array.shape[0], 0, labels, labels_fp, words_df)

    phone_set= list(phone_enc_dict.keys())
    phone_cols= ['phone_' + p for p in phone_set]

    phones_df = pd.DataFrame(data=phone_features, columns=['phone_' + p for p in phone_set])


    phones_df = phones_df.rename(columns={'phone__':'phone_x'})

    stuf = phones_df.columns

    ph_cols = phones_df.columns
    ph_cols_ = []
    for p in ph_cols: 
        if not p == 'phone_': 
            ph_cols_.append(p)
    ph_cols = ph_cols_

    # Now edit the formant df. 
    vs = ['A', 'E', 'I', 'O', 'U']
    vs_lower = [v.lower() for v in vs]
    vs += vs_lower
    vs_lower += ['b','d','g','v','z','y','l','dh','m','n','ng','r','w','jh','zh']
    vowels = [p for p in ph_cols if (p.split('_')[1][0].lower() in vs_lower or p.split('_')[1][:2].lower() in vs_lower)]
    vs = np.sum(phones_df[vowels].values, axis=-1)
    vs = np.clip(vs, 0, 1)

    word_feat_dfs = pd.DataFrame(data=word_onsets, columns=['word_onset'])
    word_feat_dfs['onset_labels'] = onset_labels
    word_feat_dfs['word_offset'] = word_offsets
    word_feat_dfs['offset_labels'] = offset_labels
    sent_feats_df = pd.DataFrame(data=sent_feats, columns=['sent_onset', 'sent_offset', 'boxcar'])
    
    
    return word_feat_dfs, sent_feats_df, phones_df, stress_df, complexity_df, vowels


def extract_mxu_vid_feats(subject): 
    # Step 1. Load Many's directory
    import os
    from scipy.io import loadmat
    from scipy.stats import zscore
    if subject == 'NP30':
        m_fp = '/userdata/smetzger/data/NP30_B12_mc/many_artics_new.mat'
    elif subject =='NP11': 
        m_fp = '/userdata/smetzger/data/NP11_B4_mc_update/many_artics_new.mat'
    if subject == 32: 
        print('No file yet')

    many_artics = loadmat(m_fp)
    print(many_artics.keys())
    df = pd.DataFrame({ 
        't':np.squeeze(many_artics['t_arr']), 
        'mouth_height':np.squeeze(many_artics['h_arr']), 
        'mouth_width':np.squeeze(many_artics['w_arr']),
        'nose2chin':np.squeeze(many_artics['c_arr']),
        'll':np.squeeze(many_artics['ll_arr']),
        'ul':np.squeeze(many_artics['ul_arr'])
    })

    df = df.sort_values('t')
            # plt.show()
            # plt.plot(np.diff(np.squeeze(trace['t'])))

    #Okay now we need to download the motion artifact

    artic_df = df

    artic_df['t'] = [np.round(t, 2) for t in artic_df['t']]

    mint = min(artic_df['t'])
    maxt = max(artic_df['t'])

    mint, maxt

    artic_df['t']

    newt = np.arange(0*100, maxt*100)/100
    newm = np.empty((newt.shape[0], 5)) # New motion
    newm[:] = np.nan

    mint

    maxt

    artic_df.head()

    for k, t in enumerate(newt): #[:1000]: 
        if t in list(artic_df['t'].values): 
            # try:
            newm[k, 0] = artic_df.loc[artic_df['t'] == np.round(t, 2)]['mouth_height'].values[0]
            newm[k, 1] = artic_df.loc[artic_df['t'] == np.round(t, 2)]['mouth_width'].values[0]
            newm[k, 2] = artic_df.loc[artic_df['t'] == np.round(t, 2)]['nose2chin'].values[0]
            newm[k, 3] = artic_df.loc[artic_df['t'] == np.round(t, 2)]['ll'].values[0]
            newm[k, 4] = artic_df.loc[artic_df['t'] == np.round(t, 2)]['ul'].values[0]
            # except Exception: 
            #     print('t', t)



    artic_df_new = pd.DataFrame({'t':newt, 'mouth_height': newm[:, 0], 
                                'mouth_width': newm[:, 1],
                                'nose2chin': newm[:, 2], 
                                'll': newm[:, 3], 
                                 'ul':newm[:, 4]
                                })

    artic_df_new = artic_df_new.interpolate()
    
    return artic_df_new