
import numpy as np
import pandas as pd
import os
from os.path import join
import pandas as pd
from scipy.ndimage import gaussian_filter1d

def extract_spikes(ks_out_fp, ks_sr, ks_offset=0, spike_cutoff=500, 
                  prcsd_sr=100, spike_time_fp =None, many_clusters=None): 
    
    """
    Inputs: 
    
    ks_out_fp: 
        Where the kilosort output files are.
        
    ks_sr: 
    the kilosort sampling rate (30e3 by default)
    
    many_clusters: the clusters many gives to merge. 
 
    """
    # Step 1. load the KS data
    print('using ks out fp = ', ks_out_fp)
    
    
    full_arr_filepath = ks_out_fp
    
    # Get the spike times and the units.
    
    if not  spike_time_fp is None: 
        times = np.load(join(spike_time_fp, 'spike_times.npy'))
    else: 
        times = np.load(join(full_arr_filepath, 'spike_times.npy'))
    
    
    
    units = np.load(join(full_arr_filepath, 'spike_clusters.npy'))
    df = pd.DataFrame({'times':np.squeeze(times),
                      'units':np.squeeze(units)})

    df.tail()
    spikes = pd.read_csv(join(full_arr_filepath, 'cluster_group.tsv'), sep='\t')
    spikes.head()
    print(spikes.head())
    try:
        gs = (spikes.loc[spikes['KSLabel'].isin(['good', 'mua'])]) # Good spikes. 
    except Exception: 
        print(spikes['group'].values)
        gs = (spikes.loc[spikes['group'].isin(['mua', 'good'])])
        print('gslen', len(gs))
    from collections import defaultdict

    # Make a numpy array with all the stuff in it. 
    spikedict = defaultdict()
    good_units = []
    mintime = 1e15
    maxtime = -1e15
    many_list = []
    if not many_clusters is None:
        for m in many_clusters: 
            for mm in m: 
                many_list.append(m)

    for u in gs['cluster_id'].values:
        times = df.loc[df['units'] == u]['times'].values
        print(len(times))
        if len(times) > spike_cutoff: 
            good_units.append(u)
            mintime = min(mintime, min(times))
            maxtime = max(maxtime, max(times))

        
  
    start_time = 0 
    print('maxtime', maxtime)
    end_time = maxtime/ks_sr + ks_offset
    print(end_time, 'total time')
    print('good units', good_units)
    
    
    if many_clusters is None: 
        print('many clusters is none')
        big_spiking_arr = np.zeros((int((end_time-start_time)*prcsd_sr), len(good_units)))


        for k, u in enumerate(good_units):
            times = df.loc[df['units']==u]['times'].values/30e3
            for t in times: 
                ind = int((t + ks_offset)*prcsd_sr)-1 # KS offset gets baked in here. 
                # try:
                big_spiking_arr[ind, k] += 1
                # except Exception: 
                #     continue
                    
    else: 
        
        big_spiking_arr = np.zeros((int((end_time-start_time)*prcsd_sr), len(good_units)))
        clusters = many_clusters
        print('going thru clustr process')
        used_cluster_list = []
        all_clusters = []
        used_units = []
        
        
        for cl in clusters: 
            all_clusters.extend(cl)
            
        print('all clusters', all_clusters)
            
            
        k = 0
        for u in (good_units):
            print('unit', u)
            if not (u in used_cluster_list):
                if u in all_clusters: 
                    cluster_to_use = None
                    for c in clusters: 
                        if u in c: 
                            cluster_to_use = c
                            break
                    print('u', u, 'cluster_to_use', cluster_to_use)
                else: 
                    cluster_to_use = [u]
                
                times = df.loc[df['units'].isin(cluster_to_use)]['times'].values/30e3
                for t in times: 
                    ind = int((t + ks_offset)*prcsd_sr) # KS offset gets baked in here. 
                    # try:
                    big_spiking_arr[ind, k] += 1
                    # except Exception: 
                    #     continue
                used_cluster_list.extend(cluster_to_use)
                used_units.append(cluster_to_use[0])
                k += 1
                        
        good_units = used_units
        big_spiking_arr = big_spiking_arr[:, :k]
                    
        
                
    return big_spiking_arr, good_units