# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 19:54:14 2023

@author: many
"""

# %% Read CSV file into Python as Pandas DataFrame
import os
from pathlib import Path
import pandas as pd

import numpy as np
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3
%matplotlib qt5

method = 'sen14-sen4'
# method = 'sen4'
ana_dir = Path(os.getenv('NP_ROOT'))/'analysis'/'sent_resp'/method/'venn'

df = pd.read_parquet(ana_dir/'clusTb.parquet')

df['inter'] = np.nanmax([df['delay'], df['init']], axis=0)

print(df)

# %% 

def compute_frac(p1, p2, p3, th):
    n1 =    np.mean((p1 >= th) & (p2 < th)  & (p3 < th))
    n3 =    np.mean((p1 < th)  & (p2 < th)  & (p3 >= th))
    n13 =   np.mean((p1 >= th) & (p2 < th)  & (p3 >= th))
    n2 =    np.mean((p1 < th)  & (p2 >= th) & (p3 < th))
    n12 =   np.mean((p1 >= th) & (p2 >= th) & (p3 < th))
    n23 =   np.mean((p1 < th)  & (p2 >= th) & (p3 >= th))
    n123 =  np.mean((p1 >= th) & (p2 >= th) & (p3 >= th))
    
    n1 = np.round(n1*100, decimals=1)
    n3 = np.round(n3*100, decimals=1)
    n13 = np.round(n13*100, decimals=1)
    n2 = np.round(n2*100, decimals=1)
    n12 = np.round(n12*100, decimals=1)
    n23 = np.round(n23*100, decimals=1)
    n123 = np.round(n123*100, decimals=1)
    
    return n1, n3, n13, n2, n12, n23, n123

# %%

regions = ['mPrCG', 'vPrCG', 'IFG', 'STG']
th =  1. # threshold of significance level

for region in regions:
    reg_df = df[df['region'] == region]
    
    plt.figure(2, figsize=(3,2.8))
    
    # 
    nn = compute_frac(reg_df['stim'], reg_df['inter'], reg_df['prod'], th)
    plt.cla()
    v3 = venn3(subsets=nn, set_labels=('Listening', 'Speaking', 'Delay / Initiation'), set_colors=('b', 'r', 'g'))
    plt.title('{}, total = {}%'.format(region, np.round(np.sum(nn), 1)))
    fig_path = ana_dir / '{}_{}_stim_inter_prod.png'.format(method, region)
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    
    # 
    nn = compute_frac(reg_df['stim'], reg_df['delay'], reg_df['prod'], th)
    plt.cla()
    v3 = venn3(subsets=nn, set_labels=('Listening', 'Speaking', 'Delay'), set_colors=('b', 'r', 'g'))
    plt.title('{}, total = {}%'.format(region, np.round(np.sum(nn), 1)))
    fig_path = ana_dir / '{}_{}_stim_delay_prod.png'.format(method, region)
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    
    # 
    nn = compute_frac(reg_df['stim'], reg_df['init'], reg_df['prod'], th)
    plt.cla()
    v3 = venn3(subsets=nn, set_labels=('Listening', 'Speaking', 'Initiation'), set_colors=('b', 'r', 'g'))
    plt.title('{}, total = {}%'.format(region, np.round(np.sum(nn), 1)))
    fig_path = ana_dir / '{}_{}_stim_init_prod.png'.format(method, region)
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
