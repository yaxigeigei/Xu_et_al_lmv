# -*- coding: utf-8 -*-
"""
Created on Wed May 15 20:19:44 2024

@author: DuXu
"""

import pandas as pd
import numpy as np
import array
from pathlib import Path
from scipy.io import savemat

from sca.models import SCA, WeightedPCA
from sca.util import get_sample_weights, get_accuracy
import torch

sca_dir = Path(r"C:\chang_lab\project_np\analysis\pop_dynamics\sca")

clus_df = pd.read_pickle(sca_dir/"df"/"clus.pkl")
resp_df = pd.read_pickle(sca_dir/"df"/"resp.pkl")

# %% Extract arrays from dataframes

def cat_elem(df):
    df = df.map(lambda x: np.array(x) if isinstance(x, array.array) else x)
    arr = df.iloc[:,1:].to_numpy()
    arr = np.array([[elem for elem in row] for row in arr])
    arr = np.transpose(arr, (2, 0, 1))
    arr = np.reshape(arr, (-1, arr.shape[2]), order='F')
    return arr

resp = cat_elem(resp_df)

X = np.copy(resp);
X = X - np.nanmean(X, axis=0)[None,:]
X = X / np.nanstd(X, axis=0)[None,:]
X[np.where(np.isnan(X))] = 0

is_resp = np.any(clus_df.loc[:,"atten":"prod"], axis=1)

# X = X[:,is_resp]
# clus_df = clus_df.loc[is_resp,:]

# %% Set up SCA

# Set the default CUDA device to GPU
torch.set_default_device(0)

def batch_sca(X, n_components=np.arange(6, 21)):
    # 
    mdls = [SCA(n_components=i) for i in n_components]
    Z = [mdl.fit_transform(X) for mdl in mdls]
    params = [mdl.params for mdl in mdls]
    
    for z, p in zip(Z, params): p['Z'] = z
    
    d = {'X': X, 'n_components': n_components}
    for i, p in zip(n_components, params): d['decomp_' + str(i)] = p
    
    return mdls, d

out_dir = sca_dir/"computed_sca" 
out_dir.mkdir(exist_ok=True)

# %% With units from all regions

sca, sca_dict = batch_sca(X[:,is_resp])
savemat(out_dir / "sca_all.mat", sca_dict)

# %% With units from each region

for region in clus_df["region"].unique():
    is_select = is_resp & (clus_df["region"]==region)
    is_select = is_select.astype("bool")
    sca, sca_dict = batch_sca(X[:,is_select])
    savemat(out_dir / "sca_{}.mat".format(region), sca_dict)

# %% With mPrCG units excluding linker units

region = 'mPrCG'
linker = 'bridge'
is_select = is_resp & (clus_df["region"]==region) & (clus_df["hcGroup"]!=linker)
is_select = is_select.astype("bool")
sca, sca_dict = batch_sca(X[:,is_select])
savemat(out_dir / "sca_{}_no-{}.mat".format(region,linker), sca_dict)

# %% Without units from a specific region

region = 'mPrCG'
is_select = is_resp & (clus_df["region"]!=region)
is_select = is_select.astype("bool")
sca, sca_dict = batch_sca(X[:,is_select])
savemat(out_dir / "sca_non-{}.mat".format(region), sca_dict)

# %% PCA

def batch_pca(X, n_components=np.arange(20, 21)):
    # 
    mdls = [WeightedPCA(n_components=i) for i in n_components]
    Z = [mdl.fit_transform(X) for mdl in mdls]
    params = [mdl.params for mdl in mdls]
    
    for z, p in zip(Z, params): p['Z'] = z
    
    d = {'X': X, 'n_components': n_components}
    for i, p in zip(n_components, params): d['decomp_' + str(i)] = p
    
    return mdls, d

out_dir = sca_dir/"computed_pca"
out_dir.mkdir(exist_ok=True)

# With units from all regions
pca, pca_dict = batch_pca(X[:,is_resp])
savemat(out_dir / "pca_all.mat", pca_dict)

# With units from each region
for region in clus_df["region"].unique():
    is_select = is_resp & (clus_df["region"]==region)
    is_select = is_select.astype("bool")
    pca, pca_dict = batch_pca(X[:,is_select])
    savemat(out_dir / "pca_{}.mat".format(region), pca_dict)

# # Without units from a specific region
# region = 'mPrCG'
# is_select = is_resp & (clus_df["region"]!=region)
# pca, pca_dict = batch_pca(X[:,is_select])
# savemat(out_dir / "pca_non-{}.mat".format(region), pca_dict)

# %%




# %% SCA with orthogonality constraint (this is slow to compute)

def batch_sca_orth(X, n_components=np.arange(12, 13)):
    # 
    mdls = [SCA(n_components=i, orth=True) for i in n_components]
    Z = [mdl.fit_transform(X) for mdl in mdls]
    params = [mdl.params for mdl in mdls]
    
    for z, p in zip(Z, params): p['Z'] = z
    
    d = {'X': X, 'n_components': n_components}
    for i, p in zip(n_components, params): d['decomp_' + str(i)] = p
    
    return mdls, d

out_dir = sca_dir/"computed_sca_orth"
out_dir.mkdir(exist_ok=True)

# With units from all regions
sca, sca_dict = batch_sca_orth(X[:,is_resp])
savemat(out_dir / "sca_all.mat", sca_dict)

# %% Sparse PCA

from sklearn.decomposition import MiniBatchSparsePCA

def batch_spca(X, n_components=np.arange(12, 13)):
    # 
    mdls = [MiniBatchSparsePCA(n_components=i, random_state=61, n_jobs=-1, batch_size=30, verbose=1) for i in n_components]
    Z = [mdl.fit_transform(X) for mdl in mdls]
    params = [{'U': mdl.components_.T, 'V': mdl.components_} for mdl in mdls]
    
    for z, p in zip(Z, params): p['Z'] = z
    
    d = {'X': X, 'n_components': n_components}
    for i, p in zip(n_components, params): d['decomp_' + str(i)] = p
    
    return mdls, d

out_dir = sca_dir/"computed_spca"
out_dir.mkdir(exist_ok=True)

# With units from all regions
spca, spca_dict = batch_spca(X[:,is_resp])
savemat(out_dir / "spca_all.mat", spca_dict)

# With units from each region
for region in clus_df["region"].unique():
    is_select = is_resp & (clus_df["region"]==region)
    is_select = is_select.astype("bool")
    pca, pca_dict = batch_spca(X[:,is_select])
    savemat(out_dir / "spca_{}.mat".format(region), pca_dict)

# %% Check input
from matplotlib import pyplot as plt

plt.figure()
plt.plot(np.hstack(resp_df["u410100245"]))

# %%

# Plot the loss over all iterations
plt.figure()
plt.plot(sca.losses)
plt.xlabel('Training Epoch')
plt.ylabel('Loss')
plt.title('Loss over training')

# # Plot the variance explained
# plt.figure()
# # plt.plot(sca_resp.explained_squared_activity)
# plt.plot(pca.explained_variance_)
# plt.ylim(0)

# # Plot the variance explained
# plt.figure()
# plt.plot(sca.explained_squared_activity)
# # plt.plot(pca.explained_variance_)
# plt.ylim(0)

# 
product = sca.params['V']@sca.params['V'].T
plt.figure()
plt.imshow(product, clim=[-1,1], cmap='RdBu')
plt.colorbar()

# %% 

#For SCA, the reconstruction loss and r2 automatically get assigned to the model after fitting as attributes
print('SCA r2:', sca.r2_score)
print('SCA reconst_loss:', sca.reconstruction_loss)

#For PCA, we use the get_accuracy function from sca.utils
[pca_r2_score, pca_reconstruction_loss] = get_accuracy(pca, X)
print('PCA r2:', pca_r2_score)
print('PCA reconstr_loss:', pca_reconstruction_loss)
