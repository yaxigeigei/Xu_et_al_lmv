# %% 
%load_ext autoreload
%autoreload 2

import os
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.io import savemat
import torch

from sca.models import SCA, WeightedPCA
from sca.util import get_sample_weights, get_accuracy

lmv_dir = Path(os.environ["NP_ROOT"])/"analysis_lmv"
data_dir = lmv_dir/"data"/"pkl_m2_ex3_sentence-avg"
sca_dir = lmv_dir/"pop_dynamics"/"sca"

# Load data
clus_df = pd.read_pickle(data_dir/"clus.pkl")
resp_df = pd.read_pickle(data_dir/"resp.pkl")

# Set the default CUDA device to GPU
torch.set_default_device(0)

# %% Extract arrays from dataframes
import data
from scipy.stats import zscore

resp = data.cat_elem(resp_df)
X = zscore(resp, axis=0, nan_policy='omit')

is_resp = np.any(clus_df.loc[:,"atten":"prod"], axis=1)

# X = X[:,is_resp]
# clus_df = clus_df.loc[is_resp,:]

# %% Set up SCA
from pipeline import batch_sca

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
from pipeline import batch_pca

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
from pipeline import batch_sca_orth

out_dir = sca_dir/"computed_sca_orth"
out_dir.mkdir(exist_ok=True)

# With units from all regions
sca, sca_dict = batch_sca_orth(X[:,is_resp])
savemat(out_dir / "sca_all.mat", sca_dict)

# %% Sparse PCA
from pipeline import batch_spca

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
