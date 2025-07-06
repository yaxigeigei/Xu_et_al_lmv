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

import data

lmv_dir = Path(os.environ["NP_ROOT"])/"analysis_lmv"
data_dir = lmv_dir/"data"/"pkl_m2_ex3_sentence-avg"
sca_dir = lmv_dir/"pop_dynamics"/"sca_ls"

# Load data
task_df = pd.read_pickle(data_dir/"task.pkl")
clus_df = pd.read_pickle(data_dir/"clus.pkl")
resp_df = pd.read_pickle(data_dir/"resp.pkl")
resp_df = data.numpy_elem(resp_df)

# Set the default CUDA device to GPU
torch.set_default_device(0)

# %% Extract arrays from dataframes

# Get responses and trial timestamps
resp = data.cat_elem(resp_df)
t = resp[:,0]
X = resp[:,1:]

# Make a timeseries of trial indices
trial_idx = np.concatenate([np.ones_like(x)*i for i, x in enumerate(resp_df["time"])])

# Select only the responsive units
is_resp = np.any(clus_df.loc[:,"atten":"prod"], axis=1)
resp_clus_df = clus_df.loc[is_resp,:]
X = X[:,is_resp]

# Separate responses by listening and speaking phases
is_stim = (t > task_df["stimMatchOn"][0]) & (t < task_df["stimMatchOff"][0])
is_prod = (t > task_df["prodMatchOn"][0]) & (t < task_df["prodMatchOff"][0])
t_stim = t[is_stim]
t_prod = t[is_prod]
X_stim = X[is_stim,:]
X_prod = X[is_prod,:]

# Zscore responses
X_stim = data.zscore(X_stim, axis=0)
X_prod = data.zscore(X_prod, axis=0)

# %% 
from pipeline import batch_sca

out_dir = sca_dir/"computed_sca"
out_dir.mkdir(exist_ok=True)

# %% 

region = "mPrCG"
is_mPrCG = resp_clus_df["region"] == region
is_mPrCG = is_mPrCG.astype("bool")  # convert Pandas nullable boolean dtype to numpy boolean

n_comps = np.arange(6, 21)

stim_mdls, stim_dict = batch_sca(X_stim[:,is_mPrCG], n_components=n_comps)
savemat(out_dir / f"sca_{region}_stim.mat", stim_dict)

prod_mdls, prod_dict = batch_sca(X_prod[:,is_mPrCG], n_components=n_comps)
savemat(out_dir / f"sca_{region}_prod.mat", prod_dict)

# %% 
import plot

n_comp = 12
comp_idx = np.where(n_comps == n_comp)[0][0]

stim_result = stim_dict[f"decomp_{n_comp}"]  # dict_keys(['U', 'b_u', 'V', 'b_v', 'Z'])
prod_result = prod_dict[f"decomp_{n_comp}"]

n_trials = len(task_df)
is_win = (trial_idx == 0) & is_stim
t_plot = t[is_win]
Z = stim_result["Z"].reshape(-1, n_trials, n_comp, order="F")  # Shape (time, trial, comp)

fig, axs = plt.subplots(3, 4, figsize=(12, 5))
plot.latent_activity(axs, t_plot, Z, save_path=sca_dir/f"sca_{n_comp}-comp_{region}_stim.png")

is_win = (trial_idx == 0) & is_prod
t_plot = t[is_win]
Z = prod_result["Z"].reshape(-1, n_trials, n_comp, order="F")  # Shape (time, trial, comp)

fig, axs = plt.subplots(3, 4, figsize=(12, 5))
plot.latent_activity(axs, t_plot, Z, save_path=sca_dir/f"sca_{n_comp}-comp_{region}_prod.png")

# %% 
from scipy.linalg import subspace_angles

# Get component matrices
U_stim = stim_result["U"]  # Shape (n_units, n_comp)
U_prod = prod_result["U"]  # Shape (n_units, n_comp)

# Compute principal angles between subspaces
angles = subspace_angles(U_stim, U_prod)
angles_deg = np.rad2deg(angles)

# Plot principal angles
plt.figure(figsize=(8, 4))
plt.bar(range(1, n_comp+1), angles_deg)
plt.xlabel("Component")
plt.ylabel("Principal Angle (degrees)")
plt.title(f"Principal Angles Between {region} Stim and Prod Components")
plt.xticks(range(1, n_comp+1))
plt.ylim(0, 90)
plt.tight_layout()
plt.savefig(sca_dir/f"sca_{n_comp}-comp_{region}_principal_angles.png")
plt.show()

# Print summary statistics
print(f"Mean angle: {np.mean(angles_deg):.1f}°")
print(f"Median angle: {np.median(angles_deg):.1f}°")
print(f"Min angle: {np.min(angles_deg):.1f}°")
print(f"Max angle: {np.max(angles_deg):.1f}°")

# %% 
from sklearn.cross_decomposition import CCA, PLSCanonical

# Initialize CCA with n_comp components
# mdl = CCA(n_components=n_comp)
mdl = PLSCanonical(n_components=n_comp)

# Get latent variables
Z_stim = stim_result["Z"]  # Shape (time * trials, comp)
Z_prod = prod_result["Z"]  # Shape (time * trials, comp)

# Interpolate Z_stim to match the number of samples in Z_prod
Z_stim_interp = np.zeros_like(Z_prod)
for i in range(Z_prod.shape[1]):
    Z_stim_interp[:, i] = np.interp(np.linspace(0, 1, Z_prod.shape[0]), 
                                   np.linspace(0, 1, Z_stim.shape[0]), 
                                   Z_stim[:, i])

# Fit CCA and transform the data
mdl.fit(Z_stim_interp, Z_prod)
Z_stim_transformed, Z_prod_transformed = mdl.transform(Z_stim, Z_prod)

# Reshape transformed components for plotting
Z_stim_transformed = Z_stim_transformed.reshape(-1, n_trials, n_comp, order="F")  # Shape (time, trial, comp)
Z_prod_transformed = Z_prod_transformed.reshape(-1, n_trials, n_comp, order="F")  # Shape (time, trial, comp)

# %% 

# Compute correlation coefficients between transformed components
corr_coefs = np.zeros(n_comp)
for i in range(n_comp):
    corr_coefs[i] = np.corrcoef(Z_stim_transformed[:, i], Z_prod_transformed[:, i])[0, 1]

# Print correlation coefficients
print("\nCorrelation coefficients between transformed components:")
for i, corr in enumerate(corr_coefs, 1):
    print(f"Component {i}: {corr:.3f}")
print(f"\nMean correlation: {np.mean(corr_coefs):.3f}")

# %%

# Plot transformed components
fig, axs = plt.subplots(3, 4, figsize=(12, 5))
is_win = (trial_idx == 0) & is_stim
t_plot = t[is_win]
plot.latent_activity(axs, t_plot, Z_stim_transformed, save_path=sca_dir/f"sca-pls_{n_comp}-comp_{region}_stim.png")

fig, axs = plt.subplots(3, 4, figsize=(12, 5))
is_win = (trial_idx == 0) & is_prod
t_plot = t[is_win]
plot.latent_activity(axs, t_plot, Z_prod_transformed, save_path=sca_dir/f"sca-pls_{n_comp}-comp_{region}_prod.png")

# %% 
from sklearn.cross_decomposition import CCA, PLSCanonical

# Initialize CCA with n_comp components
pls = PLSCanonical(n_components=n_comp)

# Interpolate Z_stim to match the number of samples in Z_prod
X_stim_reg = X_stim[:,is_mPrCG]
X_prod_reg = X_prod[:,is_mPrCG]
X_stim_interp = np.zeros_like(X_prod_reg)
for i in range(X_prod_reg.shape[1]):
    X_stim_interp[:, i] = np.interp(np.linspace(0, 1, X_prod_reg.shape[0]), 
                                   np.linspace(0, 1, X_stim_reg.shape[0]), 
                                   X_stim_reg[:, i])

# Fit CCA and transform the data
pls.fit(X_stim_interp, X_prod_reg)
X_stim_transformed, X_prod_transformed = pls.transform(X_stim_reg, X_prod_reg)

# Reshape transformed components for plotting
X_stim_transformed = X_stim_transformed.reshape(-1, n_trials, n_comp, order="F")  # Shape (time, trial, comp)
X_prod_transformed = X_prod_transformed.reshape(-1, n_trials, n_comp, order="F")  # Shape (time, trial, comp)

# %% 

# Compute correlation coefficients between transformed components
corr_coefs = np.zeros(n_comp)
for i in range(n_comp):
    corr_coefs[i] = np.corrcoef(X_stim_transformed[:, i], X_prod_transformed[:, i])[0, 1]

# Print correlation coefficients
print("\nCorrelation coefficients between transformed components:")
for i, corr in enumerate(corr_coefs, 1):
    print(f"Component {i}: {corr:.3f}")
print(f"\nMean correlation: {np.mean(corr_coefs):.3f}")

# %%

# Plot transformed components
fig, axs = plt.subplots(3, 4, figsize=(12, 5))
is_win = (trial_idx == 0) & is_stim
t_plot = t[is_win]
plot.latent_activity(axs, t_plot, X_stim_transformed, save_path=sca_dir/f"cca_{n_comp}-comp_{region}_stim.png")

fig, axs = plt.subplots(3, 4, figsize=(12, 5))
is_win = (trial_idx == 0) & is_prod
t_plot = t[is_win]
plot.latent_activity(axs, t_plot, X_prod_transformed, save_path=sca_dir/f"cca_{n_comp}-comp_{region}_prod.png")









# %% With units from each region

for region in clus_df["region"].unique():
    is_select = is_resp & (clus_df["region"]==region)
    is_select = is_select.astype("bool")
    sca, sca_dict = batch_sca(X[:,is_select])
    savemat(out_dir / "sca_{}.mat".format(region), sca_dict)

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
