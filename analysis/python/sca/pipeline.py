import numpy as np
from sca.models import SCA, WeightedPCA
from sca.util import get_sample_weights, get_accuracy
from sklearn.decomposition import MiniBatchSparsePCA






def batch_sca(X, n_components=np.arange(6, 21)):
    # initialize SCA models for different numbers of components
    mdls = [SCA(n_components=i) for i in n_components]

    # Fit models to data
    Z = [mdl.fit_transform(X) for mdl in mdls]
    
    # Get parameters from models
    params = [mdl.params for mdl in mdls]

    # Add projections to parameters
    for z, p in zip(Z, params): p['Z'] = z

    # Create dictionary of results
    d = {'X': X, 'n_components': n_components}
    for i, p in zip(n_components, params): d['decomp_' + str(i)] = p
    
    return mdls, d


def batch_pca(X, n_components=np.arange(20, 21)):
    # Initialize PCA models for different numbers of components
    mdls = [WeightedPCA(n_components=i) for i in n_components]

    # Fit models to data
    Z = [mdl.fit_transform(X) for mdl in mdls]

    # Get parameters from models
    params = [mdl.params for mdl in mdls]

    # Add projections to parameters
    for z, p in zip(Z, params): p['Z'] = z

    # Create dictionary of results
    d = {'X': X, 'n_components': n_components}
    for i, p in zip(n_components, params): d['decomp_' + str(i)] = p
    
    return mdls, d

    
def batch_sca_orth(X, n_components=np.arange(12, 13)):
    # 
    mdls = [SCA(n_components=i, orth=True) for i in n_components]
    Z = [mdl.fit_transform(X) for mdl in mdls]
    params = [mdl.params for mdl in mdls]
    
    for z, p in zip(Z, params): p['Z'] = z
    
    d = {'X': X, 'n_components': n_components}
    for i, p in zip(n_components, params): d['decomp_' + str(i)] = p
    
    return mdls, d


def batch_spca(X, n_components=np.arange(12, 13)):
    # 
    mdls = [MiniBatchSparsePCA(n_components=i, random_state=61, n_jobs=-1, batch_size=30, verbose=1) for i in n_components]
    Z = [mdl.fit_transform(X) for mdl in mdls]
    params = [{'U': mdl.components_.T, 'V': mdl.components_} for mdl in mdls]
    
    for z, p in zip(Z, params): p['Z'] = z
    
    d = {'X': X, 'n_components': n_components}
    for i, p in zip(n_components, params): d['decomp_' + str(i)] = p
    
    return mdls, d