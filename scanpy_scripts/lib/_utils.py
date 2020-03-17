"""
Utility functions
"""

import numpy as np
import scipy.sparse as sp
import pandas as pd
import scanpy as sc

def run_harmony(
        adata,
        batch,
        theta=2.0,
        use_rep='X_pca',
        key_added='hm',
        random_state=0,
):
    if not isinstance(batch, (tuple, list)):
        batch = [batch]
    if not isinstance(theta, (tuple, list)):
        theta = [theta]
    for b in batch:
        if b not in adata.obs.columns:
            raise KeyError(f'{b} is not a valid obs annotation.')
    if use_rep not in adata.obsm.keys():
        raise KeyError(f'{use_rep} is not a valid embedding.')
    meta = adata.obs[batch].reset_index()
    embed = adata.obsm[use_rep]

    # ===========
    import rpy2.robjects
    from rpy2.robjects.packages import importr
    harmony = importr('harmony')
    from rpy2.robjects import numpy2ri, pandas2ri
    numpy2ri.activate()
    pandas2ri.activate()
    set_seed = rpy2.robjects.r("set.seed")
    set_seed(random_state)
    hm_embed = harmony.HarmonyMatrix(
            embed, meta, batch, theta, do_pca=False, verbose=False)
    pandas2ri.deactivate()
    numpy2ri.deactivate()
    hm_embed = numpy2ri.ri2py(hm_embed)
    # ===========
    if key_added:
        adata.obsm[f'{use_rep}_{key_added}'] = hm_embed
    else:
        adata.obsm[use_rep] = hm_embed


def run_bbknn(adata, batch, use_rep='X_pca', key_added='bk', **kwargs):
    if use_rep not in adata.obsm.keys():
        raise KeyError(f'{use_rep} is not a valid embedding.')
    if use_rep != 'X_pca':
        if 'X_pca' in adata.obsm.keys():
            if 'X_pca_bkup' in adata.uns.keys():
                print("overwrite existing `.obsm['X_pca_bkup']`")
            adata.obsm['X_pca_bkup'] = adata.obsm['X_pca']
        adata.obsm['X_pca'] = adata.obsm[use_rep]
    import bbknn
    if 'neighbors' in adata.uns.keys():
        if 'neighbors_bkup' in adata.uns.keys():
            print("overwrite existing `.uns['neighbors']`")
        adata.uns['neighbors_bkup'] = adata.uns['neighbors']
    bbknn.bbknn(adata, batch_key=batch, **kwargs)
    if key_added:
        adata.uns[f'neighbors_{key_added}'] = adata.uns['neighbors']
        del adata.uns['neighbors']
        if 'neighbors_bkup' in adata.uns.keys():
            adata.uns['neighbors'] = adata.uns['neighbors_bkup']
            del adata.uns['neighbors_bkup']
    if use_rep != 'X_pca':
        del adata.obsm['X_pca']
        if 'X_pca_bkup' in adata.obsm.keys():
            adata.obsm['X_pca'] = adata.obsm['X_pca_bkup']
            del adata.obsm['X_pca_bkup']


def run_phate(adata, use_rep='X', key_added=None, knn=5, decay=40, t='auto', n_pca=100, random_state=0, **kwargs):
    import phate
    if use_rep == 'X':
        data = adata.X
    elif use_rep == 'raw':
        data = adata.raw.X
    elif use_rep in adaat.obsm.keys():
        data = adata.obsm[use_rep]
    else:
        raise KeyError(f'{use_rep} not found.')
    phate_operator = phate.PHATE(
        knn=knn, decay=decay, t=t, n_pca=n_pca, random_state=random_state, **kwargs)
    kc_phate = phate_operator.fit_transform(data)
    slot_name = f'X_phate_{key_added}' if key_added else 'X_phate'
    ad_kc.obsm[slot_name] = kc_phate


def split_by_group(adata, groupby, out_type='dict'):
    if groupby not in adata.obs.columns:
        raise KeyError(f'{groupby} is not a valid obs annotation.')
    groups = sorted(list(adata.obs[groupby].unique()))
    if out_type == 'dict':
        out_adatas = {}
        for grp in groups:
            out_adatas[grp] = adata[adata.obs[groupby] == grp, :].copy()
    elif out_type == 'list':
        out_adatas = []
        for grp in groups:
            out_adatas.append(adata[adata.obs[groupby] == grp, :].copy())
    else:
        raise ValueError(f'{out_type}: unsupported type, choose from "dict" or "list".')
    return out_adatas


def show_obs_categories(ad, columns=None):
    if columns:
        columns = [k for k in columns if k in ad.obs.columns]
    else:
        columns = [k for k in ad.obs.columns if ad.obs[k].dtype.name == 'category']
    for k in columns:
        print(k)
        print(ad.obs[k].value_counts().to_dict())


def regroup(adata, groupby, regroups):
    if groupby not in adata.obs.columns:
        raise KeyError(f'{groupby} is not a valid obs annotation.')
    groups = adata.obs[groupby].astype(str)
    new_groups = groups.copy()
    for new_grp, old_grps in regroups.items():
        if isinstance(old_grps, (list, tuple)):
            for grp in old_grps:
                new_groups[groups == grp] = new_grp
        else:
            new_groups[groups == old_grps] = new_grp
    regroup_keys = [g for g in regroups.keys() if g in set(new_groups.unique())]
    new_groups = new_groups.astype('category')
    if len(regroup_keys) == len(new_groups.cat.categories):
        new_groups = new_groups.cat.reorder_categories(regroup_keys)
    return new_groups


def subsample(adata, fraction, groupby=None, min_n=0, max_n=10000, **kwargs):
    if groupby:
        if groupby not in adata.obs.columns:
            raise KeyError(f'{groupby} is not a valid obs annotation.')
        groups = adata.obs[groupby].unique()
        n_obs_per_group = {}
        for grp in groups:
            k = adata.obs[groupby] == grp
            grp_size = sum(k)
            n_obs_per_group[grp] = int(min(max_n,
                                           max(np.ceil(grp_size * fraction),
                                               min(min_n, grp_size))))
        sampled_groups = [sc.pp.subsample(adata[adata.obs[groupby] == grp, :],
                                          n_obs=n_obs_per_group[grp],
                                          copy=True,
                                          **kwargs) for grp in groups]
        sampled_obs_names = np.concatenate(
                [sg.obs_names.values for sg in sampled_groups])
        subsampled = adata[adata.obs_names.isin(sampled_obs_names), :]
    else:
        subsampled = sc.pp.subsample(adata, fraction, **kwargs, copy=True)
    return subsampled


def pseudo_bulk(
        adata, groupby, use_rep='X', highly_variable=False, FUN=np.mean):
    """Make pseudo bulk data from grouped sc data
    """
    if adata.obs[groupby].dtype.name == 'category':
        group_attr = adata.obs[groupby].values
        groups = adata.obs[groupby].cat.categories.values
    else:
        group_attr = adata.obs[groupby].astype(str).values
        groups = np.unique(group_attr)
    n_level = len(groups)
    if highly_variable:
        if isinstance(highly_variable, (list, tuple)):
            if use_rep == 'raw':
                k_hv = adata.raw.var_names.isin(highly_variable)
            else:
                k_hv = adata.var_names.isin(highly_variable)
        else:
            k_hv = adata.var['highly_variable'].values
    if use_rep == 'X':
        x = adata.X
        features = adata.var_names.values
        if highly_variable:
            x = x[:, k_hv]
            features = features[k_hv]
    elif use_rep == 'raw':
        x = adata.raw.X
        features = adata.raw.var_names.values
        if highly_variable:
            x = x[:, k_hv]
            features = features[k_hv]
    elif use_rep in adata.layers.keys():
        x = adata.layers[use_rep]
        features = adata.var_names.values
        if highly_variable:
            x = x[:, k_hv]
            features = features[k_hv]
    elif use_rep in adata.obsm.keys():
        x = adata.obsm[use_rep]
        features = np.arange(x.shape[1])
    elif (isinstance(use_rep, np.ndarray) and
            use_rep.shape[0] == adata.shape[0]):
        x = use_rep
        features = np.arange(x.shape[1])
    else:
        raise KeyError(f'{use_rep} invalid.')
    summarised = np.zeros((n_level, x.shape[1]))
    for i, grp in enumerate(groups):
        k_grp = group_attr == grp
        if sp.issparse(x):
            summarised[i] = FUN(x[k_grp, :], axis=0)
        else:
            summarised[i] = FUN(x[k_grp, :], axis=0, keepdims=True)
    return pd.DataFrame(summarised.T, columns=groups, index=features)
