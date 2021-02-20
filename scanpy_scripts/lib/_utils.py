"""
Utility functions
"""

import numpy as np
import scipy.sparse as sp
import pandas as pd
import anndata
import scanpy as sc
from ..obj_utils import (
    _set_obsm_key,
    _restore_obsm_key,
    _backup_default_key,
    _rename_default_key,
    _delete_backup_key,
)


def lognorm_to_counts(X, n_counts=None, force=False, rounding=True):
    X_expm1 = np.expm1(X)
    if n_counts is not None:
        size_factor = n_counts / 1e4
        X_counts = (X_expm1.T * sp.csr_matrix(sp.diags(size_factor))).T
        res = np.abs(X_counts.data - np.round(X_counts.data)).sum() / X_counts.data.sum()
        if res < 1e-6:
            if rounding:
                X_counts.data = np.round(X_counts.data).astype(np.int32)
            return X_counts
        else:
            sc.logging.warn('Non-integer residuals too large, try inferring size_factor')
    x_min = np.array([X_expm1.getrow(i).data.min() for i in range(X_expm1.shape[0])])
    size_factor = 1 / x_min
    X_counts = (X_expm1.T * sp.csr_matrix(sp.diags(size_factor))).T
    res = np.abs(X_counts.data - np.round(X_counts.data)).sum() / X_counts.data.sum()
    if res < 1e-6 or force:
        if rounding:
            X_counts.data = np.round(X_counts.data).astype(np.int32)
        return X_counts
    else:
        raise ValueError('Non-integer residuals too large, failed to recover counts')


def restore_adata(adata, restore_type=['norm', 'count'], obs_cols=None, var_cols=None, obsm_keys=None):
    if not adata.raw:
        raise ValueError('adata.raw not found')

    if isinstance(restore_type, (list, tuple)):
        restore_type = restore_type[0]

    if restore_type == 'norm':
        X = adata.raw.X
    elif restore_type == 'count':
        X = lognorm_to_counts(adata.raw.X)
    else:
        raise ValueError(f'{restore_type}: invalid <restore_type>, choose between "norm" and "count"')

    obs = adata.obs[obs_cols].copy() if obs_cols else adata.obs
    var = adata.raw.var[var_cols].copy() if var_cols else adata.raw.var
    ad = anndata.AnnData(X=X, obs=obs, var=var)
    if obsm_keys:
        for ok in obsm_keys:
            if ok in adata.obsm.keys():
                ad.obsm[ok] = adata.obsm[ok]
    return ad


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
    import warnings
    from rpy2.rinterface import RRuntimeWarning
    warnings.filterwarnings("ignore", category=RRuntimeWarning)
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
    import bbknn
    _set_obsm_key(adata, 'X_pca', use_rep)
    try:
        _backup_default_key(adata.uns, 'neighbors')
        bbknn.bbknn(adata, batch_key=batch, **kwargs)
        if key_added:
            _rename_default_key(adata.uns, 'neighbors', f'neighbors_{key_added}')
        else:
            _delete_backup_key(adata.uns, 'neighbors')
    finally:
        _restore_obsm_key(adata, 'X_pca', use_rep)


def run_phate(adata, use_rep='X', key_added=None, knn=5, decay=40, t='auto', n_pca=100, random_state=0, verbose=False, **kwargs):
    import phate
    if use_rep == 'X':
        data = adata.X
    elif use_rep == 'raw':
        data = adata.raw.X
    elif use_rep in adata.obsm.keys():
        data = adata.obsm[use_rep]
    elif use_rep in adata.uns.keys():
        data = adata.uns[use_rep]['distances']
        kwargs['knn_dist'] = 'precomputed'
    else:
        raise KeyError(f'{use_rep} not found.')
    phate_operator = phate.PHATE(
        knn=knn, decay=decay, t=t, n_pca=n_pca, random_state=random_state, verbose=verbose, **kwargs)
    kc_phate = phate_operator.fit_transform(data)
    slot_name = f'X_phate_{key_added}' if key_added else 'X_phate'
    adata.obsm[slot_name] = kc_phate


def write_table(
        adata, outdir='.',
        slots=('X', 'obs', 'var', 'obsm', 'varm', 'raw.X', 'raw.var'),
        fmt='tsv', transpose=False, compression=False,
        X_dtype=None, raw_X_dtype=None, obs_columns=None, var_columns=None, obsm_keys=None, varm_keys=None, raw_var_columns=None,
):
    import os
    if not outdir:
        raise ValueError('`outdir` cannot be empty')
    outdir = outdir.rstrip('/')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    sep = ',' if fmt == 'csv' else '\t'
    suffix = '.gz' if compression else ''

    for slot in slots:
        if slot == 'X' or slot == 'raw.X':
            ad = adata.raw if slot == 'raw.X' else adata
            dtype = raw_X_dtype if slot == 'raw.X' else X_dtype
            if dtype is None:
                dtype = np.float32
            X = ad.X.T if transpose else ad.X
            if sp.issparse(X):
                X = X.toarray()
            if dtype in (int, np.int16, np.int32, np.int64):
                X = np.round(X)
            if transpose:
                df = pd.DataFrame(X, index=ad.var_names.values, columns=adata.obs_names.values, dtype=dtype)
                df.index.name = 'Gene'
            else:
                df = pd.DataFrame(X, index=adata.obs_names.values, columns=ad.var_names.values, dtype=dtype)
                df.index.name = 'Cell'
            df.to_csv(f'{outdir}/{slot}.{fmt}{suffix}', sep=sep)
        elif (slot == 'obs' or slot == 'var') and getattr(adata, slot) is not None:
            df_columns = obs_columns if slot == 'obs' else var_columns
            if df_columns is None:
                df_columns = df.columns.to_list()
            elif isinstance(df_columns, str):
                df_columns = [df_columns]
            df = getattr(adata, slot)[df_columns]
            df.index.name = 'Cell' if slot == 'obs' else 'Gene'
            df.to_csv(f'{outdir}/{slot}.{fmt}{suffix}', sep=sep)
        elif slot == 'raw.var' and adata.raw is not None and adata.raw.var is not None:
            df_columns = raw_var_columns
            if df_columns is None:
                df_columns = df.columns.to_list()
            elif isinstance(df_columns, str):
                df_columns = [df_columns]
            df = adata.raw.var[df_columns]
            df.index.name = 'Gene'
            df.to_csv(f'{outdir}/{slot}.{fmt}{suffix}', sep=sep)
        elif (slot == 'obsm' or slot == 'varm') and getattr(adata, slot) is not None:
            obj_keys = obsm_keys if slot == 'obsm' else varm_keys
            if obj_keys is None:
                obj_keys = obj.keys()
            elif isinstance(df_columns, str):
                df_columns = [df_columns]
            obj = getattr(adata, slot)
            for k in obj_keys:
                X = obj[k]
                ind = adata.obs_names.values if slot == 'obsm' else adata.var_names.values
                col = [f'{k}_{i}' for i in range(X.shape[1])]
                df =  pd.DataFrame(X, index=ind, columns=col, dtype=np.float32)
                df.index.name = 'Cell' if slot == 'obsm' else 'Gene'
                df.to_csv(f'{outdir}/{slot}.{k}.{fmt}{suffix}', sep=sep)


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
        sampled_obs_names = []
        for grp in groups:
            k = adata.obs[groupby] == grp
            grp_size = sum(k)
            ds_grp_size = int(min(
                max_n, max(np.ceil(grp_size * fraction), min(min_n, grp_size))))
            sampled_obs_names.extend(
                list(adata.obs_names[k][np.random.choice(grp_size, ds_grp_size, replace=False)]))

        subsampled = adata[adata.obs_names.isin(sampled_obs_names)].copy()
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
