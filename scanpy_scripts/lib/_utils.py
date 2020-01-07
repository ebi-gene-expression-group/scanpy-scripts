"""
Utility functions
"""

import numpy as np
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


def run_seurat_integration(
        adata,
):
    pass


def split_by_group(adata, groupby):
    if groupby not in adata.obs.columns:
        raise KeyError(f'{groupby} is not a valid obs annotation.')
    groups = adata.obs[groupby].unique().to_list()
    adata_dict = {}
    for grp in groups:
        adata_dict[grp] = adata[adata.obs[groupby] == grp, :]
    return adata_dict


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
    return new_groups.astype('category')


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
            k_hv = adata.var_names.isin(highly_variable)
        else:
            k_hv = adata.var['highly_variable'].values
    if use_rep == 'X':
        x = adata.X
        features = adata.var_names.values
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
        summarised[i] = FUN(x[k_grp, :], axis=0, keepdims=True)
    return pd.DataFrame(summarised.T, columns=groups, index=features)


def LR_annotate(
        adata,
        train_label,
        train_x,
        use_rep='X',
        highly_variable=True,
        subset=None,
        penalty='l1',
        sparcity=0.2,
        return_label=False,
):
    """Annotate cells using logistic regression
    """
    if subset is not None:
        train_label = train_label[subset]
        train_x = train_x[subset, :]
    from sklearn.linear_model import LogisticRegression
    lr = LogisticRegression(penalty=penalty, C=sparcity)
    lr.fit(train_x, train_label)
    if use_rep == 'X':
        if highly_variable:
            predict_x = adata.X[:, adata.var.highly_variable]
        else:
            predict_x = adata.X
    elif use_rep in adata.layers.keys():
        predict_x = adata.layers[use_rep]
    elif use_rep in adata.obsm.keys():
        predict_x = adata.obsm[use_rep]
    if return_label:
        return lr.predict(predict_x)
    else:
        return lr.predict_proba(predict_x)


def annotate(adata, groupby, annot, annotation_matrix):
    """Annotate a clustering based on a matrix of values

    + groupby: the clustering to be annotated
    + annot : the annotationo to be aligned/transferred
    + annotation_matrix: a matrix with shape (n,m) where n equals the number
                         of clusters and m equals the number of annotation
                         categories.
    """
    if groupby not in adata.obs.columns:
        raise KeyError(f'"{groupby}" not found.')
    if annot not in adata.obs.columns:
        raise KeyError(f'"{annot}" not found.')
    groupby_annot = f'{groupby}_annot'
    clustering = adata.obs[groupby]
    annotation = adata.obs[annot]
    k_nan = annotation == 'nan'
    cell_types = np.unique(annotation[~k_nan].values)
    n_celltype = cell_types.shape[0]
    n_cluster = np.unique(clustering).shape[0]
    n_cell = clustering.shape[0]

    def dedup(names):
        names_copy = names.copy()
        count = {}
        for i, name in enumerate(names):
            n = count.get(name, 0) + 1
            count[name] = n
            if n > 1 or name in names[(i+1):]:
                names_copy[i] = f'{name}{n}'
            else:
                names_copy[i] = name
        return names_copy

    if (annotation_matrix.shape[0] == n_cluster and
            annotation_matrix.shape[1] == n_celltype):
        k_max = np.argmax(annotation_matrix, axis=1)
        cluster_annot = cell_types[k_max]
        adata.obs[groupby_annot] = adata.obs[groupby].cat.rename_categories(
                dedup(cluster_annot))
    else:
        raise ValueError(f'[annotation_matrix] not in compatible shape')
