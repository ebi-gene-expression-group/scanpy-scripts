import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.preprocessing import normalize
from collections import Counter
from collections import defaultdict
import anndata
import scanpy as sc
from ._diffexp import extract_de_table


def calc_marker_stats(ad, groupby, genes=None, use_rep='raw', inplace=False, partial=False):
    if ad.obs[groupby].dtype.name != 'category':
        raise ValueError('"%s" is not categorical' % groupby)
    n_grp = ad.obs[groupby].cat.categories.size
    if n_grp < 2:
        raise ValueError('"%s" must contain at least 2 categories' % groupby)
    if use_rep == 'raw':
        X = ad.raw.X
        var_names = ad.raw.var_names.values
    else:
        X = ad.X
        var_names = ad.var_names.values
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    if genes:
        v_idx = var_names.isin(genes)
        X = X[:, v_idx]
        var_names = var_names[v_idx]

    X = normalize(X, norm='max', axis=0)
    k_nonzero = X.sum(axis=0).A1 > 0
    X = X[:, np.where(k_nonzero)[0]]
    var_names = var_names[k_nonzero]

    n_var = var_names.size
    x = np.arange(n_var)

    grp_indices = {k: g.index.values for k, g in ad.obs.reset_index().groupby(groupby, sort=False)}
    
    frac_df = pd.DataFrame({k: (X[idx, :]>0).mean(axis=0).A1 for k, idx in grp_indices.items()}, index=var_names)
    mean_df = pd.DataFrame({k: X[idx, :].mean(axis=0).A1 for k, idx in grp_indices.items()}, index=var_names)

    if partial:
        stats_df = None
    else:
        frac_order = np.apply_along_axis(np.argsort, axis=1, arr=frac_df.values)
        y1 = frac_order[:, n_grp-1]
        y2 = frac_order[:, n_grp-2]
        y3 = frac_order[:, n_grp-3] if n_grp > 2 else y2
        top_frac_grps = frac_df.columns.values[y1]
        top_fracs = frac_df.values[x, y1]
        frac_diffs = top_fracs - frac_df.values[x, y2]
        max_frac_diffs = top_fracs - frac_df.values[x, y3]
 
        mean_order = np.apply_along_axis(np.argsort, axis=1, arr=mean_df.values)
        y1 = mean_order[:, n_grp-1]
        y2 = mean_order[:, n_grp-2]
        y3 = mean_order[:, n_grp-3] if n_grp > 2 else y2
        top_mean_grps = mean_df.columns.values[y1]
        top_means = mean_df.values[x, y1]
        mean_diffs = top_means - mean_df.values[x, y2]
        max_mean_diffs = top_means - mean_df.values[x, y3]

        stats_df = pd.DataFrame({
            'top_frac_group': top_frac_grps, 'top_frac': top_fracs, 'frac_diff': frac_diffs, 'max_frac_diff': max_frac_diffs,
            'top_mean_group': top_mean_grps, 'top_mean': top_means, 'mean_diff': mean_diffs, 'max_mean_diff': max_mean_diffs
        }, index=var_names)
        stats_df['top_frac_group'] = stats_df['top_frac_group'].astype('category')
        stats_df['top_frac_group'].cat.reorder_categories(list(ad.obs[groupby].cat.categories), inplace=True)

    if inplace:
        if use_rep == 'raw':
            ad.raw.varm[f'frac_{groupby}'] = frac_df
            ad.raw.varm[f'mean_{groupby}'] = mean_df
            if not partial:
                ad.raw.var = pd.concat([ad.raw.var, stats_df], axis=1)
        else:
            ad.varm[f'frac_{groupby}'] = frac_df
            ad.varm[f'mean_{groupby}'] = mean_df
            if not partial:
                ad.var = pd.concat([ad.raw.var, stats_df], axis=1)
    else:
        return frac_df, mean_df, stats_df


def filter_marker_stats(ad, use_rep='raw', min_frac_diff=0.1, min_mean_diff=0.1, max_next_frac=0.9, max_next_mean=0.95, single=False, how='or'):
    columns = ['top_frac_group', 'top_frac', 'frac_diff', 'max_frac_diff', 'top_mean_group', 'top_mean', 'mean_diff', 'max_mean_diff']
    if isinstance(ad, anndata.AnnData):
        stats_df = ad.raw.var[columns] if use_rep == 'raw' else ad.var[columns]
    elif isinstance(ad, pd.DataFrame):
        stats_df = ad[columns]
    else:
        raise ValueError('Invalid input, must be an AnnData or DataFrame')
    frac_diff = stats_df.frac_diff if single else stats_df.max_frac_diff
    mean_diff = stats_df.mean_diff if single else stats_df.max_mean_diff
    same_group = stats_df.top_frac_group == stats_df.top_mean_group
    meet_frac_requirement = (frac_diff >= min_frac_diff) & (stats_df.top_frac - frac_diff <= max_next_frac)
    meet_mean_requirement = (mean_diff >= min_mean_diff) & (stats_df.top_mean - mean_diff <= max_next_mean)
    if how == 'or':
        filtered = stats_df.loc[same_group & (meet_frac_requirement | meet_mean_requirement)]
    else:
        filtered = stats_df.loc[same_group & (meet_frac_requirement & meet_mean_requirement)]
    if single:
        filtered = filtered.sort_values(['top_frac_group', 'mean_diff', 'frac_diff'], ascending=[True, False, False])
    else:
        filtered = filtered.sort_values(['top_frac_group', 'mean_diff', 'frac_diff'], ascending=[True, False, False])
        #filtered = filtered.sort_values(['top_frac_group', 'max_mean_diff', 'max_frac_diff'], ascending=[True, False, False])
    filtered['top_frac_group'] = filtered['top_frac_group'].astype('category')
    filtered['top_frac_group'].cat.reorder_categories(list(stats_df['top_frac_group'].cat.categories), inplace=True)
    return filtered


def top_markers(df, top_n=5, groupby='top_frac_group'):
    return df.groupby(groupby).head(top_n).index.to_list()


def test_markers(ad, mks, groupby, n_genes=100, use_raw=True, **kwargs):
    genes = top_markers(mks, top_n=n_genes)
    aux_ad = anndata.AnnData(
        X=ad.raw.X if use_raw else ad.X,
        obs=ad.obs.copy(),
        var=ad.raw.var.copy() if use_raw else ad.var.copy()
    )
    aux_ad = aux_ad[:, genes].copy()
    sc.tl.rank_genes_groups(aux_ad, groupby=groupby, n_genes=n_genes, use_raw=False, **kwargs)
    de_tbl = extract_de_table(aux_ad.uns['rank_genes_groups'])
    return mks.reset_index().rename(columns={'index': 'genes', 'top_frac_group': 'cluster'}).merge(
        de_tbl[['cluster', 'genes', 'logfoldchanges', 'pvals', 'pvals_adj']], how='left'
    )
