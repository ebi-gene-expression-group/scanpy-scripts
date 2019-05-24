"""
Provides exported functions
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scanpy as sc

from ._read import read_10x
from ._filter import filter_anndata
from ._norm import normalize
from ._hvg import hvg
from ._neighbors import neighbors
from ._umap import umap
from ._fdg import fdg
from ._tsne import tsne
from ._louvain import louvain
from ._leiden import leiden
from ._diffexp import diffexp, diffexp_paired, extract_de_table
from ..cmd_utils import _read_obj as read_obj
from ..cmd_utils import _write_obj as write_obj


def expression_colormap():
    """Returns a nice color map for highlighting gene expression
    """
    reds = plt.cm.Reds(np.linspace(0, 1, 118))
    greys = plt.cm.Greys_r(np.linspace(0.7, 0.8, 10))
    palette = np.vstack([greys, reds])
    return LinearSegmentedColormap.from_list('expression', palette)


def cross_table(adata, x, y, normalise=None, highlight=False):
    """Make a cross table comparing two categorical annotations
    """
    x_attr = adata.obs[x]
    y_attr = adata.obs[y]
    assert not _is_numeric(x_attr.values), f'Can not operate on numerical {x}'
    assert not _is_numeric(y_attr.values), f'Can not operate on numerical {y}'
    crs_tbl = pd.crosstab(x_attr, y_attr)
    if normalise == 'x':
        x_sizes = x_attr.groupby(x_attr).size().values
        crs_tbl = (crs_tbl.T / x_sizes * 100).round(2).T
    elif normalise == 'y':
        y_sizes = x_attr.groupby(y_attr).size().values
        crs_tbl = (crs_tbl / y_sizes * 100).round(2)
    if highlight:
        return crs_tbl.style.background_gradient(cmap='viridis', axis=0)
    return crs_tbl


def _is_numeric(x):
    return x.dtype.kind in ('i', 'f')


def run_harmony(
        adata,
        batch,
        theta=2.0,
        key='X_pca',
        key_added='hm',
        tmp_dir='.harmony',
        script='harmonise.R'
):
    if not isinstance(batch, (tuple, list)):
        batch = [batch]
    if not isinstance(theta, (tuple, list)):
        theta = [theta]
    for b in batch:
        if b not in adata.obs.columns:
            raise KeyError(f'{b} is not a valid obs annotation.')
    if key not in adata.obsm.keys():
        raise KeyError(f'{key} is not a valid embedding.')
    meta = adata.obs[batch]
    embed = adata.obsm[key]

    import os
    os.makedirs(tmp_dir, exist_ok=True)
    meta_fn = os.path.join(tmp_dir, 'meta.tsv')
    embed_fn = os.path.join(tmp_dir, 'embed.tsv')
    meta.to_csv(meta_fn, sep='\t', index=True, header=True)
    pd.DataFrame(embed, index=adata.obs_names).to_csv(
        embed_fn, sep='\t', index=True, header=True)
    grouping = ','.join(batch)
    thetas = ','.join(list(map(str, theta)))
    hm_embed_fn = os.path.join(tmp_dir, 'hm_embed.tsv')

    import subprocess as sbp
    cmd = f'Rscript {script} {embed_fn} {meta_fn} {grouping} {thetas} {hm_embed_fn}'
    sbp.call(cmd.split())
    hm_embed = pd.read_csv(hm_embed_fn, header=0, index_col=0, sep='\t')
    adata.obsm[f'{key}_{key_added}'] = hm_embed.values
    os.remove(meta_fn)
    os.remove(embed_fn)
    os.remove(hm_embed_fn)


def split_by_group(adata, groupby):
    if groupby not in adata.obs.columns:
        raise KeyError(f'{groupby} is not a valid obs annotation.')
    groups = adata.obs[groupby].unique().to_list()
    adata_dict = {}
    for grp in groups:
        adata_dict[grp] = adata[adata.obs[groupby] == grp, :]
    return adata_dict


def subsample(adata, fraction, groupby=None, min_n=0, **kwargs):
    if groupby:
        if groupby not in adata.obs.columns:
            raise KeyError(f'{groupby} is not a valid obs annotation.')
        groups = adata.obs[groupby].unique()
        n_obs_per_group = {}
        for grp in groups:
            k = adata.obs[groupby] == grp
            grp_size = sum(k)
            n_obs_per_group[grp] = int(max(np.ceil(grp_size * fraction),
                                           min(min_n, grp_size)))
        sampled_groups = [sc.pp.subsample(adata[adata.obs[groupby] == grp, :],
                                          n_obs=n_obs_per_group[grp],
                                          copy=True,
                                          **kwargs) for grp in groups]
        subsampled = sampled_groups[0].concatenate(sampled_groups[1:])
        subsampled.var = adata.var.copy()
    else:
        subsampled = sc.pp.subsample(adata, fraction, **kwargs, copy=True)
    return subsampled


def pseudo_bulk(
        adata, groupby, use_rep='X', highly_variable=False, FUN=np.mean):
    """Make pseudo bulk data from grouped sc data
    """
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
    else:
        raise KeyError(f'{use_rep} not found.')
    summarised = np.zeros((n_level, x.shape[1]))
    for i, grp in enumerate(groups):
        k_grp = group_attr == grp
        summarised[i] = FUN(x[k_grp, :], axis=0, keepdims=True)
    return pd.DataFrame(summarised.T, columns=groups, index=features)


def plot_df_heatmap(
        df,
        cmap='viridis',
        title=None,
        figsize=(7, 7),
        rotation=90,
        save=None,
        **kwargs,
):
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(df, cmap=cmap, aspect='auto', **kwargs)
    if 0 < rotation < 90:
        horizontalalignment = 'right'
    plt.xticks(
        range(len(df.columns)),
        df.columns,
        rotation=rotation,
        horizontalalignment=horizontalalignment,
    )
    plt.yticks(range(len(df.index)), df.index)
    if title:
        fig.suptitle(title)
    fig.colorbar(im)
    if save:
        plt.savefig(fname=save, bbox_inches='tight', pad_inches=0.1)


def plot_qc(adata, groupby=None):
    qc_metrics = ('n_counts', 'n_genes', 'percent_mito')
    for qmt in qc_metrics:
        if qmt not in adata.obs.columns:
            raise ValueError(f'{qmt} not found.')
    if groupby:
        sc.pl.violin(adata, keys=qc_metrics, groupby=groupby, rotation=45)
    sc.pl.violin(adata, keys=qc_metrics, multi_panel=True, rotation=45)
    sc.pl.scatter(
        adata, x='n_counts', y='n_genes', color='percent_mito', alpha=0.5)
    sc.pl.scatter(adata, x='n_counts', y='n_genes', color=groupby, alpha=0.5)
    sc.pl.scatter(
        adata, x='n_counts', y='percent_mito', color=groupby, alpha=0.5)


def simple_default_pipeline(
        adata,
        qc_only=False,
        transform_x=True,
        scale=False,
        min_genes=200,
        min_cells=3,
        max_counts=25000,
        max_mito=20,
        batch=None,
        combat=False,
        combat_args=None,
        hvg_flavor='cell_ranger',
        transform_rdim=True,
        n_neighbors=15,
        n_pcs=40,
):
    if qc_only:
        adata.var['mito'] = adata.var_names.str.startswith('MT-')
        qc_tbls = sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mito'], percent_top=None)
        adata.obs['n_counts'] = qc_tbls[0]['total_counts'].values
        adata.obs['n_genes'] = qc_tbls[0]['n_genes_by_counts'].values
        adata.obs['percent_mito'] = qc_tbls[0]['pct_counts_mito'].values
        adata.var['n_cells'] = qc_tbls[1]['n_cells_by_counts'].values
        plot_qc(adata, batch)
    else:
        if transform_x:
            sc.pp.filter_cells(adata, min_genes=min_genes)
            sc.pp.filter_genes(adata, min_cells=min_cells)
            k = ((adata.obs['n_counts'] <= max_counts) &
                 (adata.obs['percent_mito'] <= max_mito))
            adata._inplace_subset_obs(k)
            if 'counts' not in adata.layers:
                adata.layers['counts'] = adata.X.copy()
            sc.pp.normalize_total(adata, target_sum=1e4, fraction=0.9)
            sc.pp.log1p(adata)
            adata.raw = adata
            if combat:
                sc.pp.combat(adata, **combat_args)
            if batch and batch in adata.obs.columns:
                by_batch = (batch, 1)
            else:
                by_batch = None
            if hvg_flavor == 'cell_ranger':
                hvg(adata, flavor='cell_ranger', n_top_genes=2000, by_batch=by_batch)
            else:
                hvg(adata, flavor='seurat', by_batch=by_batch)
            sc.pl.highly_variable_genes(adata)
            if scale:
                sc.pp.scale(adata, max_value=10)
        if transform_rdim:
            sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
            neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
            umap(adata)
            sc.tl.diffmap(adata, n_comps=15)
            leiden(adata, resolution=(0.1, 0.4, 0.7))
    return adata
