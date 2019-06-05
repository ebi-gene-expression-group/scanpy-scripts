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
from ._pca import pca
from ._neighbors import neighbors
from ._umap import umap
from ._fdg import fdg
from ._tsne import tsne
from ._louvain import louvain
from ._leiden import leiden
from ._diffexp import diffexp, diffexp_paired, extract_de_table
from ._diffmap import diffmap
from ._dpt import dpt
from ._paga import paga
from ..cmd_utils import _read_obj as read_obj
from ..cmd_utils import _write_obj as write_obj
from ..cmd_utils import switch_layer


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
        use_rep='X_pca',
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
    if use_rep not in adata.obsm.keys():
        raise KeyError(f'{use_rep} is not a valid embedding.')
    meta = adata.obs[batch]
    embed = adata.obsm[use_rep]

    import os
    os.makedirs(tmp_dir, exist_ok=True)
    # ===========
    # from rpy2.robjects.package import importr
    # harmony = importr('harmony')
    # from rpy2.robjects import numpy2ri, pandas2ri
    # numpy2ri.activate()
    # hm_embed = harmony.HarmonyMatrix(
    #     embed, adata.obs[[batch]].reset_index(), batch, theta, do_pca=False)
    # numpy2ri.deactivate()
    # ===========
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
    hm_embed = pd.read_csv(hm_embed_fn, header=0, index_col=0, sep='\t').values
    if key_added:
        adata.obsm[f'{use_rep}_{key_added}'] = hm_embed
    else:
        adata.obsm[use_rep] = hm_embed
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
        batch=None,
        filter_params={'min_genes': 200, 'min_cells': 3, 'max_counts': 25000, 'max_mito': 20},
        norm_params={'target_sum': 1e4, 'fraction': 0.9},
        combat_args={'key': None},
        hvg_params={'flavor': 'seurat', 'by_batch': None},
        scale_params={'max_value': 10},
        pca_params={'n_comps': 50, 'svd_solver': 'arpack', 'use_highly_variable': True},
        harmony_params={'batch': None, 'theta': 2.0, 'script': 'harmonise.R'},
        nb_params={'n_neighbors': 15, 'n_pcs': 40},
        umap_params={},
        tsne_params={},
        diffmap_params={'n_comps': 15},
        leiden_params={
            'resolution': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]},
        fdg_params={'layout': 'fa'},
):
    """
    Scanpy pipeline
    """
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
        if filter_params is not None and isinstance(filter_params, dict):
            if 'min_genes' in filter_params:
                sc.pp.filter_cells(adata, min_genes=filter_params['min_genes'])
            if 'min_cells' in filter_params:
                sc.pp.filter_genes(adata, min_cells=filter_params['min_cells'])
            if 'min_counts' in filter_params:
                k = adata.obs['n_counts'] >= filter_params['min_counts']
                adata._inplace_subset_obs(k)
            if 'max_counts' in filter_params:
                k = adata.obs['n_counts'] <= filter_params['max_counts']
                adata._inplace_subset_obs(k)
            if 'max_mito' in filter_params:
                k = adata.obs['percent_mito'] <= filter_params['max_mito']
                adata._inplace_subset_obs(k)
            if 'counts' not in adata.layers.keys():
                adata.layers['counts'] = adata.X
        if norm_params is not None and isinstance(norm_params, dict):
            if 'counts' in adata.layers.keys():
                adata.X = adata.layers['counts']
            sc.pp.normalize_total(adata, **norm_params)
            sc.pp.log1p(adata)
            adata.raw = adata
        if (combat_args is not None and (
                isinstance(combat_args, dict) and
                combat_args.get('key', None) and
                combat_args['key'] in adata.obs.keys())):
            adata.layers['X'] = adata.X
            adata.X = adata.raw.X
            sc.pp.combat(adata, **combat_args)
        if hvg_params is not None and isinstance(hvg_params, dict):
            hvg(adata, **hvg_params)
        if scale_params is not None and isinstance(scale_params, dict):
            sc.pp.scale(adata, **scale_params)
        if pca_params is not None and isinstance(pca_params, dict):
            pca(adata, **pca_params)
        if (harmony_params is not None and (
                isinstance(harmony_params, dict) and
                harmony_params.get('batch', None) and
                harmony_params['batch'] in adata.obs.keys())):
            run_harmony(adata, **harmony_params)
        if nb_params is not None and isinstance(nb_params, dict):
            neighbors(adata, **nb_params)
        if umap_params is not None and isinstance(umap_params, dict):
            umap(adata, **umap_params)
        if tsne_params is not None and isinstance(tsne_params, dict):
            tsne(adata, **tsne_params)
        if diffmap_params is not None and isinstance(diffmap_params, dict):
            diffmap(adata, **diffmap_params)
        if leiden_params is not None and isinstance(leiden_params, dict):
            leiden(adata, **leiden_params)
        if fdg_params is not None and isinstance(fdg_params, dict):
            fdg(adata, **fdg_params)
    return adata
