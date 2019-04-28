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
from ._tsne import tsne
from ._louvain import louvain
from ._leiden import leiden
from ._diffexp import diffexp
from ..cmd_utils import _read_obj as read_obj
from ..cmd_utils import _write_obj as write_obj


def expression_colormap():
    """Returns a nice color map for highlighting gene expression
    """
    reds = plt.cm.Reds(np.linspace(0, 1, 128))
    greys = plt.cm.Greys_r(np.linspace(0.7, 0.8, 20))
    palette = np.vstack([greys, reds])
    return LinearSegmentedColormap.from_list('expression', palette)


def cross_table(adata, x, y, normalise=None, highlight=False):
    """Make a cross table comparing two categorical annotations
    """
    x_attr = adata.obs[x]
    y_attr = adata.obs[y]
    assert not _is_numeric(x_attr.values), "Can not operate on numeric variable {}".format(x)
    assert not _is_numeric(y_attr.values), "Can not operate on numeric variable {}".format(y)
    crs_tbl = pd.crosstab(x, y)
    if normalise == 'x':
        x_sizes = x_attr.groupby(x_attr).size().values
        crs_tbl = (crs_tbl.T / x_sizes * 100).round(2).T
    elif normalise == 'y':
        y_sizes = x_attr.groupby(y_attr).size().values
        crs_tbl = (crs_tbl / y_sizes * 100).round(2)
    if highlight:
        return crs_tbl.style.background_gradient(cmap='viridis', axis=0)
    return crs_tbl


def run_harmony():
    # import subprocess as sbp
    pass


def _is_numeric(x):
    return x.dtype.kind in ('i', 'f')


def pseudo_bulk(adata, groupby, FUN=np.mean):
    """Make pseudo bulk data from grouped sc data
    """
    group_attr = adata.obs[groupby].astype(str).values
    groups = np.unique(group_attr)
    n_level = len(groups)
    summarised = np.zeros((n_level, adata.X.shape[1]))
    for i, grp in enumerate(groups):
        k_grp = group_attr == grp
        summarised[i] = FUN(adata.X[k_grp, :], axis=0, keepdims=True)
    return summarised


def plot_qc(adata):
    qc_metrics = ('n_counts', 'n_genes', 'percent_mito')
    for qmt in qc_metrics:
        if qmt not in adata.obs.columns:
            raise ValueError('{} not found.'.format(qmt))
    sc.pl.violin(adata, keys=qc_metrics, groupby='Sample', rotation=45)
    sc.pl.violin(adata, keys=qc_metrics, multi_panel=True, rotation=45)
    sc.pl.scatter(adata, x='n_counts', y='n_genes', color='percent_mito', alpha=0.5)
    sc.pl.scatter(adata, x='n_counts', y='n_genes', color='Sample', alpha=0.5)
    sc.pl.scatter(adata, x='n_counts', y='percent_mito', color='Sample', alpha=0.5)


def simple_default_pipeline(
        h5ad_fn,
        qc_only=False,
        min_genes=200,
        min_cells=3,
        max_counts=25000,
        max_mito=20,
        batch=None,
        n_neighbors=15,
        n_pcs=40,
):
    adata = read_obj(h5ad_fn)
    if qc_only:
        adata.var['mito'] = adata.var_name.str.startswith('MT-')
        qc_tbls = sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mito'], percent_top=None)
        adata.obs['n_counts'] = qc_tbls[0]['total_counts'].values
        adata.obs['n_genes'] = qc_tbls[0]['n_genes_by_counts'].values
        adata.obs['percent_mito'] = qc_tbls[0]['pct_counts_mito'].values
        adata.var['n_cells'] = qc_tbls[1]['n_cells_by_counts'].values
        plot_qc(adata)
    else:
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        k = ((adata.obs['n_counts'] <= max_counts) &
             (adata.obs['percent_mito'] <= max_mito))
        adata = adata[k, :]
        adata.layers['counts'] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4, fraction=0.9)
        sc.pp.log1p(adata)
        adata.raw = adata
        if batch and batch in adata.obs.columns:
            sc.pp.combat(adata, key=batch)
        hvg(adata, flavor='cell_ranger')
        sc.pl.highly_variable_genes(adata)
        sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
        neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        umap(adata)
        sc.tl.draw_graph(adata)
        leiden(adata, resolution=(0.1, 0.4, 0.7))
    return adata
