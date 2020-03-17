"""
Provides exported functions
"""

import numpy as np
import scanpy as sc
from scanpy.plotting._tools.scatterplots import plot_scatter
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt

from ._read import read_10x, read_10x_atac
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
from ._paga import paga, plot_paga
from ._doublets import run_scrublet, test_outlier
from ._plot import (
    expression_colormap,
    clear_colors,
    cross_table,
    set_figsize,
    plot_df_heatmap,
    plot_qc,
    plot_metric_by_rank,
    plot_embedding,
    plot_diffexp,
    highlight,
    dotplot2,
)
from ._utils import (
    run_harmony,
    run_bbknn,
    run_phate,
    split_by_group,
    regroup,
    subsample,
    pseudo_bulk,
    show_obs_categories
)
from ..cmd_utils import _read_obj as read_obj
from ..cmd_utils import _write_obj as write_obj
from ..cmd_utils import switch_layer


def calculate_qc(adata):
    if 'mito' not in adata.var.columns:
        adata.var['mito'] = adata.var_names.str.startswith('MT-')
    if 'ribo' not in adata.var.columns:
        adata.var['ribo'] = adata.var_names.str.startswith('RPL') | adata.var_names.str.startswith('RPS')
    if 'hb' not in adata.var.columns:
        adata.var['hb'] = adata.var_names.str.startswith('HB')
    qc_tbls = sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mito', 'ribo', 'hb'], percent_top=[50])
    adata.obs['n_counts'] = qc_tbls[0]['total_counts'].values
    adata.obs['log1p_n_counts'] = np.log1p(adata.obs['n_counts'])
    adata.obs['n_genes'] = qc_tbls[0]['n_genes_by_counts'].values
    adata.obs['log1p_n_genes'] = np.log1p(adata.obs['n_genes'])
    adata.obs['percent_mito'] = qc_tbls[0]['pct_counts_mito'].values
    adata.obs['percent_ribo'] = qc_tbls[0]['pct_counts_ribo'].values
    adata.obs['percent_hb'] = qc_tbls[0]['pct_counts_hb'].values
    adata.obs['percent_top50'] = qc_tbls[0]['pct_counts_in_top_50_genes'].values
    adata.var['n_counts'] = qc_tbls[1]['total_counts'].values
    adata.var['n_cells'] = qc_tbls[1]['n_cells_by_counts'].values


def _scale_factor(x):
    xmin = np.min(x)
    xmax = np.max(x)
    return 5.0 / (xmax - xmin)


def fit_gaussian(x, n=10, threshold=0.05, xmin=None, xmax=None, plot=False, nbins=500, hist_bins=100):
    xmin = x.min() if xmin is None else xmin
    xmax = x.max() if xmax is None else xmax
    gmm = GaussianMixture(n_components=n, random_state=0)
    x_fit = x[(x>=xmin) & (x<=xmax)]
    f = _scale_factor(x_fit)
    x_fit = x_fit * f
    gmm.fit(x_fit.reshape(-1, 1))
    while not gmm.converged_:
        gmm.fit(x_fit.reshape(-1, 1), warm_start=True)
    x0 = np.linspace(x.min(), x.max(), num=nbins)
    y_pdf = np.zeros((n, nbins))
    y_cdf = np.zeros((n, nbins))
    for i in range(n):
        y_pdf[i] = norm.pdf(x0 * f, loc=gmm.means_[i, 0], scale=n  *gmm.covariances_[i, 0, 0]) * gmm.weights_[i]
        y_cdf[i] = norm.cdf(x0 * f, loc=gmm.means_[i, 0], scale=n  *gmm.covariances_[i, 0, 0]) * gmm.weights_[i]
    y0 = y_pdf.sum(axis=0)
    y1 = y_cdf.sum(axis=0)
    x_peak = x0[np.argmax(y0)]
    try:
        x_left = x0[(y0 < threshold) & (x0 < x_peak)].max()
    except:
        sc.logging.warn('Failed to find lower bound, using min value instead.')
        x_left = x0.min()
    try:
        x_right = x0[(y0 < threshold) & (x0 > x_peak)].min()
    except:
        sc.logging.warn('Failed to find upper bound, using max value instead.')
        x_right = x0.max()
    if plot:
        fig, ax1 = plt.subplots()
        _ = ax1.hist(x, bins=hist_bins)
        ax2 = ax1.twinx()
        ax2.plot(x0, y0, c='k')
        ax2.hlines(y=threshold, xmin=x.min(), xmax=x.max(), linewidth=1, linestyle='dotted')
        ax2.vlines(x=[xmin, xmax], ymin=y0.min(), ymax=y0.max(), linewidth=1, linestyle='dashed')
        if not np.isnan(x_left):
            ax2.vlines(x=x_left, ymin=y0.min(), ymax=y0.max(), linewidth=1)
        if not np.isnan(x_right):
            ax2.vlines(x=x_right, ymin=y0.min(), ymax=y0.max(), linewidth=1)
    return x_left, x_right, gmm


def auto_qc_filter(
        adata,
        metrics=['n_counts', 'n_genes', 'percent_mito', 'percent_ribo', 'percent_hb', 'percent_top50'],
        min_count=1000, min_gene=100, min_mito=0.01, max_mito=20, min_ribo=0, max_ribo=100
):
    k_pass = np.ones(adata.n_obs).astype(bool)

    if 'n_counts' in metrics:
        x_low, x_high, _ = fit_gaussian(adata.obs['log1p_n_counts'].values, xmin=np.log1p(min_count))
        min_count = int(np.expm1(x_low))
        max_count = int(np.expm1(x_high))
        k_count = (adata.obs['n_counts'] >= min_count) & (adata.obs['n_counts'] <= max_count)
        sc.logging.warn(f'n_counts: [{min_count}, {max_count}], {k_count.sum()} pass')
        k_pass = k_pass & k_count

    if 'n_genes' in metrics:
        x_low, x_high, _ = fit_gaussian(adata.obs['log1p_n_genes'].values, xmin=np.log1p(min_gene))
        min_gene = int(np.expm1(x_low))
        max_gene = int(np.expm1(x_high))
        k_gene = (adata.obs['n_genes'] >= min_gene) & (adata.obs['n_genes'] <= max_gene)
        sc.logging.warn(f'n_genes: [{min_gene}, {max_gene}], {k_gene.sum()} pass')
        k_pass = k_pass & k_gene

    if 'percent_mito' in metrics:
        max_mito = 20
        if (adata.obs['percent_mito'].values > 0).sum() > 0:
            x_low, x_high, _ = fit_gaussian(np.log1p(adata.obs['percent_mito'].values), xmin=np.log1p(min_mito), xmax=np.log1p(max_mito))
            max_mito = np.expm1(x_high)
        k_mito = (adata.obs['percent_mito'] <= max_mito)
        sc.logging.warn(f'percent_mito: [0, {max_mito}], {k_mito.sum()} pass')
        k_pass = k_pass & k_mito

    if 'percent_ribo' in metrics:
        x_low, x_high, _ = fit_gaussian(np.log1p(adata.obs['percent_ribo'].values), xmin=np.log1p(min_ribo), xmax=np.log1p(max_ribo))
        min_ribo = np.expm1(x_low)
        max_ribo = np.expm1(x_high)
        k_ribo = (adata.obs['percent_ribo'] >= min_ribo) & (adata.obs['percent_ribo'] <= max_ribo)
        sc.logging.warn(f'percent_ribo: [{min_ribo}, {max_ribo}], {k_ribo.sum()} pass')
        k_pass = k_pass & k_ribo

    if 'percent_hb' in metrics:
        max_hb = 1.0
        k_hb = adata.obs['percent_hb'] <= max_hb
        sc.logging.warn(f'percent_hb: [0, 10], {k_hb.sum()} pass')
        k_pass = k_pass & k_hb

    if 'percent_top50' in metrics:
        x_low, x_high, _ = fit_gaussian(adata.obs['percent_top50'].values)
        max_top50 = x_high
        min_top50 = x_low
        k_top50 = (adata.obs['percent_top50'] <= max_top50) & (adata.obs['percent_top50'] >= min_top50)
        sc.logging.warn(f'percent_top50: [{min_top50}, {max_top50}], {k_top50.sum()} pass')
        k_pass = k_pass & k_top50

    sc.logging.warn(f'{k_pass.sum()} pass')
    return k_pass


def get_good_sized_batch(batches, min_size=10):
    x = batches.value_counts()
    return x.index[x >= min_size].to_list()


def simple_default_pipeline(adata, batch=None, filter_only=False, post_norm_only=False, post_filter_only=False, post_pca_only=False, filter_kw={}):
    if not (post_filter_only or post_norm_only or post_pca_only):
        calculate_qc(adata)
        k_cell = auto_qc_filter(adata, **filter_kw)
        if batch and batch in adata.obs.columns:
            batches = get_good_sized_batch(adata.obs.loc[k_cell, batch])
            k_cell = k_cell & adata.obs[batch].isin(batches)
        adata = adata[k_cell, :].copy()
        adata.var['n_counts'] = adata.X.sum(axis=0).A1
        adata.var['n_genes'] = (adata.X > 0).sum(axis=0).A1
        k_gene = adata.var['n_genes'] >= 3
        adata = adata[:, k_gene].copy()
        if filter_only:
            return adata

    if not (post_norm_only or post_pca_only):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    if not post_pca_only:
        adata.raw = adata
        hvg(adata, flavor='seurat')
        zero_center = (adata.n_obs <= 20000)
        sc.pp.scale(adata, zero_center=zero_center, max_value=10)
        sc.pp.pca(adata, n_comps=30, zero_center=zero_center, svd_solver='arpack', use_highly_variable=True)

    if batch and batch in adata.obs.columns:
        run_harmony(adata, batch=batch)
        neighbors(adata, use_rep='X_pca_hm', key_added='hm', n_pcs=20, n_neighbors=15)
        umap(adata, use_graph='neighbors_hm', key_added='hm')
        leiden(adata, use_graph='neighbors_hm', resolution=[0.1, 0.3, 0.5, 0.7, 0.9], key_added='hm')
    else:
        sc.pp.neighbors(adata, use_rep='X_pca', n_pcs=20, n_neighbors=15)
        sc.tl.umap(adata)
        leiden(adata, resolution=[0.1, 0.3, 0.5, 0.7, 0.9])

    return adata


def custom_pipeline(
        adata,
        qc_only=False,
        plot=True,
        batch=None,
        filter_params={'min_genes': 200, 'min_cells': 3, 'max_counts': 25000, 'max_mito': 20, 'min_mito': 0},
        norm_params={'target_sum': 1e4, 'fraction': 0.9},
        combat_args={'key': None},
        hvg_params={'flavor': 'seurat', 'by_batch': None},
        scale_params={'max_value': 10},
        pca_params={'n_comps': 50, 'svd_solver': 'arpack', 'use_highly_variable': True},
        harmony_params={'batch': None, 'theta': 2.0},
        nb_params={'n_neighbors': 15, 'n_pcs': 20},
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
        calculate_qc(adata)
        if plot:
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
            if 'min_mito' in filter_params:
                k = adata.obs['percent_mito'] >= filter_params['min_mito']
                adata._inplace_subset_obs(k)
            if 'max_mito' in filter_params:
                k = adata.obs['percent_mito'] <= filter_params['max_mito']
                adata._inplace_subset_obs(k)
            if 'min_ribo' in filter_params:
                k = adata.obs['percent_ribo'] >= filter_params['min_ribo']
                adata._inplace_subset_obs(k)
            if 'max_ribo' in filter_params:
                k = adata.obs['percent_ribo'] <= filter_params['max_ribo']
                adata._inplace_subset_obs(k)
            if 'min_hb' in filter_params:
                k = adata.obs['percent_hb'] >= filter_params['min_hb']
                adata._inplace_subset_obs(k)
            if 'max_hb' in filter_params:
                k = adata.obs['percent_hb'] <= filter_params['max_hb']
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
                harmony_params.get('batch', None))):
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
