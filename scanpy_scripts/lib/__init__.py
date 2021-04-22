"""
Provides exported functions
"""

import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt

if sc.__version__.startswith('1.4'):
    from scanpy.plotting._tools.scatterplots import plot_scatter
else:
    plot_scatter = sc.pl.embedding

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
from ._leiden import leiden, leiden_shredding
from ._diffexp import diffexp, diffexp_paired, extract_de_table
from ._diffmap import diffmap
from ._dpt import dpt
from ._paga import paga, plot_paga
from ._doublets import run_scrublet
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
    dotplot_combined_coexpression,
    plot_genes,
    plot_markers,
)
from ._utils import (
    lognorm_to_counts,
    restore_adata,
    run_harmony,
    run_bbknn,
    run_phate,
    split_by_group,
    regroup,
    subsample,
    pseudo_bulk,
    show_obs_categories,
    write_table,
)
from ._annot import (
    LR_train,
    LR_predict,
    annotate,
)
from ._markers import (
    calc_marker_stats,
    filter_marker_stats,
    test_markers,
)
from ..cmd_utils import _read_obj as read_obj
from ..cmd_utils import _write_obj as write_obj
from ..cmd_utils import switch_layer


class QcLowPassError(ValueError):
    pass


def calculate_qc(adata, mt_prefix='MT-', rb_prefix=('RPL', 'RPS'), hb_prefix='HB', ercc_prefix=None, extra_qc_vars=[], flag_only=False):
    if 'mito' not in adata.var.columns and mt_prefix:
        adata.var['mito'] = adata.var_names.str.startswith(mt_prefix)
    if 'ribo' not in adata.var.columns and rb_prefix and isinstance(rb_prefix, (list, tuple)) and len(rb_prefix) == 2:
        adata.var['ribo'] = adata.var_names.str.startswith(rb_prefix[0]) | adata.var_names.str.startswith(rb_prefix[1])
    if 'hb' not in adata.var.columns and hb_prefix:
        adata.var['hb'] = adata.var_names.str.startswith(hb_prefix)
    if 'ercc' not in adata.var.columns and ercc_prefix:
        adata.var['ercc'] = adata.var_names.str.startswith(ercc_prefix)
    if flag_only:
        return
    qc_vars = []
    for metric in ('mito', 'ribo', 'hb', 'ercc'):
        if metric in adata.var.columns:
            qc_vars.append(metric)
    for qv in extra_qc_vars:
        if qv in adata.var.keys():
            qc_vars.append(qv)
    qc_tbls = sc.pp.calculate_qc_metrics(
        adata, qc_vars=qc_vars, percent_top=[50])
    adata.obs['n_counts'] = qc_tbls[0]['total_counts'].values
    adata.obs['log1p_n_counts'] = np.log1p(adata.obs['n_counts'])
    adata.obs['n_genes'] = qc_tbls[0]['n_genes_by_counts'].values
    adata.obs['log1p_n_genes'] = np.log1p(adata.obs['n_genes'])
    for metric in qc_vars:
        adata.obs[f'percent_{metric}'] = qc_tbls[0][f'pct_counts_{metric}'].values
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
        min_count=1000, min_gene=100, min_mito=0.01, max_mito=20, min_ribo=0, max_ribo=100, max_hb=1.0, min_top50=0, max_top50=100,
        min_pass_rate = 0.6, force=False
):
    k_pass = np.ones(adata.n_obs).astype(bool)

    if 'n_counts' in metrics:
        try:
            x_low, x_high, _ = fit_gaussian(adata.obs['log1p_n_counts'].values, xmin=np.log1p(min_count))
        except ValueError:
            x_low, x_high = np.log1p(min_count), adata.obs['log1p_n_counts'].max()
        min_count = int(np.expm1(x_low))
        max_count = int(np.expm1(x_high))
        k_count = (adata.obs['n_counts'] >= min_count) & (adata.obs['n_counts'] <= max_count)
        sc.logging.warn(f'n_counts: [{min_count}, {max_count}], {k_count.sum()} pass')
        if k_count.sum() < adata.n_obs * min_pass_rate and not force:
            raise QcLowPassError('n_counts')
        k_pass = k_pass & k_count

    if 'n_genes' in metrics:
        try:
            x_low, x_high, _ = fit_gaussian(adata.obs['log1p_n_genes'].values, xmin=np.log1p(min_gene))
        except ValueError:
            x_low, x_high = np.log1p(min_gene), adata.obs['log1p_n_genes'].max()
        min_gene = int(np.expm1(x_low))
        max_gene = int(np.expm1(x_high))
        k_gene = (adata.obs['n_genes'] >= min_gene) & (adata.obs['n_genes'] <= max_gene)
        sc.logging.warn(f'n_genes: [{min_gene}, {max_gene}], {k_gene.sum()} pass')
        if k_gene.sum() < adata.n_obs * min_pass_rate and not force:
            raise QcLowPassError('n_genes')
        k_pass = k_pass & k_gene

    if 'percent_mito' in metrics:
        max_mito = max_mito
        if (adata.obs['percent_mito'].values > 0).sum() > 0:
            x_low, x_high, _ = fit_gaussian(np.log1p(adata.obs['percent_mito'].values), xmin=np.log1p(min_mito), xmax=np.log1p(max_mito))
            max_mito = np.expm1(x_high)
        k_mito = (adata.obs['percent_mito'] <= max_mito)
        sc.logging.warn(f'percent_mito: [0, {max_mito}], {k_mito.sum()} pass')
        if k_mito.sum() < adata.n_obs * min_pass_rate and not force:
            raise QcLowPassError('percent_mito')
        k_pass = k_pass & k_mito

    if 'percent_ribo' in metrics:
        x_low, x_high, _ = fit_gaussian(np.log1p(adata.obs['percent_ribo'].values), xmin=np.log1p(min_ribo), xmax=np.log1p(max_ribo))
        min_ribo = np.expm1(x_low)
        max_ribo = np.expm1(x_high)
        k_ribo = (adata.obs['percent_ribo'] >= min_ribo) & (adata.obs['percent_ribo'] <= max_ribo)
        sc.logging.warn(f'percent_ribo: [{min_ribo}, {max_ribo}], {k_ribo.sum()} pass')
        if k_mito.sum() < adata.n_obs * min_pass_rate and not force:
            raise QcLowPassError('percent_ribo')
        k_pass = k_pass & k_ribo

    if 'percent_hb' in metrics:
        k_hb = adata.obs['percent_hb'] <= max_hb
        sc.logging.warn(f'percent_hb: [0, {max_hb}], {k_hb.sum()} pass')
        k_pass = k_pass & k_hb
        if k_hb.sum() < adata.n_obs * min_pass_rate and not force:
            raise QcLowPassError('percent_hb')

    if 'percent_top50' in metrics:
        x_low, x_high, _ = fit_gaussian(adata.obs['percent_top50'].values, xmin=min_top50, xmax=max_top50)
        max_top50 = x_high
        min_top50 = x_low
        k_top50 = (adata.obs['percent_top50'] <= max_top50) & (adata.obs['percent_top50'] >= min_top50)
        sc.logging.warn(f'percent_top50: [{min_top50}, {max_top50}], {k_top50.sum()} pass')
        if k_top50.sum() < adata.n_obs * min_pass_rate and not force:
            raise QcLowPassError('percent_top50')
        k_pass = k_pass & k_top50

    sc.logging.warn(f'{k_pass.sum()} pass')
    return k_pass


def get_good_sized_batch(batches, min_size=10):
    x = batches.value_counts()
    return x.index[x >= min_size].to_list()


def remove_genes(adata, var_flags):
    if isinstance(var_flags, (tuple, list)):
        pass
    elif isinstance(var_flags, str):
        var_flags = [var_flags]
    else:
        var_flags = []
    if var_flags:
        mask = np.zeros(adata.n_vars).astype(bool)
        for vf in var_flags:
            mask = mask | adata.var[vf]
        adata = adata[:, ~mask].copy()
    return adata


def simple_default_pipeline(
        adata, batch=None, filter_only=False, post_norm_only=False, post_filter_only=False, post_pca_only=False, do_clustering=True,
        zero_center=None, batch_method='harmony', do_combat=False, random_state=0, clustering_resolution=[0.1, 0.3, 0.5, 0.7, 0.9],
        filter_kw={}, hvg_kw={}, pca_kw={}, nb_kw={}, umap_kw={}, hm_kw={}, bk_kw={}):
    if not (post_filter_only or post_norm_only or post_pca_only):
        if not np.all(pd.Series(
                    ['n_counts', 'n_genes', 'percent_mito', 'percent_ribo', 'percent_hb', 'percent_top50']
                ).isin(list(adata.obs.columns))):
            calculate_qc(adata)
        if (adata.obs['n_counts'] == 0).sum() > 0:
            adata = adata[adata.obs['n_counts']>0].copy()
        k_cell = auto_qc_filter(adata, **filter_kw)
        if batch:
            if isinstance(batch, (list, tuple)):
                for b in batch:
                    if b in adata.obs.columns:
                        batches = get_good_sized_batch(adata.obs.loc[k_cell, b])
                        k_cell = k_cell & adata.obs[b].isin(batches)
            elif isinstance(batch, str):
                batches = get_good_sized_batch(adata.obs.loc[k_cell, batch])
                k_cell = k_cell & adata.obs[batch].isin(batches)
        adata = adata[k_cell, :].copy()
        if filter_only:
            return adata

    if not (post_norm_only or post_pca_only):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    if not post_pca_only:
        if adata.raw is None:
            adata.raw = adata
        else:
            adata.X = adata.raw.X[:, adata.raw.var_names.isin(adata.var_names)].copy()
        if 'n_counts' not in adata.var.keys():
            n_counts = np.expm1(adata.X).sum(axis=0)
            adata.var['n_counts'] = n_counts.A1 if sp.issparse(adata.X) else n_counts
        n_cells = (adata.X > 0).sum(axis=0)
        adata.var['n_cells'] = n_cells.A1 if sp.issparse(adata.X) else n_cells
        k_gene = adata.var['n_cells'] >= 3
        if (~k_gene).sum() > 0:
            adata = adata[:, k_gene].copy()
        if batch and do_combat:
            sc.pp.combat(adata, key=batch)
        hvg(adata, flavor='seurat', **hvg_kw)
        zero_center = (adata.n_obs <= 20000) if zero_center is None else zero_center
        sc.pp.scale(adata, zero_center=zero_center, max_value=10)
        n_comps = min(30, adata.n_obs-1, adata.var.highly_variable.sum()-1)
        if n_comps < 2:
            raise ValueError('n_obs or n_vars too small for pca calculation')
        pca(adata, n_comps=n_comps, zero_center=zero_center, svd_solver='arpack', use_highly_variable=True, **pca_kw)

    n_neighbors = nb_kw.pop('n_neighbors', 15)
    n_pcs = nb_kw.pop('n_pcs', 20)
    if batch:
        if batch_method == 'bbknn':
            key_added = bk_kw.pop('key_added', 'bk')
            run_bbknn(adata, batch=batch, key_added=key_added, **bk_kw)
            umap(adata, use_graph=f'neighbors_{key_added}', key_added=key_added, random_state=random_state, **umap_kw)
            if do_clustering:
                leiden(adata, use_graph=f'neighbors_{key_added}', resolution=clustering_resolution, key_added=key_added)
        else:
            key_added = hm_kw.pop('key_added', 'hm')
            run_harmony(adata, batch=batch, key_added=key_added, **hm_kw)
            neighbors(adata, use_rep=f'X_pca_{key_added}', key_added=key_added, n_pcs=n_pcs, n_neighbors=n_neighbors)
            umap(adata, use_graph=f'neighbors_{key_added}', key_added=key_added, random_state=random_state, **umap_kw)
            if do_clustering:
                leiden(adata, use_graph=f'neighbors_{key_added}', resolution=clustering_resolution, key_added=key_added)
    else:
        sc.pp.neighbors(adata, use_rep='X_pca', n_pcs=n_pcs, n_neighbors=n_neighbors)
        sc.tl.umap(adata, random_state=random_state, **umap_kw)
        if do_clustering:
            leiden(adata, resolution=clustering_resolution)

    return adata


def subcluster(adata, groupby, groups, res, new_key, ad_aux=None, **kwargs):
    kwargs['post_norm_only'] = True
    kwargs['do_clustering'] = False
    if isinstance(res, (list, tuple)):
        res = res[0]
    if 'batch' in kwargs:
        #key_added = 'bk' if kwargs['batch_method'] == 'bbknn' else 'hm'
        graph = 'neighbors_bk' if kwargs['batch_method'] == 'bbknn' else 'neighbors_hm'
    else:
        #key_added = None
        graph = 'neighbors'
    k_groups = adata.obs[groupby].isin(groups)
    if ad_aux is None:
        ad = adata[k_groups].copy()
        ad_aux = simple_default_pipeline(ad, **kwargs)
        return_ad = True
    else:
        return_ad = False
    leiden(ad_aux, use_graph=graph, resolution=res, key_added='aux')
    adata.obs[new_key] = adata.obs[groupby].astype(str)
    adata.obs.loc[k_groups, new_key] = '_'.join(groups) + ',' + ad_aux.obs['leiden_aux'].astype(str)
    if return_ad:
        return ad_aux


def save_pipeline_object(
    ad,
    out_prefix=None,
    batch_method=None,
    obs_keys=[],
    obsm_keys=[],
    uns_keys=[],
):
    if batch_method is None:
        obs_keys=['leiden_r0_1', 'leiden_r0_3', 'leiden_r0_5', 'leiden_r0_7', 'leiden_r0_9']
        obsm_keys=['X_pca']
    elif batch_method == 'harmony':
        obs_keys=['leiden_hm_r0_1', 'leiden_hm_r0_3', 'leiden_hm_r0_5', 'leiden_hm_r0_7', 'leiden_hm_r0_9']
        obsm_keys=['X_pca', 'X_pca_hm']
        uns_keys=['neighbors']
    elif batch_method == 'bbknn':
        obs_keys=['leiden_bk_r0_1', 'leiden_bk_r0_3', 'leiden_bk_r0_5', 'leiden_bk_r0_7', 'leiden_bk_r0_9']
        obsm_keys=['X_pca', 'X_pca_bk']
        uns_keys=['neighbors']

    ad1 = ad.copy()
    for k in obs_keys:
        if k in ad1.obs.keys():
            del ad1.obs[k]
    for k in obsm_keys:
        if k in ad1.obsm.keys():
            del ad1.obsm[k]
    for k in uns_keys:
        if k in ad1.uns.keys():
            del ad1.uns[k]
    clear_colors(ad1)
    if ad1.raw:
        ad1.X = ad1.raw.X
        ad1.raw = None
    if out_prefix:
        ad1.write(f'{out_prefix}.processed.h5ad', compression='lzf')
    return ad1


def integrate(ads, ad_types=None, ad_prefices=None, annotations=None, batches=None, join='inner', n_hvg=4000, pool_only=False):
    n_ad = len(ads)
    if ad_types is not None:
        if isinstance(ad_types, str):
            ad_types = [ad_types] * n_ad
        elif not isinstance(ad_types, (tuple, list)) or len(ad_types) != n_ad:
            raise ValueError('invalid `ad_types` provided')
    else:
        ad_types = ['auto'] * n_ad
    if annotations is not None:
        if isinstance(annotations, str):
            annotations = [annotations] * n_ad
        elif not isinstance(annotations, (tuple, list)) or len(annotations) != n_ad:
            raise ValueError('invalid `annotations` provided')
    if batches is not None:
        if isinstance(batches, str):
            batches = [batches] * n_ad
        elif not isinstance(batches, (tuple, list)) or len(batches) != n_ad:
            raise ValueError('invalid `batches` provided')

    norm_ads = []
    for i, ad in enumerate(ads):
        ad_type = ad_types[i]
        if ad_type not in ('raw_norm', 'counts', 'norm'):
            if ad.raw and sp.issparse(ad.raw.X):
                ad_type = 'raw_norm'
            elif sp.issparse(ad.X) and np.abs(ad.X.data - ad.X.data.astype(int)).sum() < (1e-6 * ad.X.data.size):
                ad_type = 'counts'
            elif sp.issparse(ad.X) and np.abs(np.expm1(ad.X[0:10, :]).sum(axis=1).A1 - np.expm1(ad.X[0:10, :]).sum(axis=1).A1.mean()).sum() < 1e-2:
                ad_type = 'norm'
            else:
                raise ValueError(f'Cannot determine the type of anndata at position {i}')
            print(ad_type)

        if ad_type == 'raw_norm':
            norm_ad = anndata.AnnData(X=ad.raw.X, obs=ad.obs.copy(), var=ad.raw.var.copy())
        elif ad_type == 'norm':
            norm_ad = anndata.AnnData(X=ad.X, obs=ad.obs.copy(), var=ad.var.copy())
        else:
            norm_ad = anndata.AnnData(X=ad.X, obs=ad.obs.copy(), var=ad.var.copy())
            sc.pp.normalize_total(norm_ad, target_sum=1e4)
            sc.pp.log1p(norm_ad)

        post_norm_count = np.expm1(norm_ad.X[0:10, :]).sum(axis=1).A1.mean().astype(int)
        if post_norm_count != 10000:
            norm_ad.X = norm_ad.X / (post_norm_count / 1e4)

        prefix = ad_prefices[i] if ad_prefices else str(i)
        if batches and batches[i] in norm_ad.obs.columns and batches[i] != 'batch':
            if 'batch' in norm_ad.obs.columns:
                del norm_ad.obs['batch']
            norm_ad.obs.rename(columns={batches[i]: 'batch'}, inplace=True)
        if annotations and annotations[i] in norm_ad.obs.columns and annotations[i] != 'annot':
            if 'annot' in norm_ad.obs.columns:
                del norm_ad.obs['annot']
            norm_ad.obs.rename(columns={annotations[i]: 'annot'}, inplace=True)
            norm_ad.obs['annot'] = f'{prefix}_' + norm_ad.obs['annot'].astype(str)
        norm_ads.append(norm_ad)
        del norm_ad

    pooled = anndata.AnnData.concatenate(*norm_ads, batch_key='dataset', join=join)

    if ad_prefices and len(ad_prefices) == n_ad:
        pooled.obs['dataset'] = pooled.obs['dataset'].astype(str)
        for i, ad_prefix in enumerate(ad_prefices):
            k = pooled.obs['dataset'] == str(i)
            pooled.obs.loc[k, 'dataset'] = ad_prefix
        pooled.obs['dataset'] = pooled.obs['dataset'].astype('category')

    calculate_qc(pooled, flag_only=True)

    if pool_only:
        return pooled

    pooled1 = simple_default_pipeline(
        pooled, post_norm_only=True, do_clustering=False,
        batch='dataset' if batches is None else ['dataset', 'batch'],
        hvg_kw={'by_batch': ('dataset', 1), 'n_hvg': n_hvg},
        pca_kw={'remove_genes': ('mito', 'ribo')}
    )
    return pooled1
