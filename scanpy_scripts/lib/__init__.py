"""
Provides exported functions
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LinearSegmentedColormap
import scanpy as sc
from scanpy.plotting._tools.scatterplots import plot_scatter

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
    cross_table,
    set_figsize,
    plot_df_heatmap,
    plot_qc,
    plot_metric_by_rank,
    plot_embedding,
    plot_diffexp,
)
from ._utils import (
    run_harmony,
    run_bbknn,
    run_seurat_integration,
    split_by_group,
    regroup,
    subsample,
    pseudo_bulk,
    LR_annotate,
    annotate,
)
from ..cmd_utils import _read_obj as read_obj
from ..cmd_utils import _write_obj as write_obj
from ..cmd_utils import switch_layer


def simple_default_pipeline(
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
        if 'mito' not in adata.var.columns:
            adata.var['mito'] = adata.var_names.str.startswith('MT-')
        if 'ribo' not in adata.var.columns:
            adata.var['ribo'] = adata.var_names.str.startswith('RPL') | adata.var_names.str.startswith('RPS')
        if 'hb' not in adata.var.columns:
            adata.var['hb'] = adata.var_names.str.startswith('HB')
        qc_tbls = sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mito', 'ribo', 'hb'], percent_top=None)
        adata.obs['n_counts'] = qc_tbls[0]['total_counts'].values
        adata.obs['n_genes'] = qc_tbls[0]['n_genes_by_counts'].values
        adata.obs['percent_mito'] = qc_tbls[0]['pct_counts_mito'].values
        adata.obs['percent_ribo'] = qc_tbls[0]['pct_counts_ribo'].values
        adata.obs['percent_hb'] = qc_tbls[0]['pct_counts_hb'].values
        adata.var['n_cells'] = qc_tbls[1]['n_cells_by_counts'].values
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
