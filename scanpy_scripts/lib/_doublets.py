"""
Remove doublets
"""

import numpy as np
from scipy.stats import norm
import scanpy as sc
import scrublet as scr
from statsmodels.stats.multitest import multipletests

def test_outlier(x, upper_mad_only=True):
    med = np.median(x)
    if upper_mad_only:
        mad = np.median(x[x>med] - med) * 1.4826
    else:
        mad = np.median(np.abs(x - med)) * 1.4826
    pvals = 1 - norm.cdf(x, loc=med, scale=mad)
    bh_pvals = multipletests(pvals, method='fdr_bh')[1]
    return pvals, bh_pvals


def run_scrublet(adata, resolution_function=None):
    old_verbosity = sc.settings.verbosity
    sc.settings.verbosity = 1
    if resolution_function is None:
        resolution_function = lambda x: np.maximum(np.maximum(np.log10(x)-1, 0)**2, 0.1) 
    scrub = scr.Scrublet(adata.X)
    ds, pd = scrub.scrub_doublets(verbose=False)
    adata.obs['scrublet_score'] = ds

    adata_copy = adata.copy()
    sc.pp.filter_genes(adata_copy, min_cells=3)
    sc.pp.normalize_total(adata_copy, target_sum=1e4, fraction=0.9)
    sc.pp.log1p(adata_copy)
    sc.pp.highly_variable_genes(adata_copy, min_mean=0.0125, max_mean=3, min_disp=0.5, subset=True)
    sc.pp.scale(adata_copy, zero_center=False)
    sc.pp.pca(adata_copy, svd_solver='arpack', zero_center=False)
    sc.pp.neighbors(adata_copy, n_pcs=30)
    sc.tl.umap(adata_copy)
    sc.tl.leiden(adata_copy, resolution=1)
    for clst in np.unique(adata_copy.obs['leiden']):
        clst_size = sum(adata_copy.obs['leiden'] == clst)
        sc.tl.leiden(adata_copy, restrict_to=('leiden', [clst]), resolution=resolution_function(clst_size), key_added='leiden_R')
        adata_copy.obs['leiden'] = adata_copy.obs['leiden_R']
    clst_meds = []
    for clst in np.unique(adata_copy.obs['leiden']):
        k = adata_copy.obs['leiden'] == clst
        clst_med = np.median(adata_copy.obs.loc[k, 'scrublet_score'])
        adata_copy.obs.loc[k, 'cluster_scrublet_score'] = clst_med
        clst_meds.append(clst_med)
    clst_meds = np.array(clst_meds)
    pvals, bh_pvals = test_outlier(clst_meds)
    for i, clst in enumerate(np.unique(adata_copy.obs['leiden'])):
        k = adata_copy.obs['leiden'] == clst
        adata_copy.obs.loc[k, 'pval'] = pvals[i]
        adata_copy.obs.loc[k, 'bh_pval'] = bh_pvals[i]
    sc.settings.verbosity = old_verbosity
    adata.obs['cluster_scrublet_score'] = adata_copy.obs['cluster_scrublet_score']
    adata.obs['doublet_bh_pval'] = adata_copy.obs['bh_pval']
    del adata_copy
