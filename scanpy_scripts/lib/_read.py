"""
Provides read_10x()
"""

import pandas as pd
import scanpy as sc
import logging


def read_10x(
        input_10x_h5,
        input_10x_mtx,
        genome='hg19',
        var_names='gene_symbols',
        extra_obs=None,
        extra_var=None,
        gene_startswith='mito:MT-'
):
    """
    Wrapper function for sc.read_10x_h5() and sc.read_10x_mtx(), mainly to
    support adding extra metadata
    """
    if input_10x_h5 is not None:
        adata = sc.read_10x_h5(input_10x_h5, genome=genome)
    elif input_10x_mtx is not None:
        adata = sc.read_10x_mtx(input_10x_mtx, var_names=var_names)

    if extra_obs:
        obs_tbl = pd.read_csv(extra_obs, sep='\t', header=0, index_col=0)
        adata.obs = adata.obs.merge(
            obs_tbl,
            how='left',
            left_index=True,
            right_index=True,
            suffixes=(False, False),
        )

    if extra_var:
        var_tbl = pd.read_csv(extra_var, sep='\t', header=0, index_col=0)
        adata.var = adata.var.merge(
            var_tbl,
            how='left',
            left_index=True,
            right_index=True,
            suffixes=(False, False),
        )

    gene_name = 'index'
    if var_names == 'gene_ids'
        gene_name = 'gene_symbols'
    
    gene_startswith_list = gene_startswith.split(",") # mito:MT-,ERCC:ERCC
    for gene_sw in gene_startswith_list:
        if len(gene_sw.split(":"))==2:
            var_cat = gene_sw.split(":")[0]
            starts_with = gene_sw.split(":")[1]
            try:
                gene_names = getattr(adata.var, gene_name)
                k_cat = gene_names.str.startswith(starts_with)
                if k_cat.sum() > 0:
                adata.var[var_cat] = k_cat
            else:
                logging.warning('No {} genes found, skip calculating expression of {} genes'.format(starts_with, var_cat))
            except AttributeError:
                logging.warning(
                    'Specified gene column [%s] not found, skip calculating '
                    'expression of defined genes', gene_name)

    if 'n_genes' not in adata.obs.columns:
        sc.pp.filter_cells(adata, min_genes=0)
    if 'n_counts' not in adata.obs.columns:
        sc.pp.filter_cells(adata, min_counts=0)
    if 'n_cells' not in adata.var.columns:
        sc.pp.filter_genes(adata, min_cells=0)
    if 'n_counts' not in adata.var.columns:
        sc.pp.filter_genes(adata, min_counts=0)

    return adata
