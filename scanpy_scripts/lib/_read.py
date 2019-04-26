"""
Provides read_10x()
"""

import pandas as pd
import scanpy as sc


def read_10x(
        input_10x_h5,
        input_10x_mtx,
        genome='hg19',
        var_names='gene_symbols',
        extra_obs=None,
        extra_var=None
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
    return adata
