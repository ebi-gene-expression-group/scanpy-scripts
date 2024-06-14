"""
Provides read_10x()
"""

import pandas as pd
import scanpy as sc


def read_10x(
    input_10x_h5,
    input_10x_mtx,
    genome="hg19",
    var_names="gene_symbols",
    extra_obs=None,
    extra_var=None,
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
        obs_tbl = pd.read_csv(extra_obs, sep="\t", header=0, index_col=0)
        adata.obs = adata.obs.merge(
            obs_tbl,
            how="left",
            left_index=True,
            right_index=True,
            suffixes=(False, False),
        )

    if extra_var:
        var_tbl = pd.read_csv(extra_var, sep="\t", header=0, index_col=0)
        mixed_columns = columns_with_multiple_dtypes(var_tbl)

        # Convert mixed dtype columns to 'object' type to preserve all information
        for column in mixed_columns:
            var_tbl[column] = var_tbl[column].astype('string')
    
        adata.var = adata.var.merge(
            var_tbl,
            how="left",
            left_index=True,
            right_index=True,
            suffixes=(False, False),
        )
    return adata


def columns_with_multiple_dtypes(df):
    mixed_dtype_columns = []
    for column in df.columns:
        # Get unique dtypes in the column
        unique_dtypes = df[column].apply(type).unique()
        if len(unique_dtypes) > 1:
            mixed_dtype_columns.append(column)
    return mixed_dtype_columns

