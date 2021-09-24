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
        adata.obs = _fix_booleans(adata.obs.merge(
            obs_tbl,
            how='left',
            left_index=True,
            right_index=True,
            suffixes=(False, False),
        ))

    if extra_var:
        var_tbl = pd.read_csv(extra_var, sep='\t', header=0, index_col=0)
        adata.var = _fix_booleans(adata.var.merge(
            var_tbl,
            how='left',
            left_index=True,
            right_index=True,
            suffixes=(False, False),
        ))
    return adata

def _fix_booleans(df):
  for var in df.columns:
    if (df[var].dtype.kind == 'O' and 
        df[var].dtype.name == 'object' and 
        set(pd.Categorical(df[var][df[var] != 'nan'])).issubset(set(['True', 'False']))
        ):
      d = {'False': True, 'False': False}
      df[var] = df[var].map(d)
  return(df)
