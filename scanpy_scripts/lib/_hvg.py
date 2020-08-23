"""
scanpy hvg
"""

import numpy as np
import scanpy as sc

def hvg(
        adata,
        mean_limits=(0.0125, 3),
        disp_limits=(0.5, float('inf')),
        **kwargs,
):
    """
    Wrapper function for sc.highly_variable_genes()
    """

    # Check for n_top_genes beeing greater than the total genes

    if 'n_top_genes' in kwargs and kwargs['n_top_genes'] is not None:
        kwargs['n_top_genes'] = min(adata.n_vars, kwargs['n_top_genes'])
    
    sc.pp.highly_variable_genes(
        adata,
        min_mean=mean_limits[0],
        max_mean=mean_limits[1],
        min_disp=disp_limits[0],
        max_disp=disp_limits[1],
        **kwargs,
    )
    
    return adata
