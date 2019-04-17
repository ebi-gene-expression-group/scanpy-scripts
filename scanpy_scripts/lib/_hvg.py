"""
scanpy hvg
"""

import numpy as np
import scanpy as sc

def hvg(
        adata,
        mean_limits=(0.0125, 3),
        disp_limits=(0.5, float('inf')),
        subset=False,
        by_batch=None,
        **kwargs,
):
    """
    Wrapper function for sc.highly_variable_genes(), mainly to support searching
    by batch and pooling.
    """
    if by_batch:
        batch_name = by_batch[0]
        min_n = by_batch[1]
        k_hvg = np.zeros(adata.shape[1], dtype=int)
        for batch_value in adata.obs[batch_name].unique():
            k_v = adata.obs[batch_name] == batch_value
            hvg = sc.pp.highly_variable_genes(
                adata[k_v, :],
                min_mean=mean_limits[0],
                max_mean=mean_limits[1],
                min_disp=disp_limits[0],
                max_disp=disp_limits[1],
                inplace=False,
                **kwargs,
            )['highly_variable']
            k_hvg += hvg
        if subset:
            adata = adata[:, k_hvg >= min_n]
        else:
            adata.var['highly_variable'] = k_hvg >= min_n
    else:
        sc.pp.highly_variable_genes(
            adata,
            min_mean=mean_limits[0],
            max_mean=mean_limits[1],
            min_disp=disp_limits[0],
            max_disp=disp_limits[1],
            subset=subset,
            **kwargs,
        )
    return adata
