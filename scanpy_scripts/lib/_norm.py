"""
scanpy norm
"""

import scanpy as sc


def normalize(adata, save_raw='yes', **kwargs):
    """
    Wrapper function for sc.pp.normalize_per_cell() and sc.pp.log1p(), mainly
    for supporting different ways of saving raw data.
    """
    if save_raw == 'counts':
        adata.raw = adata
    sc.pp.normalize_total(adata, **kwargs)
    sc.pp.log1p(adata)
    if save_raw == 'yes':
        adata.raw = adata

    return adata
