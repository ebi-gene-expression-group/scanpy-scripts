"""
scanpy norm
"""

import scanpy as sc


def normalize(adata, log_transform=True, **kwargs):
    """
    Wrapper function for sc.pp.normalize_per_cell() and sc.pp.log1p(), mainly
    for supporting different ways of saving raw data.
    """
    sc.pp.normalize_total(adata, **kwargs)
    if log_transform:
        sc.pp.log1p(adata)

    return adata
