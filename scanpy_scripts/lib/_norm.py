"""
scanpy norm
"""

import scanpy as sc


def normalize(adata, save_raw='yes', no_log_transform=False, **kwargs):
    """
    Wrapper function for sc.pp.normalize_per_cell() and sc.pp.log1p(), mainly
    for supporting different ways of saving raw data.
    """
    if save_raw == 'counts':
        adata.raw = adata
    sc.pp.normalize_total(adata, **kwargs)
    if not no_log_transform:
        sc.pp.log1p(adata)
    if save_raw == 'yes':
        adata.raw = adata

    return adata
