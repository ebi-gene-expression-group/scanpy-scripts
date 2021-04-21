"""
scanpy norm
"""

import scanpy as sc


def normalize(adata, save_raw='yes', save_raw_location='raw', log_transform=True, **kwargs):
    """
    Wrapper function for sc.pp.normalize_per_cell() and sc.pp.log1p(), mainly
    for supporting different ways of saving raw data.
    """
    if save_raw == 'counts':
        if save_raw_location == 'raw':
            adata.raw = adata
        else:
            adata.layers[save_raw_location] = adata.X
    sc.pp.normalize_total(adata, **kwargs)
    if log_transform:
        sc.pp.log1p(adata)
    if save_raw == 'yes':
        if save_raw_location == 'raw':
            adata.raw = adata
        else:
            adata.layers[save_raw_location] = adata.X

    return adata
