"""
scanpy norm
"""

import scanpy as sc
import math


def normalize(adata, log_transform=True, **kwargs):
    """
    Wrapper function for sc.pp.normalize_per_cell() and sc.pp.log1p(), mainly
    for supporting different ways of saving raw data.
    """
    sc.pp.normalize_total(adata, **kwargs)
    if log_transform:
        # Natural logarithm is the default by scanpy, if base is not set
        base = math.e
        sc.pp.log1p(adata, base=base)
        # scanpy is not setting base in uns['log1p'] keys, but later on asking for it
        if "log1p" in adata.uns_keys() and "base" not in adata.uns["log1p"]:
            # Note that setting base to None doesn't solve the problem at other modules that check for base later on
            # as adata.uns["log1p"]["base"] = None gets dropped at either anndata write or read.
            adata.uns["log1p"]["base"] = base

    return adata
