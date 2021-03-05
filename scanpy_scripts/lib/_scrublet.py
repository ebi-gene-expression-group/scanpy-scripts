"""
scanpy external scrublet
"""

import scanpy.external as sce
import numpy as np

# Wrapper for bbknn allowing use of non-standard slot

def scrublet(adata, filter=False, **kwargs):
    """
    Wrapper function for sce.pp.scrublet(), to allow filtering of resulting object 
    """

    sce.pp.scrublet(adata, **kwargs)    

    # Filter out predited doublets

    if filter:
        adata._inplace_subset_obs(np.invert(adata.obs['predicted_doublet']))
    
    return adata
