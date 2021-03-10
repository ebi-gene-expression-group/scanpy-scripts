"""
scanpy external scrublet
"""

import scanpy.external as sce
import numpy as np
from ..obj_utils import write_obs

# Wrapper for bbknn allowing use of non-standard slot

def scrublet(adata, filter=False, export_table=None, **kwargs):
    """
    Wrapper function for sce.pp.scrublet(), to allow filtering of resulting object 
    """

    sce.pp.scrublet(adata, **kwargs)    

    # Do any export before optional filtering

    if export_table:
        write_obs(adata, ['doublet_score', 'predicted_doublet'], export_table)
    
    # Filter out predited doublets

    if filter:
        adata._inplace_subset_obs(np.invert(adata.obs['predicted_doublet']))
    
    return adata
