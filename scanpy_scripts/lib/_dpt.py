"""
scanpy dpt
"""

import numpy as np
import scanpy as sc
from ..obj_utils import (
    _rename_default_key,
)


def dpt(
        adata,
        root=None,
        use_diffmap='X_diffmap',
        key_added=None,
        **kwargs,
):
    """
    Wrapper function for sc.tl.dpt
    """
    if root is None or not (isinstance(root, (list, tuple)) and len(root) == 2):
        root = (None, None)
    if 'iroot' not in adata.uns.keys() and root[0] is None:
        raise ValueError('Annotate your data with root cell first, i.e. '
                         'boolean vector `.uns["iroot"]` is required.')
    if root[0] is not None:
        adata.uns['iroot'] = np.random.choice(
            np.flatnonzero(adata.obs[root[0]] == root[1]))

    sc.tl.dpt(adata, **kwargs)
    if key_added:
        dpt_key = f'dpt_pseudotime_{key_added}'
        _rename_default_key(adata.obs, 'dpt_pseudotime', dpt_key)
    
    return adata
