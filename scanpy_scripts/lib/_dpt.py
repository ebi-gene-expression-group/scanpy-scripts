"""
scanpy dpt
"""

import numpy as np
import scanpy as sc
from ..cmd_utils import (
    _set_default_key,
    _restore_default_key,
)


def dpt(
        adata,
        root=None,
        use_graph='neighbors',
        use_diffmap='X_diffmap',
        key_added=None,
        **kwargs,
):
    """
    Wrapper function for sc.tl.dpt, for support named slot
    """
    if root is None or not (isinstance(root, (list, tuple)) and len(root) == 2):
        root = (None, None)
    if 'iroot' not in adata.uns.keys():
        if root[0] is None:
            raise ValueError('Annotate your data with root cell first, i.e. '
                             'boolean vector `.uns["iroot"]` is required.')
        adata.uns['iroot'] = np.flatnonzero(adata.obs[root[0]] == root[1])[0]

    _set_default_key(adata, 'uns', 'neighbors', use_graph)
    _set_default_key(adata, 'obsm', 'X_diffmap', use_diffmap)
    sc.tl.dpt(adata, **kwargs)
    _restore_default_key(adata, 'uns', 'neighbors', use_graph)
    _restore_default_key(adata, 'obsm', 'X_diffmap', use_diffmap)

    if key_added:
        dpt_key = f'dpt_pseudotime_{key_added}'
        adata.obs[dpt_key] = adata.obs['dpt_pseudotime']
        del adata.obs['dpt_pseudotime']
    return adata
