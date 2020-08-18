"""
scanpy diffmap
"""

import scanpy as sc
from ..obj_utils import (
    _rename_obsm_key,
    write_embedding,
)


def diffmap(
        adata,
        key_added=None,
        export_embedding=None,
        **kwargs,
):
    """
    Wrapper function for sc.tl.diffmap, for supporting named slot
    """
    sc.tl.diffmap(adata, **kwargs)
    
    diffmap_key = 'X_diffmap'
    if key_added:
        diffmap_key = f'{diffmap_key}_{key_added}'
        _rename_obsm_key(adata, 'X_diffmap', diffmap_key)

    if export_embedding is not None:
        write_embedding(adata, diffmap_key, export_embedding, key_added=key_added)
    return adata
