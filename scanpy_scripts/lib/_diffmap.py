"""
scanpy diffmap
"""

import scanpy as sc
from ..cmd_utils import (
    _set_default_key,
    _restore_default_key,
    _rename_default_key,
    write_embedding,
)


def diffmap(
        adata,
        use_graph='neighbors',
        key_added=None,
        export_embedding=None,
        **kwargs,
):
    """
    Wrapper function for sc.tl.diffmap, for supporting named slot
    """
    _set_default_key(adata, 'uns', 'neighbors', use_graph)
    sc.tl.diffmap(adata, **kwargs)
    _restore_default_key(adata, 'uns', 'neighbors', use_graph)

    diffmap_key = 'X_diffmap'
    if key_added:
        diffmap_key = f'{diffmap_key}_{key_added}'
        _rename_default_key(adata, 'obsm', 'X_diffmap', diffmap_key)

    if export_embedding is not None:
        write_embedding(adata, diffmap_key, export_embedding, key_added=key_added)
    return adata
