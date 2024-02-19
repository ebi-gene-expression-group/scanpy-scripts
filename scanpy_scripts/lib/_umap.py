"""
scanpy umap
"""

import scanpy as sc
from ..obj_utils import (
    _set_default_key,
    _restore_default_key,
    _backup_obsm_key,
    _rename_obsm_key,
    _delete_obsm_backup_key,
    write_embedding,
)


def umap(
    adata,
    key_added=None,
    random_state=0,
    export_embedding=None,
    **kwargs,
):
    """
    Wrapper function for sc.tl.umap, for supporting named slot of umap embeddings
    """
    if not isinstance(random_state, (list, tuple)):
        _backup_obsm_key(adata, "X_umap")

        sc.tl.umap(adata, random_state=random_state, **kwargs)

        umap_key = "X_umap"
        if key_added:
            umap_key = f"X_umap_{key_added}"
            _rename_obsm_key(adata, "X_umap", umap_key)
        else:
            _delete_obsm_backup_key(adata, "X_umap")

        if export_embedding is not None:
            write_embedding(adata, umap_key, export_embedding, key_added=key_added)
    else:
        for i, rseed in enumerate(random_state):
            if key_added is None:
                umap_key = f"r{rseed}"
            elif not isinstance(key_added, (list, tuple)):
                umap_key = f"{key_added}_r{rseed}"
            elif len(key_added) == len(random_state):
                umap_key = key_added[i]
            else:
                raise ValueError(
                    "`key_added` can only be None, a scalar, or an "
                    "iterable of the same length as `random_state`."
                )
            umap(
                adata,
                key_added=umap_key,
                random_state=rseed,
                **kwargs,
            )
    return adata
