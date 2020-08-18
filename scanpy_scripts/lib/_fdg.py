"""
scanpy fdg
"""

import scanpy as sc
from ..obj_utils import (
    _backup_obsm_key,
    _delete_obsm_backup_key,
    _rename_obsm_key,
    write_embedding,
)


def fdg(
        adata,
        layout='fa',
        key_added_ext=None,
        random_state=0,
        export_embedding=None,
        **kwargs
):
    """
    Wrapper function for sc.tl.draw_graph, for supporting named slot of fdg
    embeddings.
    """
    sc.tl.draw_graph(
        adata,
        layout=layout,
        key_added_ext=key_added_ext,
        random_state=random_state,
        **kwargs,
    )

    fdg_key  = f'X_draw_graph_{key_added_ext or layout}'

    if export_embedding is not None:
        write_embedding(adata, fdg_key, export_embedding, key_added=key_added_ext)
    return adata
