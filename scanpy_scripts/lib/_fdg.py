"""
scanpy fdg
"""

import scanpy as sc
from ..cmd_utils import (
    _rename_default_key,
    write_embedding,
)


def fdg(
        adata,
        use_graph='neighbors',
        layout='fa',
        key_added=None,
        random_state=0,
        export_embedding=None,
        **kwargs
):
    """
    Wrapper function for sc.tl.draw_graph, for supporting named slot of fdg
    embeddings.
    """
    adj_mat = None
    if use_graph != 'neighbors':
        if use_graph not in adata.uns.keys():
            raise KeyError(f'"{use_graph}" is not a valid key of `.uns`.')
        adj_mat = adata.uns[use_graph]['connectivities']
    sc.tl.draw_graph(
        adata,
        layout=layout,
        random_state=random_state,
        adjacency=adj_mat,
        **kwargs,
    )

    fdg_key = f'X_draw_graph_{layout}'
    if key_added:
        fdg_key = f'X_draw_graph_{layout}_{key_added}'
        _rename_default_key(adata, 'obsm', f'X_draw_graph_{layout}', fdg_key)

    if export_embedding is not None:
        write_embedding(adata, fdg_key, export_embedding, key_added=key_added)
    return adata
