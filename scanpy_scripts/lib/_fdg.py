"""
scanpy fdg
"""

import scanpy as sc
from ..cmd_utils import write_embedding


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
        adata.obsm[fdg_key] = adata.obsm[f'X_draw_graph_{layout}']
        del adata.obsm[f'X_draw_graph_{layout}']
    if export_embedding is not None:
        if key_added:
            if export_embedding.endswith('.tsv'):
                export_embedding = export_embedding[0:-4]
            export_embedding = f'{export_embedding}_{key_added}.tsv'
        write_embedding(adata, fdg_key, export_embedding)
    return adata
