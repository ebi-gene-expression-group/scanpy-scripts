"""
scanpy umap
"""

import scanpy as sc


def umap(
        adata,
        use_graph='neighbors',
        key_added=None,
        init_pos='spectral',
        random_state=0,
        **kwargs
):
    """
    Wrapper function for sc.tl.umap, for supporting named slot of umap embeddings
    """
    if use_graph != 'neighbors':
        adata.uns['neighbors'] = adata.uns[use_graph]
    sc.tl.umap(adata, init_pos=init_pos, random_state=random_state, **kwargs)
    if key_added:
        umap_key = 'X_umap_{}'.format(key_added)
    else:
        umap_key = 'X_umap_{}_r{}'.format(init_pos, random_state)
    adata.obsm[umap_key] = adata.obsm['X_umap']
    del adata.obsm['X_umap']
    return adata
