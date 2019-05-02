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
    if not isinstance(random_state, (list, tuple)):
        sc.tl.umap(adata, init_pos=init_pos, random_state=random_state, **kwargs)
        if key_added:
            umap_key = f'X_umap_{key_added}'
            adata.obsm[umap_key] = adata.obsm['X_umap']
            del adata.obsm['X_umap']
    else:
        for i, rseed in enumerate(random_state):
            if key_added is None:
                umap_key = f'r{rseed}'
            elif not isinstance(key_added, (list, tuple)):
                umap_key = f'{key_dded}_r{rseed}'
            elif len(key_added) == len(random_state):
                umap_key = key_added[i]
            else:
                raise ValueError('`key_added` can only be None, a scalar, or an '
                                 'iterable of the same length as `random_state`.')
            umap(
                adata,
                use_graph=use_graph,
                init_pos=init_pos,
                key_added=umap_key,
                random_state=rseed,
                **kwargs,
            )
    return adata
