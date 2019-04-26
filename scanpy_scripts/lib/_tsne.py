"""
scanpy tsne
"""

import scanpy as sc


def tsne(adata, key_added=None, random_state=0, **kwargs):
    """
    Wrapper function for sc.tl.tsne, for supporting named slot of tsne embeddings
    """
    sc.tl.tsne(adata, random_state=random_state, **kwargs)
    if key_added:
        tsne_key = 'X_tsne_{}'.format(key_added)
    else:
        tsne_key = 'X_tsne_r{}'.format(random_state)
    adata.obsm[tsne_key] = adata.obsm['X_tsne']
    del adata.obsm['X_tsne']
    return adata
