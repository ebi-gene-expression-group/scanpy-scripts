"""
scanpy tsne
"""

import scanpy as sc
from ..cmd_utils import (
    _rename_default_key,
    write_embedding,
)


def tsne(
        adata,
        key_added=None,
        random_state=0,
        export_embedding=None,
        **kwargs,
):
    """
    Wrapper function for sc.tl.tsne, for supporting named slot of tsne embeddings
    """
    if not isinstance(random_state, (list, tuple)):
        sc.tl.tsne(adata, random_state=random_state, **kwargs)

        tsne_key = 'X_tsne'
        if key_added:
            tsne_key = f'X_tsne_{key_added}'
            _rename_default_key(adata, 'obsm', 'X_tsne', tsne_key)

        if export_embedding is not None:
            write_embedding(adata, tsne_key, export_embedding, key_added=key_added)
    else:
        for i, rseed in enumerate(random_state):
            if key_added is None:
                tsne_key = f'r{rseed}'
            elif not isinstance(key_added, (list, tuple)):
                tsne_key = f'{key_added}_r{rseed}'
            elif len(key_added) == len(random_state):
                tsne_key = key_added[i]
            else:
                raise ValueError('`key_added` can only be None, a scalar, or '
                                 'an iterable of the same length as '
                                 '`random_state`.')
            tsne(
                adata,
                key_added=tsne_key,
                random_state=rseed,
                **kwargs,
            )
    return adata
