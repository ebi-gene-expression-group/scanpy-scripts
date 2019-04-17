"""
scanpy neighbors
"""

import scanpy as sc


def neighbors(adata, n_neighbors=15, key_added=None, **kwargs):
    """
    Wrapper function for sc.pp.neighbors(), for supporting multiple n_neighbors
    """
    if not isinstance(n_neighbors, (list, tuple)):
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, **kwargs)
        if key_added:
            adata.uns[key_added] = adata.uns['neighbors']
            del adata.uns['neighbors']
    else:
        for i, n_nb in enumerate(n_neighbors):
            if key_added is None:
                key = 'neighbors_{}'.format(n_nb)
            elif not isinstance(key_added, (list, tuple)):
                key = '{}_{}'.format(key_added, n_nb)
            elif len(key_added) == len(n_neighbors):
                key = key_added[i]
            else:
                raise ValueError('`key_added` can only be None, a scalar, or an '
                                 'iterable of the same length as `n_neighbors`.')
            neighbors(
                adata,
                neighbors=n_nb,
                key_added=key,
                **kwargs,
            )
    return adata
