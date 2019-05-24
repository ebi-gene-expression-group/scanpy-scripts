"""
scanpy leiden
"""

import scanpy as sc


def leiden(adata, resolution, use_graph=None, key_added=None, **kwargs):
    """
    Wrapper function for sc.tl.leiden, for supporting multiple resolutions.
    """
    if kwargs.get('restrict_to', None) and not kwargs['restrict_to'][0]:
        kwargs['restrict_to'] = None
    adj_mat = None
    if use_graph:
        if use_graph not in adata.uns:
            raise KeyError(f'"{use_graph}" is not a valid key of `.uns`.')
        adj_mat = adata.uns[use_graph]['connectivities']
    if not isinstance(resolution, (list, tuple)):
        if key_added is not None and not key_added.startswith('leiden_'):
            key_added = f'leiden_{key_added}'
        sc.tl.leiden(
            adata,
            resolution=resolution,
            adjacency=adj_mat,
            key_added=key_added,
            **kwargs
        )
    else:
        for i, res in enumerate(resolution):
            res_key = str(res).replace('.', '_')
            if key_added is None:
                graph_key = ('_' + use_graph) if use_graph else ''
                key = f'leiden{graph_key}_r{res_key}'
            elif not isinstance(key_added, (list, tuple)):
                key = f'leiden_{key_added}_r{res_key}'
            elif len(key_added) == len(resolution):
                key = key_added[i]
            else:
                raise ValueError('`key_added` can only be None, a scalar, or an '
                                 'iterable of the same length as `resolution`.')
            leiden(
                adata,
                resolution=res,
                use_graph=use_graph,
                key_added=key,
                **kwargs,
            )
    return adata
