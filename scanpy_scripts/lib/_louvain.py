"""
scanpy louvain
"""

import scanpy as sc
from ..obj_utils import write_obs


def louvain(
        adata,
        resolution,
        neighbors_key=None,
        obsp=None,
        key_added=None,
        export_cluster=None,
        **kwargs
):
    """
    Wrapper function for sc.tl.louvain, for supporting multiple resolutions.
    """
    keys = []
    if kwargs['restrict_to'] and not kwargs['restrict_to'][0]:
        kwargs['restrict_to'] = None
    
    if not isinstance(resolution, (list, tuple)):
        if key_added is not None and not key_added.startswith('louvain_'):
            key_added = f'louvain_{key_added}'
        elif key_added is None:
            key_added = 'louvain'
        sc.tl.louvain(
            adata,
            resolution=resolution,
            key_added=key_added,
            neighbors_key=neighbors_key,
            obsp=obsp,
            **kwargs
        )
        keys.append(key_added)
    else:
        for i, res in enumerate(resolution):

            res_key = str(res).replace('.', '_')

            if key_added is None:
                graph_key = ('_' + f'{neighbors_key or obsp}') if neighbors or obsp else ''
                key = f'louvain{graph_key}_r{res_key}'    
            elif not isinstance(key_added, (list, tuple)):
                key = f'louvain_{key_added}_r{res_key}'
            elif len(key_added) == len(resolution):
                key = key_added[i]
            else:
                raise ValueError('`key_added` can only be None, a scalar, or an '
                                 'iterable of the same length as `resolution`.')
            keys.extend(louvain(
                adata,
                resolution=res,
                neighbors_key=neighbors_key,
                obsp=obsp,
                key_added=key,
                **kwargs,
            ))

    if export_cluster:
        write_obs(adata, keys, export_cluster)

    return keys
