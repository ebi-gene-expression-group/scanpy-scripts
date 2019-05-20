"""
scanpy dpt
"""

import numpy as np
import scanpy as sc


def dpt(
        adata,
        root=None,
        use_graph='neighbors',
        use_diffmap='X_diffmap',
        key_added=None,
        **kwargs,
):
    """
    Wrapper function for sc.tl.dpt, for support named slot
    """
    if root is None or not (isinstance(root, (list, tuple)) and len(root) == 2):
        root = (None, None)
    if 'iroot' not in adata.uns.keys():
        if root[0] is None:
            raise ValueError('Annotate your data with root cell first, i.e. '
                             'boolean vector `.uns["iroot"]` is required.')
        adata.uns['iroot'] = np.flatnonzero(adata.obs[root[0]] == root[1])[0]
    if use_graph != 'neighbors':
        if use_graph not in adata.uns.keys():
            raise KeyError(f'{use_graph} not found in `.uns`')
        if 'neighbors' in adata.uns.keys():
            adata.uns['neighbors_bkup'] = adata.uns['neighbors']
        adata.uns['neighbors'] = adata.uns[use_graph]
    if use_diffmap != 'X_diffmap':
        if use_diffmap not in adata.obsm.keys():
            raise KeyError(f'{use_diffmap} not found in `.obsm`')
        if 'X_diffmap' in adata.uns.keys():
            adata.obsm['X_diffmap_bkup'] = adata.obsm['X_diffmap']
        adata.obsm['X_diffmap'] = adata.obsm[use_diffmap]
    sc.tl.dpt(adata, **kwargs)
    if key_added:
        dpt_key = f'dpt_pseudotime_{key_added}'
        adata.obs[dpt_key] = adata.obs['dpt_pseudotime']
        del adata.obs['dpt_pseudotime']
    if use_graph != 'neighbors':
        del adata.uns['neighbors']
        if 'neighbors_bkup' in adata.uns.keys():
            adata.uns['neighbors'] = adata.uns['neighbors_bkup']
            del adata.uns['neighbors_bkup']
    if use_diffmap != 'X_diffmap':
        del adata.obsm['X_diffmap']
        if 'X_diffmap_bkup' in adata.obsm.keys():
            adata.obsm['X_diffmap'] = adata.obsm['X_diffma_bkup']
            del adata.obsm['X_diffmap_bkup']
    return adata
