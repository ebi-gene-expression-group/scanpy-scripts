"""
scanpy paga
"""

import scanpy as sc

def paga(
        adata,
        use_graph='neighbors',
        key_added=None,
        **kwargs,
):
    """
    Wrapper function for sc.tl.paga, for supporting named slot
    """
    if use_graph != 'neighbors':
        if use_graph not in adata.uns.keys():
            raise KeyError(f'{use_graph} not found in `.uns`')
        if 'neighbors' in adata.uns.keys():
            adata.uns['neighbors_bkup'] = adata.uns['neighbors']
        adata.uns['neighbors'] = adata.uns[use_graph]
    sc.tl.paga(adata, **kwargs)
    if key_added:
        uns_key = f'paga_{key_added}'
        adata.uns[uns_key] = adata.uns['paga']
        del adata.uns['paga']
    if use_graph != 'neighbors':
        del adata.uns['neighbors']
        if 'neighbors_bkup' in adata.uns.keys():
            adata.uns['neighbors'] = adata.uns['neighbors_bkup']
            del adata.uns['neighbors_bkup']
    return adata
