"""
scanpy pca
"""

import scanpy as sc
from ..cmd_utils import write_embedding

def pca(adata, key_added=None, export_embedding=None, **kwargs):
    """
    Wrapper function for sc.pp.pca, for supporting named slot
    """
    if key_added:
        if 'X_pca' in adata.obsm.keys():
            adata.obsm['X_pca_bkup'] = adata.obsm['X_pca']
        sc.pp.pca(adata, **kwargs)
        pca_key = f'X_pca_{key_added}'
        adata.obsm[pca_key] = adata.obsm['X_pca']
        del adata.obsm['X_pca']
        if 'X_pca_bkup' in adata.obsm.keys():
            adata.obsm['X_pca'] = adata.obsm['X_pca_bkup']
            del adata.obsm['X_pca_bkup']
    else:
        sc.pp.pca(adata, **kwargs)
        pca_key = 'X_pca'

    if export_embedding is not None:
        write_embedding(adata, pca_key, export_embedding, key_added=key_added)
    return adata
