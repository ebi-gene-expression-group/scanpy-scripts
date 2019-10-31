"""
scanpy pca
"""

import scanpy as sc
from ..obj_utils import write_embedding

def pca(adata, key_added=None, remove_cc=False, export_embedding=None, **kwargs):
    """
    Wrapper function for sc.pp.pca, for supporting named slot
    """
    remove_cc_from_hvg = (
        remove_cc and 'cc' in adata.var.columns
        and (adata.var['cc'].dtype.kind == 'b' or isinstance(adata.var['cc'][0], bool))
    )
    if remove_cc_from_hvg:
        adata.var['hvg_full'] = adata.var['highly_variable'].copy()
        adata.var['highly_variable'] = adata.var['highly_variable'] & ~adata.var['cc']

    # omit "svd_solver" to let scanpy choose automatically
    if 'svd_solver' in kwargs and kwargs['svd_solver'] == 'auto':
        del kwargs['svd_solver']

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

    if remove_cc_from_hvg:
        adata.var['highly_variable'] = adata.var['hvg_full']
        del adata.var['hvg_full']

    if export_embedding is not None:
        write_embedding(adata, pca_key, export_embedding, key_added=key_added)
    return adata
