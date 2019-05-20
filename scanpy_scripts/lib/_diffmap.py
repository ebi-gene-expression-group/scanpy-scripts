"""
scanpy diffmap
"""

import scanpy as sc
from ..cmd_utils import write_embedding


def diffmap(
        adata,
        use_graph='neighbors',
        key_added=None,
        export_embedding=None,
        **kwargs,
):
    """
    Wrapper function for sc.tl.diffmap, for supporting named slot
    """
    if use_graph != 'neighbors':
        if use_graph not in adata.uns.keys():
            raise KeyError(f'{use_graph} not found in `.uns`')
        if 'neighbors' in adata.uns.keys():
            adata.uns['neighbors_bkup'] = adata.uns['neighbors']
        adata.uns['neighbors'] = adata.uns[use_graph]
    sc.tl.diffmap(adata, **kwargs)
    diffmap_key = f'X_diffmap'
    if key_added:
        diffmap_key = f'X_diffmap_{key_added}'
        adata.obsm[diffmap_key] = adata.obsm['X_diffmap']
        del adata.obsm['X_diffmap']
    if export_embedding is not None:
        if key_added:
            if export_embedding.endswith('.tsv'):
                export_embedding = export_embedding[0:-4]
            export_embedding = f'{export_embedding}_{key_added}.tsv'
        write_embedding(adata, diffmap_key, export_embedding)
    return adata
