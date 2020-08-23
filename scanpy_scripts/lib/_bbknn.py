"""
scanpy external bbknn
"""

import scanpy.external as sce

from ..obj_utils import (
    _backup_default_key,
    _delete_backup_key,
    _rename_default_key,
)

# Wrapper for bbknn allowing use of non-standard slot

def bbknn(adata, key=None, key_added=None, **kwargs):
    """
    Wrapper function for sce.pp.bbknn(), for supporting non-standard neighbors slot
    """

    _backup_default_key(adata.uns, 'neighbors')
    sce.pp.bbknn(adata, batch_key = key, **kwargs)    

    if key_added:
        _rename_default_key(adata.uns, 'neighbors',  f'neighbors_{key_added}')
    else:
        _delete_backup_key(adata.uns, 'neighbors')
    
    return adata
