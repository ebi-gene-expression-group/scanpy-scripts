"""
scanpy paga
"""

import scanpy as sc
from ..cmd_utils import (
    _backup_default_key,
    _delete_backup_key,
    _rename_default_key,
    _set_default_key,
    _restore_default_key,
)


def paga(
        adata,
        use_graph='neighbors',
        key_added=None,
        **kwargs,
):
    """
    Wrapper function for sc.tl.paga, for supporting named slot
    """
    _set_default_key(adata.uns, 'neighbors', use_graph)
    _backup_default_key(adata.uns, 'paga')

    sc.tl.paga(adata, **kwargs)

    _restore_default_key(adata.uns, 'neighbors', use_graph)

    if key_added:
        paga_key = f'paga_{key_added}'
        _rename_default_key(adata.uns, 'paga', paga_key)
    else:
        _delete_backup_key(adata.uns, 'paga')

    return adata


def plot_paga(
        adata,
        prefix='paga',
        use_key='paga',
        **kwargs,
):
    """Make PAGA plot
    """
    _set_default_key(adata.uns, 'paga', use_key)
    if prefix == 'paga':
        sc.pl.paga(adata, **kwargs)
    elif prefix == 'paga_compare':
        sc.pl.paga_compare(adata, **kwargs)
    else:
        raise ValueError(f'{prefix}: unknown plot type for PAGA')
    _restore_default_key(adata.uns, 'paga', use_key)
