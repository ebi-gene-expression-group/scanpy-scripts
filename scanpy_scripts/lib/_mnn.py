"""
scanpy external mnn
"""

import scanpy.external as sce
import numpy as np 
import click

# Wrapper for mnn allowing use of non-standard slot

def mnn_correct(adata, key=None, key_added=None, var_subset=None, layer=None, **kwargs):
    """
    Wrapper function for sce.pp.mnn_correct(), for supporting non-standard neighbors slot
    """

    # mnn will use .X, so we need to put other layers there for processing

    if layer:
        adata.layers['X_backup'] = adata.X
        adata.X = adata.layers[layer]

    # mnn_correct() wants batches in separate adatas

    batches = np.unique(adata.obs[key])
    alldata = []
    for batch in batches:
        alldata.append( adata[adata.obs[key] == batch,] )

    # Process var_subset into a list of strings that can be provided to
    # mnn_correct()

    if var_subset is not None and len(var_subset) > 0 and var_subset[0] is not None:
   
        subset = []
 
        for name, values in var_subset :
            if name in adata.var:
                if adata.var[name].dtype == 'bool':
                    values = [ True if x.lower() == "true" else x for x in values  ]
            else:
                raise click.ClickException(f'Var "{name}" unavailable')

            ind = [ x in values for x in adata.var[name]  ]
            subset = subset + adata.var.index[ ind  ].to_list()

        var_subset = set(subset)
        print('Will use %d selected genes for MNN' % len(var_subset))

    else:
        var_subset = None

    # Here's the main bit

    cdata = sce.pp.mnn_correct(*alldata, var_subset = var_subset, do_concatenate = True, index_unique = None, **kwargs)    
    
    # If user has specified key_added = X then they want us to overwrite .X,
    # othwerwise copy the .X to a named layer of the original object. In either
    # case make sure obs and var are the same as the original.

    if key_added is None or key_added != 'X':

        mnn_key = 'mnn'
        if layer:
            mnn_key = f"{mnn_key}_{layer}"
        
            # Layers is set (so we're not storing computed results in the .X,
            # and we had to overwrite .X to run mnn), and key_added shows we're
            # not storing in the .X, so we need to restore from the backup.

            adata.X = adata.layers['X_backup']

        if key_added:
            mnn_key = f"{mnn_key}_{key_added}"
   
        adata.layers[mnn_key] = cdata[0][adata.obs.index, adata.var.index].X

    else:
        adata.X = cdata[0][adata.obs.index, adata.var.index].X

    # Delete the backup of .X if we needed one 

    if layer:
        del adata.layers['X_backup'] 

    return adata
