"""
scanpy combat
"""

import scanpy as sc

# Wrapper for mnn allowing use of non-standard slot

def combat(adata, key=None, key_added=None, layer=None, **kwargs):
    """
    Wrapper function for scanpy.pp.combat(), for supporting non-standard slots
    """

    # If layer is set then we have to move the contents of that layer into
    # .X for analysis. We back up the original .X, but only if the user hasn't
    # specified to overwrite it anyway.

    if layer:
        if key_added and key_added != 'X':
            adata.layers['X_backup'] = adata.X

        adata.X = adata.layers[layer]
    
    # If we're storing results in .X (whether from .X or from a layer), run in
    # place to save copying objects.

    if key_added and key_added == 'X':
        sc.pp.combat(adata, key = key, **kwargs)    
    
    # If we're storing in 'layers' (key_added is not set, or is not X, then
    # don't run in place, and put the matrix in the specified layer. 

    else:

        cdata = sc.pp.combat(adata, key=key, inplace = False, **kwargs)    
            
        combat_key = 'combat'
        if layer:
            combat_key = f"{combat_key}_{layer}"
        
            # If we ran from a layer, restore the .X we had to overwrite
            
            adata.X = adata.layers['X_backup']    
            del adata.layers['X_backup'] 
            
        if key_added:
            combat_key = f"{combat_key}_{key_added}"

        adata.layers[combat_key] = cdata

    return adata
