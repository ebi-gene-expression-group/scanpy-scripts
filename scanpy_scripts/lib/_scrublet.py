"""
scanpy external scrublet
"""

import scanpy as sc
import scanpy.external as sce
import numpy as np
from ..obj_utils import write_obs
import anndata
import pandas as pd

# Wrapper for scrublet allowing text export and filtering

def scrublet(adata, adata_sim=None, filter=False, batch_key=None, export_table=None, **kwargs):
    """
    Wrapper function for sce.pp.scrublet(), to allow filtering of resulting object
    """

    # Do we need to read an object with the doublet simulations?

    if adata_sim:
        adata_sim = sc.read(adata_sim)

    # Scrublet shouldn't be run on multi-batch data, so we run the batches
    # separately and copy the stats back to the input object

    alldata = []
    if batch_key is not None:
        if batch_key in adata.obs.keys():
            print("batch key %s is in obs" % batch_key)
        else:
            print("batch key %s is NOT in obs" % batch_key)
        
        batches = np.unique(adata.obs[batch_key])
    
        # Run Scrublet independently on batches and return just the
        # scrublet-relevant parts of the objects to add to the input object

        def get_adata_scrub_parts(ad):
            return {'obs': ad.obs, 'uns': ad.uns['scrublet']}

        scrubbed = [ get_adata_scrub_parts(sce.pp.scrublet(adata[adata.obs[batch_key] == batch,], adata_sim=adata_sim, copy = True, **kwargs)) for batch in batches ] 
        scrubbed_obs = pd.concat([ scrub['obs'] for scrub in scrubbed])        

        # Now reset the obs to get the scrublet scores
        
        adata.obs = scrubbed_obs.loc[adata.obs_names.values]

        # Save the .uns from each batch separately
    
        adata.uns['scrublet'] = dict(zip(batches, [ scrub['uns'] for scrub in scrubbed ]))

    else:
        sce.pp.scrublet(adata, adata_sim=adata_sim, **kwargs)
    
    # Do any export before optional filtering

    if export_table:
        write_obs(adata, ["doublet_score", "predicted_doublet"], export_table)

    # Filter out predited doublets

    if filter:
        adata._inplace_subset_obs(np.invert(adata.obs["predicted_doublet"]))

    return adata

# Run the doublet simulation.

def scrublet_simulate_doublets(adata, **kwargs):
    adata_sim = sce.pp.scrublet_simulate_doublets(adata, **kwargs)
    adata._init_as_actual(
        X=adata_sim.X, obs=adata_sim.obs, obsm=adata_sim.obsm, uns=adata.uns
    )

# Just absorb the extra plotting args before passing to
# scanpy.external.pl.scrublet_score_distribution


def plot_scrublet(
    adata, scale_hist_obs="log", scale_hist_sim="linear", fig_size=(8, 3), **kwargs
):
    """
    Wrapper function for sce.pl.scrublet_score_distribution(), to allow
    plotting of score distribution
    """
    sce.pl.scrublet_score_distribution(
        adata,
        scale_hist_obs=scale_hist_obs,
        scale_hist_sim=scale_hist_sim,
        figsize=fig_size,
        **kwargs
    )
