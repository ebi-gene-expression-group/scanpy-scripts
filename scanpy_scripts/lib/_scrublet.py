"""
scanpy external scrublet
"""

import scanpy as sc
import scanpy.external as sce
import numpy as np
from ..obj_utils import write_obs

# Wrapper for scrublet allowing text export and filtering

def scrublet(adata, adata_sim=None, filter=False, export_table=None, **kwargs):
    """
    Wrapper function for sce.pp.scrublet(), to allow filtering of resulting object
    """

    # Do we need to read an object with the doublet simulations?

    if adata_sim:
        adata_sim = sc.read(adata_sim)

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
