"""
scanpy external scrublet
"""

import scanpy.external as sce
import numpy as np
from ..obj_utils import write_obs
import matplotlib
import matplotlib.pyplot as plt

# Wrapper for bbknn allowing use of non-standard slot


def scrublet(adata, filter=False, export_table=None, **kwargs):
    """
    Wrapper function for sce.pp.scrublet(), to allow filtering of resulting object
    """

    sce.pp.scrublet(adata, **kwargs)

    # Do any export before optional filtering

    if export_table:
        write_obs(adata, ["doublet_score", "predicted_doublet"], export_table)

    # Filter out predited doublets

    if filter:
        adata._inplace_subset_obs(np.invert(adata.obs["predicted_doublet"]))

    return adata


# Just absorb the extra plotting args before passing to
# scanpy.external.pl.scrublet_score_distribution


def plot_scrublet(
    adata,
    scale_hist_obs="log",
    scale_hist_sim="linear",
    fig_size=(8, 3),
    dpi=40,
    save=None,
    **kwargs
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
    )
    print("Writing scrublet plot to %s" % save)
    plt.savefig(save, dpi=40)
