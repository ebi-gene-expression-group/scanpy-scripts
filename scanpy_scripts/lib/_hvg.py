"""
scanpy hvg
"""

import numpy as np
import scanpy as sc


def hvg(
    adata,
    mean_limits=(0.0125, 3),
    disp_limits=(0.5, float("inf")),
    **kwargs,
):
    """
    Wrapper function for sc.highly_variable_genes()
    """

    # Check for n_top_genes beeing greater than the total genes

    if "n_top_genes" in kwargs and kwargs["n_top_genes"] is not None:
        kwargs["n_top_genes"] = min(adata.n_vars, kwargs["n_top_genes"])

    always_hv_genes = None
    if "always_hv_genes_file" in kwargs and kwargs["always_hv_genes_file"] is not None:
        with open(kwargs["always_hv_genes_file"], "r") as f:
            always_hv_genes = f.read().splitlines()

    never_hv_genes = None
    if "never_hv_genes_file" in kwargs and kwargs["never_hv_genes_file"] is not None:
        with open(kwargs["never_hv_genes_file"], "r") as f:
            never_hv_genes = f.read().splitlines()

    # to avoid upsetting the scanpy function with unexpected keyword arguments
    del kwargs["always_hv_genes_file"]
    del kwargs["never_hv_genes_file"]

    sc.pp.highly_variable_genes(
        adata,
        min_mean=mean_limits[0],
        max_mean=mean_limits[1],
        min_disp=disp_limits[0],
        max_disp=disp_limits[1],
        **kwargs,
    )

    return switch_hvgs(adata, always_hv_genes, never_hv_genes)


def switch_hvgs(adata, always_hv_genes=None, never_hv_genes=None):
    """
    Function to switch on/off highly variable genes based on a list of genes.

    >>> adata = sc.datasets.pbmc3k()
    >>> sc.pp.normalize_total(adata)
    >>> sc.pp.log1p(adata)
    >>> sc.pp.highly_variable_genes(adata)
    >>> adata = switch_hvgs(adata, always_hv_genes=['MIR1302-10', 'FAM138A'], never_hv_genes=['ISG15', 'TNFRSF4'])
    >>> adata.var.loc['ISG15'].highly_variable
    False
    >>> adata.var.loc['TNFRSF4'].highly_variable
    False
    >>> adata.var.loc['MIR1302-10'].highly_variable
    True
    >>> adata.var.loc['CPSF3L'].highly_variable
    True
    """
    if always_hv_genes is not None:
        adata.var.highly_variable = np.logical_or(
            adata.var.highly_variable, adata.var_names.isin(always_hv_genes)
        )

    if never_hv_genes is not None:
        adata.var.highly_variable = np.logical_and(
            adata.var.highly_variable, ~adata.var_names.isin(never_hv_genes)
        )

    return adata
