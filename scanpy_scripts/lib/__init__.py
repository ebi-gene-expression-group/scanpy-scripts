"""
Provides exported functions
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from ._read import read_10x
from ._filter import filter_anndata
from ._norm import normalize
from ._hvg import hvg
from ._neighbors import neighbors


def expression_colormap():
    """Returns a nice color map for highlighting gene expression
    """
    reds = plt.cm.Reds(np.linspace(0, 1, 128))
    greys = plt.cm.Greys_r(np.linspace(0.7, 0.8, 20))
    palette = np.vstack([greys, reds])
    return LinearSegmentedColormap.from_list('expression', palette)


def cross_table(adata, x, y, normalise=None, highlight=False):
    """Make a cross table comparing two categorical annotations
    """
    x_attr = adata.obs[x]
    y_attr = adata.obs[y]
    assert not _is_numeric(x_attr.values), "Can not operate on numeric variable {}".format(x)
    assert not _is_numeric(y_attr.values), "Can not operate on numeric variable {}".format(y)
    crs_tbl = pd.crosstab(x, y)
    if normalise == 'x':
        x_sizes = x_attr.groupby(x_attr).size().values
        crs_tbl = (crs_tbl.T / x_sizes * 100).round(2).T
    elif normalise == 'y':
        y_sizes = x_attr.groupby(y_attr).size().values
        crs_tbl = (crs_tbl / y_sizes * 100).round(2)
    if highlight:
        return crs_tbl.style.background_gradient(cmap='viridis', axis=0)
    return crs_tbl


def run_harmony():
    # import subprocess as sbp
    pass


def _is_numeric(x):
    return x.dtype.kind in ('i', 'f')


def pseudo_bulk(adata, groupby, FUN=np.mean):
    """Make pseudo bulk data from grouped sc data
    """
    group_attr = adata.obs[groupby].astype(str).values
    groups = np.unique(group_attr)
    n_level = len(groups)
    summarised = np.zeros((n_level, adata.X.shape[1]))
    for i, grp in enumerate(groups):
        k_grp = group_attr == grp
        summarised[i] = FUN(adata.X[k_grp, :], axis=0, keepdims=True)
    return summarised
