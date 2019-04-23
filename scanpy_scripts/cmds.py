"""
Provides sub-commands
"""

import scanpy as sc

from .cmd_utils import make_subcmd
from .cmd_options import (
    READ_CMD_OPTIONS,
    FILTER_CMD_OPTIONS,
    NORM_CMD_OPTIONS,
    HVG_CMD_OPTIONS,
    SCALE_CMD_OPTIONS,
    REGRESS_CMD_OPTIONS,
    PCA_CMD_OPTIONS,
    NEIGHBOR_CMD_OPTIONS,
    UMAP_CMD_OPTIONS,
    TSNE_CMD_OPTIONS,
)
from .lib._read import read_10x
from .lib._filter import filter_anndata
from .lib._norm import normalize
from .lib._hvg import hvg
from .lib._neighbors import neighbors
from .lib._umap import umap
from .lib._tsne import tsne


_I_DESC = '<input_obj>:   input file in format specfied by --input-format'
_O_DESC = '<output_obj>:  output file in format specfied by --output-format'
_IO_DESC = '\n'.join([_I_DESC, _O_DESC])


READ_CMD = make_subcmd(
    'read',
    READ_CMD_OPTIONS,
    read_10x,
    cmd_desc='Read 10x data and save in specified format.',
    arg_desc=_O_DESC,
)


FILTER_CMD = make_subcmd(
    'filter',
    FILTER_CMD_OPTIONS,
    filter_anndata,
    cmd_desc='Filter data based on specified conditions.',
    arg_desc=_IO_DESC,
)


NORM_CMD = make_subcmd(
    'norm',
    NORM_CMD_OPTIONS,
    normalize,
    cmd_desc='Normalise data per cell.',
    arg_desc=_IO_DESC,
)


HVG_CMD = make_subcmd(
    'hvg',
    HVG_CMD_OPTIONS,
    hvg,
    cmd_desc='Find highly variable genes.',
    arg_desc=_IO_DESC,
)


SCALE_CMD = make_subcmd(
    'scale',
    SCALE_CMD_OPTIONS,
    sc.pp.scale,
    cmd_desc='Scale data per gene.',
    arg_desc=_IO_DESC,
)


REGRESS_CMD = make_subcmd(
    'regress',
    REGRESS_CMD_OPTIONS,
    sc.pp.regress_out,
    cmd_desc='Regress-out observation variables.',
    arg_desc=_IO_DESC,
)


PCA_CMD = make_subcmd(
    'pca',
    PCA_CMD_OPTIONS,
    sc.pp.pca,
    cmd_desc='Dimensionality reduction by PCA.',
    arg_desc=_IO_DESC,
)

NEIGHBOR_CMD = make_subcmd(
    'neighbor',
    NEIGHBOR_CMD_OPTIONS,
    neighbors,
    cmd_desc='Compute a neighbourhood graph of observations.',
    arg_desc=_IO_DESC,
)

UMAP_CMD = make_subcmd(
    'umap',
    UMAP_CMD_OPTIONS,
    umap,
    cmd_desc='Embed the neighborhood graph using UMAP.',
    arg_desc=_IO_DESC,
)

TSNE_CMD = make_subcmd(
    'tsne',
    TSNE_CMD_OPTIONS,
    tsne,
    cmd_desc='Embed the cells using t-SNE.',
    arg_desc=_IO_DESC,
)
