"""
Provides sub-commands
"""

import scanpy as sc

from .cmd_utils import (
    make_subcmd,
    plot_embeddings,
)
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
    FDG_CMD_OPTIONS,
    LOUVAIN_CMD_OPTIONS,
    LEIDEN_CMD_OPTIONS,
    DIFFEXP_CMD_OPTIONS,
    PAGA_CMD_OPTIONS,
    DIFFMAP_CMD_OPTIONS,
    DPT_CMD_OPTIONS,
    PLOT_EMBED_CMD_OPTIONS,
)
from .lib._read import read_10x
from .lib._filter import filter_anndata
from .lib._norm import normalize
from .lib._hvg import hvg
from .lib._pca import pca
from .lib._neighbors import neighbors
from .lib._umap import umap
from .lib._tsne import tsne
from .lib._fdg import fdg
from .lib._louvain import louvain
from .lib._leiden import leiden
from .lib._diffexp import diffexp
from .lib._paga import paga
from .lib._diffmap import diffmap
from .lib._dpt import dpt


_I_DESC = '<input_obj>:   input file in format specfied by --input-format'
_O_DESC = '<output_obj>:  output file in format specfied by --output-format'
_P_DESC = '<output_fig>:  output figure in pdf or png format'
_IO_DESC = '\n'.join([_I_DESC, _O_DESC])
_IP_DESC = '\n'.join([_I_DESC, _P_DESC])


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
    pca,
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

FDG_CMD = make_subcmd(
    'fdg',
    FDG_CMD_OPTIONS,
    fdg,
    cmd_desc='Embed the neighborhood graph using force-directed graph.',
    arg_desc=_IO_DESC,
)

DIFFMAP_CMD = make_subcmd(
    'diffmap',
    DIFFMAP_CMD_OPTIONS,
    diffmap,
    cmd_desc='Embed the neighborhood graph using diffusion map.',
    arg_desc=_IO_DESC,
)

LOUVAIN_CMD = make_subcmd(
    'louvain',
    LOUVAIN_CMD_OPTIONS,
    louvain,
    cmd_desc='Find clusters by Louvain algorithm.',
    arg_desc=_IO_DESC,
)

LEIDEN_CMD = make_subcmd(
    'leiden',
    LEIDEN_CMD_OPTIONS,
    leiden,
    cmd_desc='Find clusters by Leiden algorithm.',
    arg_desc=_IO_DESC,
)

DIFFEXP_CMD = make_subcmd(
    'diffexp',
    DIFFEXP_CMD_OPTIONS,
    diffexp,
    cmd_desc='Find markers for each clusters.',
    arg_desc=_IO_DESC,
)

PAGA_CMD = make_subcmd(
    'paga',
    PAGA_CMD_OPTIONS,
    paga,
    cmd_desc='Trajectory inference by abstract graph analysis.',
    arg_desc=_IO_DESC,
)

DPT_CMD = make_subcmd(
    'dpt',
    DPT_CMD_OPTIONS,
    dpt,
    cmd_desc='Calculate diffusion pseudotime relative to the root cells.',
    arg_desc=_IO_DESC,
)

PLOT_EMBED_CMD = make_subcmd(
    'embed',
    PLOT_EMBED_CMD_OPTIONS,
    plot_embeddings,
    cmd_desc='Plot cell embeddings.',
    arg_desc=_IP_DESC,
)
