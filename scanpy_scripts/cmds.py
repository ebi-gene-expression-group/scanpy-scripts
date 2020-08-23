"""
Provides sub-commands
"""

import os
import sys
import scanpy as sc
import scanpy.external as sce

from .cmd_utils import (
    make_subcmd,
    make_plot_function,
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
from .lib._bbknn import bbknn
from .lib._mnn import mnn_correct
from .lib._combat import combat

LANG = os.environ.get('LANG', None)

if LANG is None or not (LANG.endswith('UTF-8') or
                        LANG.endswith('UTF8') or
                        LANG.endswith('utf-8') or
                        LANG.endswith('utf8')):
    print('This programme requires a UTF-8 locale, please check your $LANG setting.')
    sys.exit(0)


_I_DESC = '<input_obj>:   input file in format specfied by --input-format'
_O_DESC = '<output_obj>:  output file in format specfied by --output-format'
_P_DESC = '<output_fig>:  output figure in pdf or png format'
_IO_DESC = '\n'.join([_I_DESC, _O_DESC])
_IP_DESC = '\n'.join([_I_DESC, _P_DESC])


READ_CMD = make_subcmd(
    'read',
    read_10x,
    cmd_desc='Read 10x data and save in specified format.',
    arg_desc=_O_DESC,
)


FILTER_CMD = make_subcmd(
    'filter',
    filter_anndata,
    cmd_desc='Filter data based on specified conditions.',
    arg_desc=_IO_DESC,
)


NORM_CMD = make_subcmd(
    'norm',
    normalize,
    cmd_desc='Normalise data per cell.',
    arg_desc=_IO_DESC,
)


HVG_CMD = make_subcmd(
    'hvg',
    hvg,
    cmd_desc='Find highly variable genes.',
    arg_desc=_IO_DESC,
)


SCALE_CMD = make_subcmd(
    'scale',
    sc.pp.scale,
    cmd_desc='Scale data per gene.',
    arg_desc=_IO_DESC,
)


REGRESS_CMD = make_subcmd(
    'regress',
    sc.pp.regress_out,
    cmd_desc='Regress-out observation variables.',
    arg_desc=_IO_DESC,
)


PCA_CMD = make_subcmd(
    'pca',
    pca,
    cmd_desc='Dimensionality reduction by PCA.',
    arg_desc=_IO_DESC,
)

NEIGHBOR_CMD = make_subcmd(
    'neighbor',
    neighbors,
    cmd_desc='Compute a neighbourhood graph of observations.',
    arg_desc=_IO_DESC,
)

UMAP_CMD = make_subcmd(
    'umap',
    umap,
    cmd_desc='Embed the neighborhood graph using UMAP.',
    arg_desc=_IO_DESC,
)

TSNE_CMD = make_subcmd(
    'tsne',
    tsne,
    cmd_desc='Embed the cells using t-SNE.',
    arg_desc=_IO_DESC,
)

FDG_CMD = make_subcmd(
    'fdg',
    fdg,
    cmd_desc='Embed the neighborhood graph using force-directed graph.',
    arg_desc=_IO_DESC,
)

DIFFMAP_CMD = make_subcmd(
    'diffmap',
    diffmap,
    cmd_desc='Embed the neighborhood graph using diffusion map.',
    arg_desc=_IO_DESC,
)

LOUVAIN_CMD = make_subcmd(
    'louvain',
    louvain,
    cmd_desc='Find clusters by Louvain algorithm.',
    arg_desc=_IO_DESC,
)

LEIDEN_CMD = make_subcmd(
    'leiden',
    leiden,
    cmd_desc='Find clusters by Leiden algorithm.',
    arg_desc=_IO_DESC,
)

DIFFEXP_CMD = make_subcmd(
    'diffexp',
    diffexp,
    cmd_desc='Find markers for each clusters.',
    arg_desc=_IO_DESC,
)

PAGA_CMD = make_subcmd(
    'paga',
    paga,
    cmd_desc='Trajectory inference by abstract graph analysis.',
    arg_desc=_IO_DESC,
)

DPT_CMD = make_subcmd(
    'dpt',
    dpt,
    cmd_desc='Calculate diffusion pseudotime relative to the root cells.',
    arg_desc=_IO_DESC,
)

PLOT_EMBED_CMD = make_subcmd(
    'embed',
    make_plot_function('embedding'),
    cmd_desc='Plot cell embeddings.',
    arg_desc=_IP_DESC,
)

PLOT_STACKED_VIOLIN_CMD = make_subcmd(
    'sviol',
    make_plot_function('sviol'),
    cmd_desc='Plot stacked violin plots.',
    arg_desc=_IP_DESC,
)

PLOT_DOT_CMD = make_subcmd(
    'dot',
    make_plot_function('dot'),
    cmd_desc='Plot a dot plot of expression values.',
    arg_desc=_IP_DESC,
)

PLOT_MATRIX_CMD = make_subcmd(
    'matrix',
    make_plot_function('matrix'),
    cmd_desc='Plot a heatmap of the mean expression values per cluster.',
    arg_desc=_IP_DESC,
)

PLOT_HEATMAP_CMD = make_subcmd(
    'heat',
    make_plot_function('heat'),
    cmd_desc='Plot a heatmap of the expression values of genes.',
    arg_desc=_IP_DESC,
)

PLOT_PAGA_CMD = make_subcmd(
    'paga',
    make_plot_function('plot_paga', kind='paga'),
    cmd_desc='Plot PAGA trajectories.',
    arg_desc=_IP_DESC,
    opt_set='plot_paga'
)

COMBAT_CMD = make_subcmd(
    'combat',
    combat,
    cmd_desc='ComBat function for batch effect correction',
    arg_desc=_IO_DESC
)

HARMONY_INTEGRATE_CMD = make_subcmd(
    'harmony',
    sce.pp.harmony_integrate,
    cmd_desc='Use harmonypy [Korunsky19] to integrate different experiments.',    
    arg_desc=_IO_DESC,
)

BBKNN_CMD = make_subcmd(
    'bbknn',
    bbknn,
    cmd_desc='Batch balanced kNN [Polanski19].',
    arg_desc=_IO_DESC,
)

MNN_CORRECT_CMD = make_subcmd(
    'mnn',
    mnn_correct,
    cmd_desc='Correct batch effects by matching mutual nearest neighbors [Haghverdi18] [Kang18].',
    arg_desc=_IO_DESC,
)
