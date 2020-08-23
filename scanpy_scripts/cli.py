"""
scanpy
"""

import logging
import click
import scanpy as sc
from .click_utils import NaturalOrderGroup
from .cmds import (
    READ_CMD,
    FILTER_CMD,
    NORM_CMD,
    HVG_CMD,
    SCALE_CMD,
    REGRESS_CMD,
    PCA_CMD,
    NEIGHBOR_CMD,
    UMAP_CMD,
    TSNE_CMD,
    FDG_CMD,
    LOUVAIN_CMD,
    LEIDEN_CMD,
    DIFFEXP_CMD,
    PAGA_CMD,
    DIFFMAP_CMD,
    DPT_CMD,
    PLOT_EMBED_CMD,
    PLOT_PAGA_CMD,
    PLOT_STACKED_VIOLIN_CMD,
    PLOT_DOT_CMD,
    PLOT_MATRIX_CMD,
    PLOT_HEATMAP_CMD,
    HARMONY_INTEGRATE_CMD,
    BBKNN_CMD,
    MNN_CORRECT_CMD,
    COMBAT_CMD,
)


@click.group(cls=NaturalOrderGroup)
@click.option(
    '--debug',
    is_flag=True,
    default=False,
    help='Print debug information',
)
@click.option(
    '--verbosity',
    type=click.INT,
    default=3,
    help='Set scanpy verbosity',
)
@click.version_option(
    version='0.2.0',
    prog_name='scanpy',
)
def cli(debug=False, verbosity=3):
    """
    Command line interface to [scanpy](https://github.com/theislab/scanpy)
    """
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format=('%(asctime)s; %(levelname)s; %(filename)s; '
                '%(funcName)s(): %(message)s'),
        datefmt='%y-%m-%d %H:%M:%S',
    )
    logging.debug('debugging')
    sc.settings.verbosity = verbosity
    return 0


cli.add_command(READ_CMD)
cli.add_command(FILTER_CMD)
cli.add_command(NORM_CMD)
cli.add_command(HVG_CMD)
cli.add_command(SCALE_CMD)
cli.add_command(REGRESS_CMD)
cli.add_command(PCA_CMD)
cli.add_command(NEIGHBOR_CMD)


@cli.group(cls=NaturalOrderGroup)
def embed():
    """Embed cells into two-dimensional space."""


embed.add_command(UMAP_CMD)
embed.add_command(TSNE_CMD)
embed.add_command(FDG_CMD)
embed.add_command(DIFFMAP_CMD)


@cli.group(cls=NaturalOrderGroup)
def cluster():
    """Cluster cells into sub-populations."""


cluster.add_command(LOUVAIN_CMD)
cluster.add_command(LEIDEN_CMD)


cli.add_command(DIFFEXP_CMD)
cli.add_command(PAGA_CMD)
cli.add_command(DPT_CMD)


@cli.group(cls=NaturalOrderGroup)
def integrate():
    """Integrate cells from different experimental batches."""

integrate.add_command(HARMONY_INTEGRATE_CMD)
integrate.add_command(BBKNN_CMD)
integrate.add_command(MNN_CORRECT_CMD)
integrate.add_command(COMBAT_CMD)


@cli.group(cls=NaturalOrderGroup)
def plot():
    """Visualise data."""


plot.add_command(PLOT_EMBED_CMD)
plot.add_command(PLOT_PAGA_CMD)
plot.add_command(PLOT_STACKED_VIOLIN_CMD)
plot.add_command(PLOT_DOT_CMD)
plot.add_command(PLOT_MATRIX_CMD)
plot.add_command(PLOT_HEATMAP_CMD)
