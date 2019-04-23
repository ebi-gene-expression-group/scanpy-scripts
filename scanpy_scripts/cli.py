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
cli.add_command(UMAP_CMD)
cli.add_command(TSNE_CMD)
