#!/usr/bin/env python

from __future__ import print_function
import sys
import signal
import logging
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import argparse
import pandas as pd
from scanpy_wrapper_utils import *


def main(args):
    logging.debug(args)
    import scanpy.api as sc

    if args.input_format == 'anndata':
        adata = sc.read(file_name)
    elif args.input_format == 'loom':
        adata = sc.read_loom(file_name)
    else:
        logging.error('should not reach here')
        sys.exit(1)

    sc.pp.normalize_per_cell(adata, counts_per_cell_after=args.counts_per_cell)

    if args.output_format == 'loom':
        adata.write_loom(args.output_object_file)
    elif args.output_format == 'anndata':
        adata.write(args.output_object_file)

    logging.info('Done')
    return 0


if __name__ == '__main__':
    argparser = ScanpyArgParser('Normalise per-cell quantification data')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.parser.add_argument('-c', '--counts-per-cell',
                                  type=float,
                                  default=1e4,
                                  help='Aimed counts per cell after normalisation, default: 1e4')
    argparser.get_args()

    main(args)
