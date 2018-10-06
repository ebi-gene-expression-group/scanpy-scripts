#!/usr/bin/env python

from __future__ import print_function
import signal
import logging
from scanpy_wrapper_utils import ScanpyArgParser, read_input_object, write_output_object
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def main(args):
    logging.debug(args)
    import scanpy.api as sc

    adata = read_input_object(args.input_object_file, args.input_format)

    sc.pp.normalize_per_cell(adata, counts_per_cell_after=args.counts_per_cell)

    write_output_object(adata, args.output_object_file, args.output_format)

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
    args = argparser.get_args()

    main(args)
