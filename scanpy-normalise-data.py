#!/usr/bin/env python

from __future__ import print_function
import logging
from scanpy_wrapper_utils import ScanpyArgParser, read_input_object, write_output_object


def main(args):
    logging.debug(args)
    import scanpy.api as sc

    adata = read_input_object(args.input_object_file, args.input_format)

    if args.save_raw:
        adata.raw = sc.pp.log1p(adata, copy=True)

    sc.pp.normalize_per_cell(adata, counts_per_cell_after=args.scale_factor)

    write_output_object(adata, args.output_object_file, args.output_format)

    print(adata)
    logging.info('Done')
    return 0


if __name__ == '__main__':
    argparser = ScanpyArgParser('Normalise per-cell quantification data')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_argument('-s', '--scale-factor',
                           type=float,
                           default=1e4,
                           help='Aimed counts per cell after normalisation, default: 1e4')
    argparser.add_argument('-r', '--save-raw',
                           action='store_true',
                           default=False,
                           help='Save raw quantification in log scale before normalisation.')
    args = argparser.get_args()

    main(args)
