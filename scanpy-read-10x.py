#!/usr/bin/env python

from __future__ import print_function
import logging
import signal
from scanpy_wrapper_utils import ScanpyArgParser, write_output_object
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def main(args):
    logging.debug(args)
    import scanpy.api as sc

    adata = sc.read_10x_mtx(args.data_dir,
                            var_names=args.var_names)

    write_output_object(adata, args.output_object_file, args.output_format)

    print(adata)
    logging.info('Done')
    return 0


if __name__ == "__main__":
    argparser = ScanpyArgParser('Read10x data for ScanPy')
    argparser.add_argument('-d', '--data-dir',
                           required=True,
                           help='Directory containing the 10X matrix.mtx, '
                                'genes.tsv, and barcodes.tsv files')
    argparser.add_argument('-v', '--var-names',
                           choices=['gene_symbols', 'gene_ids'],
                           default='gene_symbols',
                           help='The variable index')
    argparser.add_output_object()
    args = argparser.get_args()

    main(args)
