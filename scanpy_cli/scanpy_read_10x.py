#!/usr/bin/env python3

import logging
import signal
import scanpy.api as sc

from scanpy_cli.wrapper_utils import (
    ScanpyArgParser,
    write_output_object,
)


def main(argv=None):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    argparser = ScanpyArgParser(argv, 'Read10x data for ScanPy')
    argparser.add_argument('-d', '--data-dir',
                           required=True,
                           help='Directory containing the 10X matrix.mtx, '
                                'genes.tsv, and barcodes.tsv files')
    argparser.add_argument('-v', '--var-names',
                           choices=['gene_symbols', 'gene_ids'],
                           default='gene_symbols',
                           help='The variable index')
    argparser.add_output_object()
    args = argparser.args

    logging.debug(args)

    adata = sc.read_10x_mtx(args.data_dir,
                            var_names=args.var_names)

    write_output_object(adata, args.output_object_file, args.output_format)
    return 0


if __name__ == "__main__":
    main()
