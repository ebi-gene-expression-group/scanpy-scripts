#!/usr/bin/env python

from __future__ import print_function
import os.path
import logging
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import pandas as pd
from scanpy_wrapper_utils import ScanpyArgParser, write_output_object


def main(args):
    logging.debug(args)
    import scanpy.api as sc

    adata = sc.read(os.path.join(args.data_dir, 'matrix.mtx'), cache=True).T  # transpose the data
    genes = pd.read_csv(os.path.join(args.data_dir, 'genes.tsv'), header=None, sep='\t')
    cells = pd.read_csv(os.path.join(args.data_dir, 'barcodes.tsv'), header=None)

    adata.var_names = genes[0]
    adata.var['gene_names'] = genes[1].values  # add the gene ids as annotation of the variables/genes
    adata.var_names_make_unique()

    adata.obs_names = cells[0]

    write_output_object(adata, args.output_object_file, args.output_format)
    return 0


if __name__ == "__main__":
    argparser = ScanpyArgParser('Read10x data for ScanPy')
    argparser.parser.add_argument('-d', '--data-dir',
                                  required=True,
                                  help='Directory containing the matrix.mtx, genes.tsv, '
                                       'and barcodes.tsv files provided by 10X.')
    argparser.add_output_object()
    args = argparser.get_args()

    main(args)
