#!/usr/bin/env python

from __future__ import print_function
import sys
import signal
import logging
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import argparse
import pandas as pd
import os.path


def parse_cmd_arguments():
    parser = argparse.ArgumentParser(description='Read10x data for ScanPy')
    parser.add_argument('-d', '--data-dir', required=True,
                        help='Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv files provided by 10X.')
    parser.add_argument('-o', '--output-object-file', required=True,
                        help='File name in which to store serialized python object.')
    parser.add_argument('-f', '--output-format', choices=['loom','anndata'],
                        help='Format for output. Could be loom or anndata', default="anndata")
    parser.add_argument('--debug', action='store_true',
                        help='print debug information')
    args = parser.parse_args()
    return args


def main(args):
    logging.debug(args)
    import scanpy.api as sc
    adata = sc.read(os.path.join(args.data_dir, 'matrix.mtx'), cache=True).T  # transpose the data
    genes = pd.read_csv(os.path.join(args.data_dir, 'genes.tsv'), header=None, sep='\t')
    adata.var_names = genes[0]
    adata.var['gene_names'] = genes[1].values  # add the gene ids as annotation of the variables/genes
    cells = pd.read_csv(os.path.join(args.data_dir, 'barcodes.tsv'), header=None)
    adata.obs_names = cells[0]
    adata.var_names_make_unique()

    if args.output_format == "anndata":
        adata.write(args.output_object_file)
    elif args.output_format == "loom":
        adata.write_loom(args.output_object_file)
    else:
        print("We shouldn't be here, somehow you managed to use a file format that "
              "passed checks but is not implemented", file=sys.stderr )
        sys.exit(1)
    return 0


if __name__ == "__main__":
    args = parse_cmd_arguments()

    if args.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.WARN
    logging.basicConfig(
            level=log_level,
            format='%(asctime)s; %(levelname)s; %(filename)s; %(funcName)s(): %(message)s',
            datefmt='%y-%m-%d %H:%M:%S')

    main(args)
