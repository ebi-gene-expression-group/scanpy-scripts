#!/usr/bin/env python

import argparse
import sys
from __future__ import print_function

parser = argparse.ArgumentParser(description='Read10x data for ScanPy')
parser.add_argument('--data-dir', '-d', dest='data-dir', required=True,
                    help='Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv files provided by 10X.')
parser.add_argument('-o', '--output-object-file', required=True,
                    help='File name in which to store serialized python object.')
parser.add_argument('-f', '--output-format', choices=['loom','anndata'],
                    help='Format for output. Could be loom or anndata', default="anndata")

args = parser.parse_args()


import pandas as pd
import scanpy.api as sc

adata = sc.read(args.data-dir + 'matrix.mtx', cache=True).T  # transpose the data
genes = pd.read_csv(args.data-dir + 'genes.tsv', header=None, sep='\t')
adata.var_names = genes[1]
adata.var['gene_ids'] = genes[0]  # add the gene ids as annotation of the variables/genes
adata.obs_names = pd.read_csv(args.data-dir + 'barcodes.tsv', header=None)[0]
adata.var_names_make_unique()

if args.output-format == "adata":
    adata.write(args.output-object-file)
elif args.output-format == "loom":
    adata.write_loom(args.output-object-file)
else:
    print("We shouldn't be here, somehow you managed to use a file format that "
          "passed checks but is not implemented", file=sys.stderr )
    exit(1)
