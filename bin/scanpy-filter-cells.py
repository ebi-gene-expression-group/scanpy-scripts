#!/usr/bin/env python

import argparse
import sys
from __future__ import print_function


def _validate_subset_thresholds():
    subset_names = args.subset_names.split(",")
    if args.low_thresholds is not None:
        low_thresholds = args.low_thresholds.split(",")
    else:
        print("--low-thresholds should have a comma separated list of numbers"
              " of the same size as --subset-names", file=sys.stderr)
        exit(1)
    if args.high_thresholds is not None:
        high_thresholds = args.high_thresholds.split(",")
    else:
        print("--high-thresholds should have a comma separated list of numbers"
              " of the same size as --subset-names", file=sys.stderr)
        exit(1)
    if len(subset_names) != len(low_thresholds) or len(subset_names) != len(high_thresholds):
        print("--high-thresholds, --low-thresholds and --subset-names should be all"
              " of the same size", file=sys.stderr)
        exit(1)
    return subset_names, low_thresholds, high_thresholds

    # TODO filter by name


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter cells for ScanPy')
    parser.add_argument('--input-object-file', '-i', required=True,
                        help='Path to anndata or loom file')
    parser.add_argument('--input-format', choices=['loom', 'anndata', 'auto-detect'],
                        help='Format for input object. Could be loom or anndata', default="auto-detect")
    parser.add_argument('-o', '--output-object-file', required=True,
                        help='File name in which to store serialized python object.')
    parser.add_argument('-f', '--output-format', choices=['loom', 'anndata'], required=True,
                        help='Format for output. Could be loom or anndata', default="anndata")
    parser.add_argument('-s', '--subset-names',
                        help='Parameters to subset on. Eg, the name of a gene, PC1, '
                             'a column name in anndata object, etc. ')
    parser.add_argument('-l', '--low-thresholds',
                        help='Low cutoffs for the parameters (default is -Inf).')
    parser.add_argument('-j', '--high-thresholds',
                        help='High cutoffs for the parameters (default is Inf).')
    parser.add_argument('-c', '--cells-use',
                        help='Comma-separated list of cell names to use as a subset. Alternatively, '
                             'text file with one cell per line.')

    args = parser.parse_args()

    import scanpy.api as sc

    if args.input_format == "anndata":
        adata = sc.read(args.input_object_file)
    elif args.input_format == "loom":
        from anndata import read_loom

        adata = read_loom(args.input_object_file)

    if args.subset_names is not None:
        names, high_t, low_t = _validate_subset_thresholds()

        from numbers import Number
        for name, h, l in zip(names, high_t, low_t):
            if not isinstance(h, Number):
                print("ERROR: For subset name %s high threshold is not numeric: %s" % (name, h),
                      file=sys.stderr)
            if not isinstance(l, Number):
                print("ERROR: For subset name %s low threshold is not numeric: %s" % (name, l),
                      file=sys.stderr)
            adata = adata[l < adata.obs[name] < h, :]

    if args.cells_use is not None:
        







