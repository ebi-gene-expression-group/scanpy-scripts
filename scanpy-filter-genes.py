#!/usr/bin/env python

from __future__ import print_function
import logging
import scanpy.api as sc
from scanpy_wrapper_utils import ScanpyArgParser, read_input_object, write_output_object


def main(args):
    logging.debug(args)

    adata = read_input_object(args.input_object_file, args.input_format)

    if args.subset_list:
        k = adata.var.index.isin(args.subset_list)
        if sum(k) == 0:
            logging.error("Unable to proceed as no gene passes filter.")
            return 1
        adata = adata[:, k]

    inf, neg_inf = float('Inf'), float('-Inf')
    for name, high, low in zip(args.parameter_names, args.high_thresholds, args.low_thresholds):
        if name == 'n_cells':
            if low > neg_inf:
                sc.pp.filter_genes(adata, min_cells=low)
            if high < inf:
                sc.pp.filter_genes(adata, max_cells=high)
        elif name == 'n_counts':
            if low > neg_inf:
                sc.pp.filter_genes(adata, min_counts=low)
            if high < inf:
                sc.pp.filter_genes(adata, max_counts=high)
        elif name not in adata.var.columns:
            msg = 'parameter-name "{}" not present in data, omitted'.format(name)
            logging.warning(msg)
        else:
            adata = adata[(adata.var[name] < high) & (adata.var[name] > low), :]

    write_output_object(adata, args.output_object_file, args.output_format)

    print(adata)
    logging.info('Done')
    return 0


if __name__ == '__main__':
    argparser = ScanpyArgParser('Filter genes by properties and/or simple stats')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_subset_parameters()
    argparser.add_subset_list()
    args = argparser.get_args()

    main(args)
