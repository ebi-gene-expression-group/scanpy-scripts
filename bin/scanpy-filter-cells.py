#!/usr/bin/env python

from __future__ import print_function
import logging
import scanpy.api as sc
from scanpy_wrapper_utils import ScanpyArgParser, read_input_object, write_output_object


def main(args):
    logging.debug(args)

    adata = read_input_object(args.input_object_file, args.input_format)

    if args.subset_list:
        adata = adata[:, adata.obs.isin(args.subset_list)]

    inf, neg_inf = float('Inf'), float('-Inf')
    for name, high, low in zip(args.parameter_names, args.high_thresholds, args.low_thresholds):
        if name == 'n_genes':
            if low > neg_inf:
                sc.pp.filter_cells(adata, min_genes=low)
            if high < neg_inf:
                sc.pp.filter_cells(adata, max_genes=high)
        elif name == 'n_counts':
            if low > neg_inf:
                sc.pp.filter_cells(adata, min_counts=low)
            if high < neg_inf:
                sc.pp.filter_cells(adata, max_counts=high)
        elif name not in adata.obs.columns:
            msg = 'parameter-name "{}" not present in data, omitted'.format(name)
            logging.warning(msg)
        else:
            adata = adata[(adata.obs[name] < high) & (adata.obs[name] > low), :]

    write_output_object(adata, args.output_object_file, args.output_format)

    logging.info('Done')
    return 0


if __name__ == "__main__":
    argparser = ScanpyArgParser('Filter genes for ScanPy')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_subset_parameters()
    argparser.add_subset_list()
    args = argparser.get_args()

    main(args)
