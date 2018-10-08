#!/usr/bin/env python

from __future__ import print_function
import signal
import logging
import pandas as pd
import scanpy.api as sc
from scanpy_wrapper_utils import ScanpyArgParser, read_input_object, write_output_object
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def read_subset_items(input_string):
    if ',' in input_string:
        list_parser = comma_separated_list('genes-use', str)
        return list_parser(input_string)
    else:
        return list(pd.read_table(input_string, header=None).iloc[:, 0].values)


def main(args):
    logging.debug(args)

    adata = read_input_object(args.input_object_file, args.input_format)

    if args.genes_use is not None:
        genes_to_use = read_subset_items(args.genes_use)
        adata = adata[:, adata.var_names.isin(genes_to_use)]

    inf, neg_inf = float('Inf'), float('-Inf')
    for name, h, l in zip(args.subset_names, args.high_thresholds, args.low_thresholds):
        if name == 'n_cells':
            if l > neg_inf:
                sc.pp.filter_genes(adata, min_cells=l)
            if h < inf:
                sc.pp.filter_genes(adata, max_cells=h)
        elif name == 'n_counts':
            if l > neg_inf:
                sc.pp.filter_genes(adata, min_counts=l)
            if h < inf:
                sc.pp.filter_genes(adata, max_counts=h)
        elif name not in adata.var.columns:
            logging.warning('subset-name "{}" not present in data, omitted'.format(name))
        else:
            adata = adata[(adata.var[name] < h) & (adata.var[name] > l), :]

    write_output_object(adata, args.output_object_file, args.output_format)

    logging.info('Done')
    return 0


if __name__ == '__main__':
    argparser = ScanpyArgParser('Filter genes for ScanPy')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_subset_parameters()
    argparser.parser.add_argument('-g', '--genes-use',
                                  help='Comma-separated list of gene names to use as a subset. '
                                       'Alternatively, text file with one gene per line.')
    args = argparser.get_args()

    main(args)
