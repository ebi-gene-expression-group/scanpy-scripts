#!/usr/bin/env python

import logging
import scanpy.api as sc

from scanpy_scripts.wrapper_utils import (
    ScanpyArgParser,
    read_input_object, write_output_object,
)


def main(argv=None):
    argparser = ScanpyArgParser(argv, 'Filter genes by properties and/or simple stats')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_subset_parameters()
    argparser.add_subset_list()
    args = argparser.args

    logging.debug(args)

    adata = read_input_object(args.input_object_file, args.input_format)

    if args.subset_list:
        k = adata.var.index.isin(args.subset_list)
        assert sum(k) > 0, "Unable to proceed as no gene passes filter."
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

    logging.info('Done')
    return 0


if __name__ == '__main__':
    main()
