#!/usr/bin/env python

from __future__ import print_function
import sys
import signal
import logging
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import argparse
import pandas as pd


def parse_cmd_arguments():
    cs_str = generate_comma_separated_list('string')
    cs_float = generate_comma_separated_list('numeric')
    cs_names = cs_str('subset-names')
    cs_high = cs_float('high-thresholds')
    cs_low = cs_float('low-thresholds')

    parser = argparse.ArgumentParser(description='Filter genes for ScanPy')

    parser.add_argument('-i', '--input-object-file', required=True,
                        help='Path to anndata or loom file')
    parser.add_argument('--input-format', choices=['loom', 'anndata', 'auto-detect'],
                        help='Format for input object. Could be loom or anndata', default='auto-detect')
    parser.add_argument('-o', '--output-object-file', required=True,
                        help='File name in which to store serialized python object')
    parser.add_argument('-f', '--output-format', choices=['loom', 'anndata'], required=True,
                        help='Format for output object. Could be loom or anndata', default='anndata')
    parser.add_argument('-s', '--subset-names', type=cs_names, default=[],
                        help='Parameters to subset on. Eg, the name of a gene, PC1, '
                        'a column name in anndata object, etc. ')
    parser.add_argument('-l', '--low-thresholds', type=cs_low, default=[],
                        help='Low cutoffs for the parameters (default is -Inf).')
    parser.add_argument('-j', '--high-thresholds', type=cs_high, default=[],
                        help='High cutoffs for the parameters (default is Inf).')
    parser.add_argument('-g', '--genes-use',
                        help='Comma-separated list of gene names to use as a subset. Alternatively, '
                        'text file with one gene per line.')
    parser.add_argument('--debug', action='store_true',
                        help='Print debug information')
    args = parser.parse_args()
    return args


def generate_comma_separated_list(dtyp):
    supported_dtypes = {'integer':int, 'numeric':float, 'string':str}
    if dtyp not in supported_dtypes:
        logging.warn('unsupported data type: "{}", default to "string"'.format(dtyp))
        dtype_name = 'string'
        dtype = str
    else:
        dtype_name = dtyp
        dtype = supported_dtypes[dtyp]

    def generate_named_comma_separated_list(arg_name, dtype=dtype, dtype_name=dtype_name):
        def parse_comma_separated_list(arg, arg_name=arg_name, dtype=dtype, dtype_name=dtype_name):
            L = arg.split(',')
            V = []
            for s in L:
                try:
                    v = dtype(s)
                except ValueError:
                    logging.error('--{} expects comma separated list of [{}] but received "{}"'.format(arg_name, dtype_name, s))
                    break
                V.append(v)
            return V
        return parse_comma_separated_list
    return generate_named_comma_separated_list


def main(args):
    logging.debug(args)
    import scanpy.api as sc

    if len(args.subset_names) != len(args.low_thresholds):
        logging.error('--low-thresholds should have a comma separated list of numerics of the same size as --subset-names')
        return 1

    if len(args.subset_names) != len(args.high_thresholds):
        logging.error('--high-thresholds should have a comma separated list of numerics of the same size as --subset-names')
        return 1

    if args.input_format == 'auto-detect':
        if args.input_object_file.endswith('.loom'):
            args.input_format = 'loom'
        elif args.input_object_file.endswith('.h5ad'):
            args.input_format = 'anndata'
        else:
            logging.error('unknown input format: "{}", please check suffix is either ".loom" or ".h5ad"'.format(args.input_object_file))
            return 1

    if args.input_format == 'anndata':
        adata = sc.read(args.input_object_file)
    elif args.input_format == 'loom':
        adata = sc.read_loom(args.input_object_file)
    else:
        logging.error('should not reach here')
        return 1

    if args.genes_use is not None:
        genes_to_use = list(pd.read_table(args.genes_use, header=None).iloc[:,0])
        logging.debug(len(genes_to_use))
        adata = adata[adata.var_names.isin(genes_to_use)]

    for name,h,l in zip(args.subset_names, args.high_thresholds, args.low_thresholds):
        if name not in adata.var.columns:
            logging.warn('subset-name "{}" not present in data, omitted'.format(name))
            continue
        adata = adata[(adata.var[name] < h) & (adata.var[name] > l), :]

    if args.output_format == 'loom':
        adata.write_loom(args.output_object_file)
    elif args.output_format == 'anndata':
        adata.write(args.output_object_file)

    return 0


if __name__ == '__main__':
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
