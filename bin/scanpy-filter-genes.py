#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
import logging
logging.basicConfig(
    level=logging.WARN,
    format='%(asctime)s; %(levelname)s; %(filename)s; %(funcName)s(): %(message)s',
    datefmt='%y-%m-%d %H:%M:%S')


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


if __name__ == '__main__':
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
    parser.add_argument('--output-object-file', '-o', required=True,
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
    args = parser.parse_args()
    logging.debug(args)

    if len(args.subset_names) != len(args.low_thresholds):
        logging.error('--low-thresholds should have a comma separated list of numerics of the same size as --subset-names')
        sys.exit(1)

    if len(args.subset_names) != len(args.high_thresholds):
        logging.error('--high-thresholds should have a comma separated list of numerics of the same size as --subset-names')
        sys.exit(1)

    if args.input_format == 'auto-detect':
        if args.input_object_file.endswith('.loom'):
            args.input_format == 'loom'
        elif args.input_object_file.endswith('.h5ad'):
            args.input_format == 'anndata'
        else:
            logging.error('unknown input format: "{}", please check suffix is either ".loom" or ".h5ad"'.format(args.input_object_file))
            sys.exit(1)

    import scanpy.api as sc

    if args.input_format == 'anndata':
        adata = sc.read(args.input_object_file)
    elif args.input_format == 'loom':
        adata = sc.read_loom(args.input_object_file)
    else:
        logging.error('should not reach here')
        sys.exit(1)

    for name,h,l in zip(args.subset_names, args.high_thresholds, args.low_thresholds):
        adata = adata[adata.obs[name] < h & adata.obs[name] > l, :]

    if args.output_format == 'loom':
        adata.write_loom(args.output_object_file)
    elif args.output_format == 'anndata':
        adata.write(args.output_object_file)
