#!/usr/bin/env python

import sys
import logging
import argparse
import scanpy as sc


def comma_separated_list(name, dtyp):
    supported_dtypes = {'integer':int, 'numeric':float, 'string':str}
    if dtyp not in supported_dtypes:
        logging.warn('unsupported data type: "{}", default to "string"'.format(dtyp))
        dtype_name = 'string'
        dtype = str
    else:
        dtype_name = dtyp
        dtype = supported_dtypes[dtyp]

    def parse_comma_separated_list(arg, arg_name=name, dtype=dtype, dtype_name=dtype_name):
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


class ScanpyArgParser(object):
    def __init__(self, description=None):
        self.parser = argparse.ArgumentParser(description=description)
        self.parser.add_argument('--debug', action='store_true',
                                 help='Print debug information.')
        self.names = set()


    def add_input_object(self):
        self.parser.add_argument('-i', '--input-object-file',
                                 required=True,
                                 help='Path to anndata or loom file.')
        self.parser.add_argument('-f', '--input-format',
                                 choices=['loom', 'anndata', 'auto-detect'],
                                 default='auto-detect',
                                 help='Format for input object: loom/anndata/[auto-detect].')
        self.names.add('input-object-file')
        self.names.add('input-format')


    def add_output_object(self):
        self.parser.add_argument('-o', '--output-object-file',
                                 required=True,
                                 help='File name in which to store serialized python object.')
        self.parser.add_argument('-F', '--output-format',
                                 choices=['loom', 'anndata', 'auto-detect'],
                                 default='auto-detect',
                                 help='Format for output object: loom/anndata/[auto-detect].')
        self.names.add('output-object-file')
        self.names.add('output-format')


    def add_subset_parameters(self):
        self.parser.add_argument('-s', '--subset-names',
                                 required=True,
                                 type=comma_separated_list('subset-names', 'string'),
                                 default=[],
                                 help='Parameters to subset on.')
        self.parser.add_argument('-l', '--low-thresholds',
                                 type=comma_separated_list('low-thresholds', 'numeric'),
                                 default=[],
                                 help='Low cutoffs for the parameters (default is -Inf).')
        self.parser.add_argument('-j', '--high-thresholds',
                                 type=comma_separated_list('high-thresholds', 'numeric'),
                                 default=[],
                                 help='High cutoffs for the parameters (default is Inf).')


    def _set_logging_level(self, debug=False):
        log_level = logging.DEBUG if debug else logging.WARN
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s; %(levelname)s; %(filename)s; %(funcName)s(): %(message)s',
            datefmt='%y-%m-%d %H:%M:%S')


    def _auto_detect_format(self, filename):
        if filename.endswith('.loom'):
            return 'loom'
        elif filename.endswith('.h5ad'):
            return 'anndata'
        else:
            logging.error('unknown input format: "{}", ',
                          'please check suffix is either ".loom" or ".h5ad"'.format(filename))
            sys.exit(1)


    def get_args(self):
        args = self.parser.parse_args()

        self._set_logging_level(args.debug)

        if 'input-format' in self.names and args.input_format == 'auto-detect':
            args.input_format = _auto_detect_format(args.input_object_file)

        if 'output-format' in self.names and args.output_format == 'auto-detect':
            args.output_format = _auto_detect_format(args.output_object_file)

        return args


def read_input_object(filename, format):
    if format == 'anndata':
        adata = sc.read(filename)
    elif format == 'loom':
        adata = sc.read_loom(filename)
    else:
        logging.error('should not reach here')
        sys.exit(1)
    return adata

def write_output_object(adata, filename, format):
    if format == "anndata":
        adata.write(filename)
    elif format == "loom":
        adata.write_loom(filename)
    else:
        logging.error('should not reach here')
        sys.exit(1)
