#!/usr/bin/env python

import sys
import logging
import argparse
import scanpy.api as sc


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
        self.events = [{'handler':'_set_logging_level',
                        'argv':'debug'}]


    def add_input_object(self):
        self.parser.add_argument('-i', '--input-object-file',
                                 required=True,
                                 help='Path to anndata or loom file.')
        self.parser.add_argument('-f', '--input-format',
                                 choices=['loom', 'anndata', 'auto-detect'],
                                 default='auto-detect',
                                 help='Format for input object: loom/anndata/[auto-detect].')
        self.events.append({'handler':'_detect_io_format',
                            'argv':['input_format','input_object_file']})


    def add_output_object(self):
        self.parser.add_argument('-o', '--output-object-file',
                                 required=True,
                                 help='File name in which to store serialized python object.')
        self.parser.add_argument('-F', '--output-format',
                                 choices=['loom', 'anndata', 'auto-detect'],
                                 default='auto-detect',
                                 help='Format for output object: loom/anndata/[auto-detect].')
        self.events.append({'handler':'_detect_io_format',
                            'argv':['output_format','output_object_file']})


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
        self.events.append({'handler':'_check_parameter_range',
                            'argv':['subset_names','low_thresholds','high_thresholds']})


    def _set_logging_level(self, argv):
        debug = getattr(self.args, argv)
        log_level = logging.DEBUG if debug else logging.WARN
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s; %(levelname)s; %(filename)s; %(funcName)s(): %(message)s',
            datefmt='%y-%m-%d %H:%M:%S')


    def _detect_io_format(self, argv):
        fmt_key,fn_key = argv
        fmt = getattr(self.args, fmt_key)
        fn = getattr(self.args, fn_key)
        if fmt != 'auto-detect':
            return
        if fn.endswith('.loom'):
            setattr(self.args, fmt_key, 'loom')
        elif fn.endswith('.h5ad'):
            setattr(self.args, fmt_key, 'anndata')
        else:
            logging.error('Unspecified unknown format: "{}", '
                          'please check suffix is either ".loom" or ".h5ad"'.format(fn))
            sys.exit(1)


    def _check_parameter_range(self, argv):
        names,lows,highs = [getattr(self.args,k) for k in argv]
        n = len(names)
        if len(lows) == 0:
            lows = [float('-Inf')] * n
        elif len(lows) != n:
            logging.error('--low-thresholds should be a comma separated list of numerics of the same size as {}'.format(argv[0]))
            sys.exit(1)
        if len(highs) == 0:
            highs = [float('Inf')] * n
        elif len(highs) != n:
            logging.error('--high-thresholds should be a comma separated list of numerics of the same size as {}'.format(argv[0]))
            sys.exit(1)


    def _handle_event(self, ev):
        handler = getattr(self, ev['handler'])
        handler(ev['argv'])


    def get_args(self):
        self.args = self.parser.parse_args()

        for ev in self.events:
            self._handle_event(ev)

        return self.args


def read_input_object(filename, fmt):
    if fmt == 'anndata':
        adata = sc.read(filename)
    elif fmt == 'loom':
        adata = sc.read_loom(filename)
    else:
        logging.error('should not reach here')
        sys.exit(1)
    return adata

def write_output_object(adata, filename, fmt):
    if fmt == "anndata":
        adata.write(filename)
    elif fmt == "loom":
        adata.write_loom(filename)
    else:
        logging.error('should not reach here')
        sys.exit(1)
