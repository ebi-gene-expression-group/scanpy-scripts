#!/usr/bin/env python
"""scanpy_wrapper_utils

* Provides
    1. A utilities class, ScanpyArgParser, for command line parsing.
    2. Utilities functions, read_input_object() and write_output_object(), for IO operations.
"""

import sys
import logging
import argparse
import scanpy.api as sc


def comma_separated_list(arg_name, data_type):
    """Generate functions that parse command arguments in the form of
    comma separated string.

    * Parameters
        + arg_name : str
        Name of the command line argument without the leading '--'.
        + data_type : type or function
        Python data type or function that convert a string to a value of specified type.

    * Returns:
        + f : function
        A function that converts a comma separated string to a list of values of specified type.
    """
    def parse_comma_separated_list(arg, name=arg_name, dtype=data_type):
        strings = arg.split(',')
        values = []
        for string in strings:
            try:
                value = dtype(string)
            except ValueError:
                msg = '--{} expects comma separated list of [{}] but received "{}"'.format(
                    name, dtype.__name__, string)
                logging.error(msg)
                sys.exit(1)
            values.append(value)
        return values

    return parse_comma_separated_list


class ScanpyArgParser():
    """An command line argument parser for scanpy wrapper scripts using argparse.ArgumentParser

    * Parameters
        + description : str
        A string passed to argparse.ArgumentParser() summarising the usage of the script.

    * Attributes
        + parser : argparse.ArgumentParser
        The actual argument parser that parses the command line input
    """
    def __init__(self, description=None):
        self.parser = argparse.ArgumentParser(description=description)
        self.parser.add_argument('--debug', action='store_true',
                                 help='Print debug information.')
        self._events = [{'handler':self._set_logging_level,
                         'argv':[]}]
        self._args = None


    def add_input_object(self):
        """Add options "-i/--input-object-file" and "-f/--input-format""""
        self.parser.add_argument('-i', '--input-object-file',
                                 required=True,
                                 help='Path to anndata or loom file.')
        self.parser.add_argument('-f', '--input-format',
                                 choices=['loom', 'anndata', 'auto-detect'],
                                 default='auto-detect',
                                 help='Format for input object: loom/anndata/[auto-detect].')
        self._events.append({'handler':self._detect_io_format,
                             'argv':['input_format', 'input_object_file']})


    def add_output_object(self):
        """Add options "-o/--output-object-file" and "-F/--output-format""""
        self.parser.add_argument('-o', '--output-object-file',
                                 required=True,
                                 help='File name in which to store serialized python object.')
        self.parser.add_argument('-F', '--output-format',
                                 choices=['loom', 'anndata', 'auto-detect'],
                                 default='auto-detect',
                                 help='Format for output object: loom/anndata/[auto-detect].')
        self._events.append({'handler':self._detect_io_format,
                             'argv':['output_format', 'output_object_file']})


    def add_subset_parameters(self):
        """Add options "-s/--subset-names", "-l/--low-thresholds" and "-j/--high-thresholds""""
        self.parser.add_argument('-s', '--subset-names',
                                 required=True,
                                 type=comma_separated_list('subset-names', str),
                                 default=[],
                                 help='Parameters to subset on.')
        self.parser.add_argument('-l', '--low-thresholds',
                                 type=comma_separated_list('low-thresholds', float),
                                 default=[],
                                 help='Low cutoffs for the parameters (default is -Inf).')
        self.parser.add_argument('-j', '--high-thresholds',
                                 type=comma_separated_list('high-thresholds', float),
                                 default=[],
                                 help='High cutoffs for the parameters (default is Inf).')
        self._events.append({'handler':self._check_parameter_range,
                             'argv':['subset_names', 'low_thresholds', 'high_thresholds']})


    def _set_logging_level(self):
        log_level = logging.DEBUG if self._args.debug else logging.WARN
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s; %(levelname)s; %(filename)s; %(funcName)s(): %(message)s',
            datefmt='%y-%m-%d %H:%M:%S')


    def _detect_io_format(self, fmt_key, fname_key):
        fmt = getattr(self._args, fmt_key)
        fname = getattr(self._args, fname_key)
        if fmt != 'auto-detect':
            return
        if fname.endswith('.loom'):
            setattr(self._args, fmt_key, 'loom')
        elif fname.endswith('.h5ad'):
            setattr(self._args, fmt_key, 'anndata')
        else:
            msg = ('Unspecified unknown format: "{}", '
                   'please check suffix is either ".loom" or ".h5ad"').format(fname)
            logging.error(msg)
            sys.exit(1)


    def _check_parameter_range(self, name_key, low_key, high_key):
        name_arg, low_arg, high_arg = [k.replace('_', '-') for k in (name_key, low_key, high_key)]
        names, lows, highs = [getattr(self._args, k) for k in (name_key, low_key, high_key)]
        n_par = len(names)
        if not lows:
            lows = [float('-Inf')] * n_par
        elif len(lows) != n_par:
            msg = ('--{} should be a comma separated list of the same size as --{}').format(
                low_arg, name_arg)
            logging.error(msg)
            sys.exit(1)
        if not highs:
            highs = [float('Inf')] * n_par
        elif len(highs) != n_par:
            msg = ('--{} should be a comma separated list of the same size as --{}').format(
                high_arg, name_arg)
            logging.error(msg)
            sys.exit(1)


    def get_args(self):
        """Return parsed command arguments"""
        if self._args is None:
            self._args = self.parser.parse_args()

            for evt in self._events:
                handler = evt['handler']
                argv = evt['argv']
                handler(*argv)

        return self._args


def read_input_object(filename, fmt):
    """Read an AnnData object from an input file

    * Parameters
        + filename : str
        Path of the input file
        + fmt : str
        Format of the input file, either "loom" or "anndata"

    * Returns
        + adata : AnnData
        An AnnData object
    """
    if fmt == 'anndata':
        adata = sc.read(filename)
    elif fmt == 'loom':
        adata = sc.read_loom(filename)
    else:
        logging.error('should not reach here')
        sys.exit(1)
    return adata

def write_output_object(adata, filename, fmt):
    """Write an AnnData object to an output ile

    * Parameters
        + adata : AnnData
        An AnnData object
        + filename : str
        Path of the output file
        + fmt : str
        Format of the output file, either "loom" or "anndata"
    """
    if fmt == "anndata":
        adata.write(filename)
    elif fmt == "loom":
        adata.write_loom(filename)
    else:
        logging.error('should not reach here')
        sys.exit(1)
