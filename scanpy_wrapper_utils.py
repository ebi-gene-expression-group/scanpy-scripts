#!/usr/bin/env python
"""scanpy_wrapper_utils

* Provides
    1. A utilities class, ScanpyArgParser, for command line parsing.
    2. Utilities functions, read_input_object() and write_output_object(), for IO operations.
"""

import sys
import os.path
import shutil
import logging
import argparse
import pandas as pd
import scanpy.api as sc


def comma_separated_list(arg_name, data_type):
    """Generate functions that parse command arguments in the form of comma separated string.

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

    * Methods
        + get_args()
        + add_argument(*args, **kwargs)
        + add_input_object()
        + add_output_object()
        + add_output_plot()
        + add_subset_parameters(params=None)
        + add_subset_list(dtype=str)
    """
    def __init__(self, description=None):
        self._parser = argparse.ArgumentParser(description=description)
        self._parser.add_argument('--debug', action='store_true',
                                  help='Print debug information.')
        self._events = [[self._set_logging_level, [], {}]]
        self._args = None


    def get_args(self):
        """Return parsed command arguments"""
        if self._args is None:
            self._args = self._parser.parse_args()

            for evt in self._events:
                handler, args, kwargs = evt
                handler(*args, **kwargs)

        return self._args


    def add_argument(self, *args, **kwargs):
        """Calls argparse.ArgumentParser.add_argument"""
        self._parser.add_argument(*args, **kwargs)


    def add_input_object(self):
        """Add options -i/--input-object-file and -f/--input-format"""
        self._parser.add_argument('-i', '--input-object-file',
                                  required=True,
                                  help='Path to anndata or loom file.')
        self._parser.add_argument('-f', '--input-format',
                                  choices=['loom', 'anndata', 'auto-detect'],
                                  default='auto-detect',
                                  help='Format for input object: loom/anndata/[auto-detect].')
        self._events.append([self._detect_io_format,
                             ['input_format', 'input_object_file'], {}])


    def add_output_object(self):
        """Add options -o/--output-object-file and -F/--output-format"""
        self._parser.add_argument('-o', '--output-object-file',
                                  required=True,
                                  help='File name in which to store serialized python object.')
        self._parser.add_argument('-F', '--output-format',
                                  choices=['loom', 'anndata', 'auto-detect'],
                                  default='auto-detect',
                                  help='Format for output object: loom/anndata/[auto-detect].')
        self._events.append([self._detect_io_format,
                             ['output_format', 'output_object_file'], {}])


    def add_output_plot(self):
        """Add options -P/--output-plot"""
        self._parser.add_argument('-P', '--output-plot',
                                  default=None,
                                  help='Save plot in the specified file')
        self._events.append([self._check_filename_extension,
                             ['output_plot'],
                             {'extensions': ['pdf', 'png', 'svg']}])


    def add_scatter_plot_options(self):
        """Add options for scatter plots, e.g. PCA/UMAP/tSNE"""
        self._parser.add_argument('--color-by',
                                  type=comma_separated_list('color-by', str),
                                  default=[],
                                  help='String or list of strings. Default: []')
        self._parser.add_argument('--use-raw',
                                  action='store_true',
                                  default=False,
                                  help='Use raw attribute of adata if present. Default: False')
        self._parser.add_argument('--edges',
                                  action='store_true',
                                  default=False,
                                  help='Show edges. Default: False.')
        self._parser.add_argument('--arrows',
                                  action='store_true',
                                  default=False,
                                  help='Show arrwos (requires having run rna_velocity() before). '
                                  'Default: False.')
        self._parser.add_argument('--no-sort-order',
                                  dest='sort_order',
                                  action='store_false',
                                  default=True,
                                  help='For continuous annotations used as color parameter, by default '
                                  'plot data points with higher values on top of others. Disable '
                                  'this behavior if set.')
        self._parser.add_argument('--groups',
                                  type=str,
                                  default=None,
                                  help='Restrict to a few categories in observation annotation.')
        self._parser.add_argument('--projection',
                                  choices=['2d', '3d'],
                                  default='2d',
                                  help='Projection of plot. Default: 2d')
        self._parser.add_argument('--components',
                                  type=str,
                                  default='1,2',
                                  help='Components to plot. To plot all available components use "all". '
                                  'Default: "1,2"')
        self._parser.add_argument('--palette',
                                  default=None,
                                  help='Colors to use for plotting categorical annotation groups. '
                                  'Can be a valid matplotlib.pyplot.colormap name. Default: None')
        self._parser.add_argument('--frameoff',
                                  dest='frameon',
                                  action='store_false',
                                  default=True,
                                  help='Do not draw a frame around the scatter plot. Draw by default.')


    def add_subset_parameters(self, params=None):
        """Add options -p/--parameter-names, -l/--low-thresholds and -j/--high-thresholds

        * Parameters
            + params : list
            If given, only parameters in this list are permitted.
        """
        self._parser.add_argument('-p', '--parameter-names',
                                  required=True,
                                  type=comma_separated_list('parameter-names', str),
                                  default=[],
                                  help='Parameters to subset on.' +
                                  (' Choose from {}'.format(params) if params else ''))
        self._parser.add_argument('-l', '--low-thresholds',
                                  type=comma_separated_list('low-thresholds', float),
                                  default=[],
                                  help='Low cutoffs for the parameters (default is -Inf).')
        self._parser.add_argument('-j', '--high-thresholds',
                                  type=comma_separated_list('high-thresholds', float),
                                  default=[],
                                  help='High cutoffs for the parameters (default is Inf).')
        self._events.append([self._check_parameter_range,
                             ['parameter_names', 'low_thresholds', 'high_thresholds'],
                             {'params':params}])

    def add_subset_list(self, dtype=str):
        """Add option -s/--subset-list

        * Parameters
            + dtype : type or function
            Python data type or function that convert a string to a value of specified type.
        """
        self._parser.add_argument('-s', '--subset-list',
                                  default='',
                                  help='Comma-separated list of entries to use as a subset. '
                                  'Alternatively, text file with one entry per line.')
        self._events.append([self._list_or_read_file,
                             ['subset_list'],
                             {'dtype':dtype}])


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
                   'please check extension is either ".loom" or ".h5ad"').format(fname)
            logging.error(msg)
            sys.exit(1)


    def _check_filename_extension(self, fname_key, extensions):
        fname_opt = fname_key.replace('_', '-')
        fname = getattr(self._args, fname_key)
        if not fname:
            return
        ext = os.path.splitext(fname)[1]
        if extensions and ext[1:] not in extensions:
            msg = ('Unsupported file format from --{}: {}. '
                   'Choose from {}').format(fname_opt, fname, extensions)
            logging.error(msg)
            sys.exit(1)
        setattr(self._args, fname_key + '_format', ext[1:])


    def _check_parameter_range(self, name_key, low_key, high_key, params=None):
        name_opt, low_opt, high_opt = [k.replace('_', '-') for k in (name_key, low_key, high_key)]
        names, lows, highs = [getattr(self._args, k) for k in (name_key, low_key, high_key)]
        if not names:
            return
        if params:
            for name in names:
                if name not in params:
                    msg = 'Unsupported subset parameter: {}. Choose from {}'.format(name, params)
                    logging.error(msg)
                    sys.exit(1)
        n_par = len(names)
        if not lows:
            lows = [float('-Inf')] * n_par
        elif len(lows) != n_par:
            msg = ('--{} should be a comma separated list of the same size as --{}').format(
                low_opt, name_opt)
            logging.error(msg)
            sys.exit(1)
        if not highs:
            highs = [float('Inf')] * n_par
        elif len(highs) != n_par:
            msg = ('--{} should be a comma separated list of the same size as --{}').format(
                high_opt, name_opt)
            logging.error(msg)
            sys.exit(1)


    def _list_or_read_file(self, key, dtype=str):
        opt = key.replace('_', '-')
        str_list_parser = comma_separated_list(opt, dtype)
        input_string = getattr(self._args, key)
        if input_string:
            if ',' in input_string:
                entries = str_list_parser(input_string)
            elif os.path.isfile(input_string):
                entries = list(pd.read_table(input_string,
                                             header=None).iloc[:, 0].values.astype(dtype))
            else:
                try:
                    entries = [dtype(input_string)]
                except ValueError:
                    msg = '--{} expects comma separated list of [{}] but received "{}"'.format(
                        opt, dtype.__name__, input_string)
                    logging.error(msg)
                    sys.exit(1)
            setattr(self._args, key, entries)


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


def export_mtx(adata, fname_prefix, var=[], obs=[], use_raw=False):
    """Export AnnData object to mtx formt

    * Parameters
        + adata : AnnData
        An AnnData object
        + fname_prefix : str
        Prefix of the exported files, full names will be <fname_prefix>_matrix.mtx,
        <fname_prefix>_genes.tsv, <fname_prefix>_barcodes.tsv
        + var : list
        A list of column names to be exported to gene table
        + obs : list
        A list of column names to be exported to barcode/cell table
    """
    import scipy.sparse as sp
    if use_raw:
        adata = adata.raw
    mat = sp.coo_matrix(adata.X)
    n_obs, n_var = mat.shape
    n_entry = len(mat.data)
    header = '%%MatrixMarket matrix coordinate real general\n%\n{} {} {}\n'.format(
        n_var, n_obs, n_entry)
    df = pd.DataFrame({'col': mat.col + 1, 'row': mat.row + 1, 'data': mat.data})
    mtx_fname = fname_prefix + '_matrix.mtx'
    gene_fname = fname_prefix + '_genes.tsv'
    barcode_fname = fname_prefix + '_barcode.tsv'
    with open(mtx_fname, 'a') as f:
        f.write(header)
        df.to_csv(f, sep=' ', header=False, index=False)

    obs = list(set(obs) & set(adata.obs.columns))
    obs_df = adata.obs[obs].reset_index(level=0)
    obs_df.to_csv(barcode_fname, sep='\t', header=False, index=False)
    var = list(set(var) & set(adata.var.columns))
    var_df = adata.var[var].reset_index(level=0)
    if len(var) == 0:
        var_df['gene'] = var_df['index']
    var_df.to_csv(gene_fname, sep='\t', header=False, index=False)


def save_output_plot(func_name, filename):
    """Save output plot to the specified location

    * Parameters
        + func_name : str
        name of the function that produced the plot
        + filename: str
        path of the saved plot
    """
    autosaved = '{}{}{}.{}'.format(sc.settings.figdir,
                                   func_name,
                                   sc.settings.plot_suffix,
                                   sc.settings.file_format_figs)
    dirname = os.path.split(filename)[0]
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)
    shutil.move(autosaved, filename)
