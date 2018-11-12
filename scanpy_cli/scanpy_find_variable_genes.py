#!/usr/bin/env python3

import logging
import matplotlib
matplotlib.use('Agg')

import scanpy.api as sc
from scanpy_cli.wrapper_utils import (
    ScanpyArgParser,
    read_input_object, write_output_object, save_output_plot,
)


def main(argv=None):
    argparser = ScanpyArgParser(argv, 'Keep variable genes')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_output_plot()
    argparser.add_argument('--flavor',
                           choices=['seurat', 'cell_ranger'],
                           default='seurat',
                           help='Choose the flavor for computing normalised dispersion. '
                                'Default: seurat')
    argparser.add_subset_parameters(params=['mean', 'disp'])
    argparser.add_argument('-b', '--n-bins',
                           type=int,
                           default=20,
                           help='Number of bins for binning the mean gene expression. '
                                'Normalisation is done with respect to each bin. Default: 20')
    argparser.add_argument('-n', '--n-top-genes',
                           type=int,
                           default=None,
                           help='Number of highly variable genes to keep, '
                                'ignoring --subset-parameters when set. Default: None')
    args = argparser.args

    logging.debug(args)

    adata = read_input_object(args.input_object_file, args.input_format)

    min_mean, max_mean, min_disp, max_disp = None, None, None, None
    for name, high, low in zip(args.parameter_names, args.high_thresholds, args.low_thresholds):
        if name == 'mean':
            min_mean = low
            max_mean = high
        elif name == 'disp':
            min_disp = low
            max_disp = high
        else:
            msg = 'Unsupported parameter name "{}", omitted'.format(name)
            logging.warning(msg)

    if args.output_plot:
        sc.set_figure_params(format=args.output_plot_format)
        filtered = sc.pp.filter_genes_dispersion(adata.X,
                                                 flavor=args.flavor,
                                                 min_mean=min_mean, max_mean=max_mean,
                                                 min_disp=min_disp, max_disp=max_disp,
                                                 n_bins=args.n_bins,
                                                 n_top_genes=args.n_top_genes,
                                                 log=(args.flavor == 'seurat'))
        sc.settings.verbosity = 0
        sc.pl.filter_genes_dispersion(filtered, show=False, save=True)
        save_output_plot('filter_genes_dispersion', args.output_plot)

    sc.pp.filter_genes_dispersion(adata,
                                  flavor=args.flavor,
                                  min_mean=min_mean, max_mean=max_mean,
                                  min_disp=min_disp, max_disp=max_disp,
                                  n_bins=args.n_bins,
                                  n_top_genes=args.n_top_genes,
                                  log=(args.flavor == 'seurat'))

    write_output_object(adata, args.output_object_file, args.output_format)

    logging.info('Done')
    return 0


if __name__ == '__main__':
    main()
