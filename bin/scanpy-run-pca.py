#!/usr/bin/env python

from __future__ import print_function
import logging
import matplotlib
matplotlib.use('Agg')
from scanpy_wrapper_utils import ScanpyArgParser, comma_separated_list
from scanpy_wrapper_utils import read_input_object, write_output_object, save_output_plot


def main(args):
    logging.debug(args)
    import scanpy.api as sc

    adata = read_input_object(args.input_object_file, args.input_format)

    sc.tl.pca(adata,
              n_comps=args.n_pcs,
              zero_center=args.zero_center,
              svd_solver=args.svd_solver,
              random_state=args.random_seed,
              chunked=args.chunked,
              chunk_size=args.chunk_size)

    write_output_object(adata, args.output_object_file, args.output_format)

    if args.output_plot:
        sc.set_figure_params(format=args.output_plot_format)
        sc.settings.verbosity = 0
        sc.pl.pca(adata, show=False, save=True,
                  color=args.color,
                  use_raw=args.use_raw,
                  edges=args.edges,
                  arrows=args.arrows,
                  sort_order=args.sort_order,
                  groups=args.groups,
                  projection=args.projection,
                  components=args.components,
                  palette=args.palette,
                  frameon=args.frameon)
        save_output_plot('pca', args.output_plot)

    logging.info('Done')
    return 0


if __name__ == '__main__':
    argparser = ScanpyArgParser('Run PCA on normalised data')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_argument('-n', '--n-pcs',
                           type=int,
                           default=50,
                           help='Number of principal components to compute. Default: 50')
    argparser.add_argument('-z', '--zero-center',
                           action='store_true',
                           dest='zero_center',
                           default=None,
                           help='Compute standard PCA from covariance matrix. '
                                'See also --no-zero-center. '
                                'If neither --zero-center nor --no-zero-center is set, '
                                'automatically set to TRUE if input is sparse, otherwise FALSE.')
    argparser.add_argument('-Z', '--no-zero-center',
                           action='store_false',
                           dest='zero_center',
                           help='Compute PCA without zero-centering using TruncatedSVD '
                                'from scikit-learn. See also --zero-center.')
    argparser.add_argument('--svd-solver',
                           choices=['arpack', 'randomised', 'auto'],
                           default='auto',
                           help='"arpack" for scipy ARPACK wrapper, or '
                                '"randomised" for randomised algorithm by Halko (2009). '
                                '"auto" chooces automatically depending on problem size. '
                                'Default: auto')
    argparser.add_argument('-s', '--random-seed',
                           type=int,
                           default=0,
                           help='Random seed for initialising optimisation')
    argparser.add_argument('-c', '--chunked',
                           action='store_true',
                           help='Perform an incremental PCA on segments of --chunk-size. '
                                'Imply --zero-center, ignore --svd-solver and --random-seed.')
    argparser.add_argument('--chunk-size',
                           type=int,
                           default=None,
                           help='Number of observations to include in each chunk. '
                                'Required if --chunked is set')
    argparser.add_output_plot()
    argparser.add_argument('--color',
                           type=comma_separated_list('color', str),
                           default=[],
                           help='String or list of strings. Default: []')
    argparser.add_argument('--use-raw',
                           action='store_true',
                           default=False,
                           help='Use raw attribute of adata if present. Default: False')
    argparser.add_argument('--edges',
                           action='store_true',
                           default=False,
                           help='Show edges. Default: False.')
    argparser.add_argument('--arrows',
                           action='store_true',
                           default=False,
                           help='Show arrwos (requires to run rna_velocity() before). '
                                'Default: False.')
    argparser.add_argument('--no-sort-order',
                           dest='sort_order',
                           action='store_false',
                           default=True,
                           help='For continuous annotations used as color parameter, by default '
                                'plot data points with higher values on top of others. Disable '
                                'this behavior if set.')
    argparser.add_argument('--groups',
                           type=str,
                           default=None,
                           help='Restrict to a few categories in observation annotation.')
    argparser.add_argument('--projection',
                           choices=['2d', '3d'],
                           default='2d',
                           help='Projection of plot. Default: 2d')
    argparser.add_argument('--components',
                           type=str,
                           default='1,2',
                           help='Components to plot. To plot all available components use "all". '
                                'Default: "1,2"')
    argparser.add_argument('--palette',
                           default=None,
                           help='Colors to use for plotting categorical annotation groups. '
                                'Can be a valid matplotlib.pyplot.colormap name. Default: None')
    argparser.add_argument('--frameoff',
                           dest='frameon',
                           action='store_false',
                           default=True,
                           help='Do not draw a frame around the scatter plot. Draw by default.')
    args = argparser.get_args()

    main(args)
