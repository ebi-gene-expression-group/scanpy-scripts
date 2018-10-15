#!/usr/bin/env python

from __future__ import print_function
import logging
from scanpy_wrapper_utils import ScanpyArgParser
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
        sc.pl.pca(adata, show=False, save=True)
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
                                'automatically set to FALSE if input is sparse, otherwise TRUE.')
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
    args = argparser.get_args()

    main(args)
