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

    sc.tl.umap(adata,
               min_dist=args.min_dist,
               spread=args.spread,
               n_components=args.n_components,
               maxiter=args.maxiter,
               alpha=args.alpha,
               gamma=args.gamma,
               negative_sample_rate=args.negative_sample_rate,
               init_pos=args.init_pos,
               random_state=args.random_seed,
               a=args.a,
               b=args.b)

    write_output_object(adata, args.output_object_file, args.output_format)

    if args.output_embeddings_file:
        adata.obsm.to_df()[['X_umap1', 'X_umap2']].to_csv(args.output_embeddings_file, index=None)

    if args.output_plot:
        sc.set_figure_params(format=args.output_plot_format)
        sc.settings.verbosity = 0
        sc.pl.umap(adata, show=False, save=True,
                   color=args.color_by,
                   use_raw=args.use_raw,
                   edges=args.edges,
                   arrows=args.arrows,
                   sort_order=args.sort_order,
                   groups=args.groups,
                   projection=args.projection,
                   components=args.components,
                   palette=args.palette,
                   frameon=args.frameon)

        save_output_plot('umap', args.output_plot)

    logging.info('Done')
    return 0


if __name__ == '__main__':
    argparser = ScanpyArgParser('Run UMAP on data with neighborhood graph computed')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_argument('--output-embeddings-file',
                           default=None,
                           help='File name in which to store a csv-format UMAP embeddings table '
                                'with coordinates by cell.')
    argparser.add_argument('--min-dist',
                           type=float,
                           default=0.5,
                           help='The effective minimum distance between embedded points. '
                                'Default: 0.5')
    argparser.add_argument('--spread',
                           type=float,
                           default=1.0,
                           help='The effective scale of embedded points. Default: 1.0')
    argparser.add_argument('-n', '--n-components',
                           type=int,
                           default=2,
                           help='The number of dimensions of the embedding. Default: 2')
    argparser.add_argument('--maxiter',
                           type=int,
                           default=None,
                           help='The number of iterations of the optimisation. Default: None')
    argparser.add_argument('--alpha',
                           type=float,
                           default=1.0,
                           help='The initial learning rate for the embedding optimisation. '
                                'Default: 1.0')
    argparser.add_argument('--gamma',
                           type=float,
                           default=1.0,
                           help='Weighting applied to negative samples in low dimensional '
                                'embedding optimisation. Default: 1.0')
    argparser.add_argument('--negative-sample-rate',
                           type=int,
                           default=5,
                           help='The number of negative edge/1-simplex samples to use per '
                                'positive edge/1-simplex sample in optimising the low dimensional '
                                'embedding. Default: 5')
    argparser.add_argument('--init-pos',
                           default='spectral',
                           help='How to initialise the low dimensional embedding. Choices are: '
                                'any key for adata.obsm, "paga", "spectral", or "random". '
                                'Default: spectral')
    argparser.add_argument('-s', '--random-seed',
                           type=int,
                           default=0,
                           help='The seed used by the random number generator. Default: 0')
    argparser.add_argument('-a',
                           type=float,
                           default=None,
                           help='More specific parameters controlling the embedding. If None, '
                                'these values are set automatically as determined by --min-dist '
                                'and --spread. Default: None')
    argparser.add_argument('-b',
                           type=float,
                           default=None,
                           help='More specific parameters controlling the embedding. If None, '
                                'these values are set automatically as determined by --min-dist '
                                'and --spread. Default: None')
    argparser.add_output_plot()
    argparser.add_scatter_plot_options()
    args = argparser.get_args()

    main(args)
