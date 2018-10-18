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

    sc.tl.tsne(adata,
               n_pcs=args.n_pcs,
               use_rep=args.use_rep,
               perplexity=args.perplexity,
               early_exaggeration=args.early_exaggeration,
               learning_rate=args.learning_rate,
               random_state=args.random_seed,
               use_fast_tsne=args.use_fast_tsne,
               n_jobs=args.n_jobs)

    write_output_object(adata, args.output_object_file, args.output_format)

    if args.output_embeddings_file:
        adata.obsm.to_df()[['X_tsne1', 'X_tsne2']].to_csv(args.output_embeddings_file, index=None)

    if args.output_plot:
        sc.set_figure_params(format=args.output_plot_format)
        sc.settings.verbosity = 0
        sc.pl.tsne(adata, show=False, save=True,
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

        save_output_plot('tsne', args.output_plot)

    logging.info('Done')
    return 0


if __name__ == '__main__':
    argparser = ScanpyArgParser('Run t-SNE on data with neighborhood graph computed')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_argument('--output-embeddings-file',
                           default=None,
                           help='File name in which to store a csv-format embeddings table with '
                                'coordinates by cell.')
    argparser.add_argument('-n', '--n-pcs',
                           type=int,
                           default=None,
                           help='The number of PCs to use.')
    argparser.add_argument('-r', '--use-rep',
                           default=None,
                           help='Use the indicated representation. If None, the representation '
                                'is chosen automatically: for .n_vars < 50, .X is used, otherwise '
                                'X_pca is used. If X_pca is not present, it\'s computed with '
                                'default parameters. Default: None')
    argparser.add_argument('--perplexity',
                           type=float,
                           default=30,
                           help='The perplexity is related to the number of nearest neighbors '
                                'that is used in other manifold learning algorithms. Larger '
                                'datasets usually require a larger perplexity. Consider selecting '
                                'a value between 5 and 50. The choice is not extremely critical '
                                'since t-SNE is quite insensitive to this parameter. '
                                'Default: 30')
    argparser.add_argument('--early-exaggeration',
                           type=float,
                           default=12.0,
                           help='Controls how tight natural clusters in the original space are '
                                'in the embedded space and how much space will be between '
                                'them. For larger values, the space between natural clusters will '
                                'be larger in the embedded space. Again, the choice of this '
                                'parameter is not very critical. Default: 12')
    argparser.add_argument('--learning-rate',
                           type=float,
                           default=1000,
                           help='The learning rate can be a critical parameter. It should be '
                                'between 100 and 1000. If the cost function increases during '
                                'initial optimization, the early exaggeration factor or the '
                                'learning rate might be too high. If the cost function gets stuck '
                                'in a bad local minimum increasing the learning rate helps '
                                'sometimes. Default: 1000')
    argparser.add_argument('--no-fast-tsne',
                           dest='use_fast_tsne',
                           default=True,
                           action='store_false',
                           help='The multicoreTSNE package is used by default. Set this flag to '
                                'not using it.')
    argparser.add_argument('--n-jobs',
                           type=int,
                           default=None,
                           help='Number of jobs. If None, use sc.settings.n_jobs')
    argparser.add_argument('-s', '--random-seed',
                           type=int,
                           default=0,
                           help='The seed used by the random number generator. Default: 0')
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
