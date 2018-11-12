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
    argparser = ScanpyArgParser(argv, 'Run t-SNE on data with neighborhood graph computed')
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
    argparser.add_scatter_plot_options()
    args = argparser.args

    logging.debug(args)

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

        save_output_plot('tsne', args.output_plot)

    logging.info('Done')
    return 0


if __name__ == '__main__':
    main()
