#!/usr/bin/env python

from __future__ import print_function
import logging
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

    if args.output_plot:
        sc.set_figure_params(format=args.output_plot_format)
        sc.settings.verbosity = 0
        sc.pl.umap(adata, show=False, save=True,
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

        save_output_plot('umap', args.output_plot)

    logging.info('Done')
    return 0


if __name__ == '__main__':
    argparser = ScanpyArgParser('Run PCA on normalised data')
    argparser.add_input_object()
    argparser.add_output_object()
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
