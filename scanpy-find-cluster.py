#!/usr/bin/env python

from __future__ import print_function
import logging
from scanpy_wrapper_utils import ScanpyArgParser, comma_separated_list
from scanpy_wrapper_utils import read_input_object, write_output_object


def main(args):
    logging.debug(args)
    import scanpy.api as sc

    adata = read_input_object(args.input_object_file, args.input_format, sparse=False)

    if args.restrict_to is not None:
        args.restrict_to = (args.restrict_to[0], args.restrict_to[1:])
    cluster_columns = []
    for res in args.resolution:
        if len(args.resolution) == 1:
            key = args.key_added
        else:
            key = '{}_r{}'.format(args.key_added, res)
        sc.tl.louvain(adata,
                      flavor=args.flavor,
                      resolution=res,
                      restrict_to=args.restrict_to,
                      key_added=key,
                      use_weights=args.use_weights,
                      random_state=args.random_seed)
        cluster_columns.append(key)

    write_output_object(adata, args.output_object_file, args.output_format)

    if args.output_text_file:
        adata.obs[cluster_columns].reset_index(level=0).rename(columns={'index':'cells'}).to_csv(
            args.output_text_file, sep='\t', index=None)

    print(adata)
    logging.info('Done')
    return 0


if __name__ == '__main__':
    argparser = ScanpyArgParser('Run Louvain clustering on data with neighborhood graph computed')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_argument('--output-text-file',
                           default=None,
                           help='File name in which to store text format set of clusters')
    argparser.add_argument('--flavor',
                           choices=['vtraag', 'igraph'],
                           default='vtraag',
                           help='Choose between two packages for computing the clustering.'
                                '"vtraag" is much more powerful, and the default.')
    argparser.add_argument('--resolution',
                           type=comma_separated_list('resolution', float),
                           default=[1.0],
                           help='For the default flavor "vtraag", you can provide a resolution '
                                '(higher resolution means finding more and smaller clusters). '
                                'Multiple resolutions can be provided as a comma-separated-list, '
                                'e.g. "0.5,1.0,1.5". Default: 1.0')
    argparser.add_argument('--restrict-to',
                           type=comma_separated_list('restrict-to', str),
                           default=None,
                           help='Restrict the clustering to the categories within the key for '
                                'sample annotation, tuple needs to contain (obs key, list of '
                                'categories).')
    argparser.add_argument('--key-added',
                           default='louvain',
                           help='Key under which to add the cluster labels. If more than one '
                           'resolution is provided, the key is used as prefix where the full key '
                           'will be "{prefix}_r{resolution}". Default: louvain')
    argparser.add_argument('--use-weights',
                           action='store_true',
                           default=False,
                           help='Use weights from knn graph.')
    argparser.add_argument('-s', '--random-seed',
                           type=int,
                           default=0,
                           help='The seed used to initialise optimisation. Default: 0')
    args = argparser.get_args()

    main(args)
