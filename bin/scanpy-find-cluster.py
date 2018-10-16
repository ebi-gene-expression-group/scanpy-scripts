#!/usr/bin/env python

from __future__ import print_function
import logging
from scanpy_wrapper_utils import ScanpyArgParser, comma_separated_list
from scanpy_wrapper_utils import read_input_object, write_output_object


def main(args):
    logging.debug(args)
    import scanpy.api as sc

    adata = read_input_object(args.input_object_file, args.input_format)

    sc.tl.louvain(adata,
                  flavor=args.flavor,
                  resolution=args.resolution,
                  restrict_to=args.restrict_to,
                  key_added=args.key_added,
                  use_weights=args.use_weights,
                  random_state=args.random_seed)

    write_output_object(adata, args.output_object_file, args.output_format)

    logging.info('Done')
    return 0


if __name__ == '__main__':
    argparser = ScanpyArgParser('Run Louvain clustering on data with neighborhood graph computed')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_argument('--flavor',
                           choices=['vtraag', 'igraph'],
                           default='vtraag',
                           help='Choose between two packages for computing the clustering.'
                                '"vtraag" is much more powerful, and the default.')
    argparser.add_argument('--resolution',
                           type=float,
                           default=1.0,
                           help='For the default flavor "vtraag", you can provide a resolution '
                                '(higher resolution means finding more and smaller clusters). '
                                'Default: 1.0')
    argparser.add_argument('--restrict-to',
                           type=comma_separated_list('restrict-to', str),
                           default=[],
                           help='Restrict the clustering to the categories within the key for '
                                'sample annotation, tuple needs to contain (obs key, list of '
                                'categories).')
    argparser.add_argument('--key-added',
                           default='louvain',
                           help='Key under which to add the cluster labels. Default: louvain')
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
