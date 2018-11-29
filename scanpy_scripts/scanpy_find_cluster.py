#!/usr/bin/env python

import logging

from scanpy_scripts.wrapper_utils import (
    ScanpyArgParser, comma_separated_list,
    read_input_object, write_output_object,
)


def main(argv=None):
    argparser = ScanpyArgParser(argv, 'Run Louvain clustering on data with neighborhood graph computed')
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
                           type=float,
                           default=1.0,
                           help='For the default flavor "vtraag", you can provide a resolution '
                                '(higher resolution means finding more and smaller clusters). '
                                'Default: 1.0')
    argparser.add_argument('--restrict-to',
                           type=comma_separated_list('restrict-to', str),
                           default=None,
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
    args = argparser.args

    logging.debug(args)
    import scanpy.api as sc

    adata = read_input_object(args.input_object_file, args.input_format)

    if args.restrict_to is not None:
        args.restrict_to = (args.restrict_to[0], args.restrict_to[1:])
    sc.tl.louvain(adata,
                  flavor=args.flavor,
                  resolution=args.resolution,
                  restrict_to=args.restrict_to,
                  key_added=args.key_added,
                  use_weights=args.use_weights,
                  random_state=args.random_seed)

    write_output_object(adata, args.output_object_file, args.output_format)

    if args.output_text_file:
        adata.obs[[args.key_added]].reset_index(level=0).rename(columns={'index':'cells'}).to_csv(
            args.output_text_file, sep='\t', index=None)

    logging.info('Done')
    return 0


if __name__ == '__main__':
    main()
