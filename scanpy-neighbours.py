#!/usr/bin/env python

from __future__ import print_function
import logging
from scanpy_wrapper_utils import ScanpyArgParser
from scanpy_wrapper_utils import read_input_object, write_output_object


def main(args):
    logging.debug(args)
    import scanpy.api as sc

    adata = read_input_object(args.input_object_file, args.input_format)

    sc.pp.neighbors(adata,
                    n_neighbors=args.n_neighbors,
                    n_pcs=args.n_pcs,
                    use_rep=args.use_rep,
                    knn=args.knn,
                    random_state=args.random_seed,
                    method=args.method,
                    metric=args.metric)

    write_output_object(adata, args.output_object_file, args.output_format)

    print(adata)
    logging.info('Done')
    return 0


if __name__ == '__main__':
    argparser = ScanpyArgParser('Compute neighborhood graph on PCA analysed data')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_argument('-N', '--n-neighbors',
                           type=int,
                           default=15,
                           help='Size of local neighbourhood used for manifold approximation. '
                                'Default: 15')
    argparser.add_argument('-n', '--n-pcs',
                           type=int,
                           default=None,
                           help='Number of principal components to use. Default: None')
    argparser.add_argument('-r', '--use-rep',
                           default=None,
                           help='Use the indicated representation. '
                                'If None, the representation is chosen automatically: '
                                'for .n_vars < 50, .X is used, otherwise "X_pca" is used.')
    argparser.add_argument('--knn',
                           dest='knn',
                           action='store_true',
                           default=True,
                           help='If True, use a hard threshold to restrict the number of neighbors '
                                'to --n-neigbors. Otherwise, use a Gaussian Kernel to assign low '
                                'weights to neighbors more distant then the --n-neighbors nearest '
                                'neighbors.')
    argparser.add_argument('--no-knn',
                           dest='knn',
                           action='store_false',
                           help='See --knn')
    argparser.add_argument('-s', '--random-seed',
                           type=int,
                           default=0,
                           help='Random seed for numpy')
    argparser.add_argument('-m', '--method',
                           choices=['umap', 'gauss'],
                           default='umap',
                           help='Use "umap" or "gauss" with adpative width for computing '
                                'connectivities. Default: "umap"')
    argparser.add_argument('-M', '--metric',
                           default='euclidean',
                           help='A known metric\'s name. Choices are "euclidean", "l2", "l1", '
                                '"manhattan", "cityblock", "braycurtis", "canberra", "chebyshev", '
                                '"correlation", "cosine", "dice", "hamming", "jaccard", '
                                '"kulsinski", "mahalanobis", "matching", "minkowski", '
                                '"rogerstanimoto", "russellrao", "seuclidean", "sokalmichener", '
                                '"sokalsneath", "sqeuclidean", "yule", "wminkowski", '
                                '"precomputed". Default "euclidean".')
    args = argparser.get_args()

    main(args)
