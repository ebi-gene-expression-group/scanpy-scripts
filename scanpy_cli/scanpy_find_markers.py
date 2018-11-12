#!/usr/bin/env python3

import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import scanpy.api as sc

from scanpy_cli.wrapper_utils import (
    ScanpyArgParser, comma_separated_list,
    read_input_object, write_output_object, save_output_plot,
)


def main(argv=None):
    argparser = ScanpyArgParser(argv, 'Keep variable genes')
    argparser.add_input_object()
    argparser.add_output_object()
    argparser.add_argument('--output-text-file',
                           default=None,
                           help='File name in which to store a text table of marker genes.')
    argparser.add_argument('-g', '--groupby',
                           default='louvain',
                           help='The key of the observations grouping to consider. '
                                'Default: louvain')
    argparser.add_argument('--no-raw',
                           dest='use_raw',
                           action='store_false',
                           default=True,
                           help='Use raw attribute if present.')
    argparser.add_argument('--groups',
                           type=comma_separated_list('groups', str),
                           default='all',
                           help='Subset of groups, e.g. "g1,g2,g3", to which comparison shall be '
                                'restricted. If not passed, a ranking will be generated for all '
                                'groups.')
    argparser.add_argument('--reference',
                           default='rest',
                           help='If "rest", compare each group to the union of the rest of the '
                                'group. If a group identifier, compare with respect to this group.')
    argparser.add_argument('-n', '--n-genes',
                           type=int,
                           default=100,
                           help='The number of genes that appear in the returned tables.')
    argparser.add_argument('-m', '--method',
                           choices=['logreg', 't-test', 'wilcoxon', 't-test_overestim_var'],
                           default='t-test_overestim_var',
                           help='If "t-test", uses t-test, if "wilcoxon", uses '
                                'Wilcoxon-Rank-Sum. If "t-test_overestim_var", overestimates '
                                'variance of each group. If "logreg" uses logistic regression. '
                                'Default: "t-test_overestim_var".')
    argparser.add_argument('--rankby_abs',
                           action='store_true',
                           default=False,
                           help='Rank genes by the absolute value of the score, not by the '
                                'score. The returned scores are never the absolute values.')
    argparser.add_output_plot()
    argparser.add_argument('--show-n-genes',
                           type=int,
                           default=10,
                           help='Number of genes to show per group in plot. Default: 10')
    argparser.add_argument('--key',
                           default=None,
                           help='Key used to store the ranking results in adata.uns.')
    args = argparser.args

    logging.debug(args)

    adata = read_input_object(args.input_object_file, args.input_format)

    if len(args.groups) == 1:
        args.groups = args.groups[0]

    sc.tl.rank_genes_groups(adata,
                            groupby=args.groupby,
                            use_raw=args.use_raw,
                            groups=args.groups,
                            reference=args.reference,
                            n_genes=args.n_genes,
                            method=args.method,
                            rankby_abs=args.rankby_abs)

    write_output_object(adata, args.output_object_file, args.output_format)

    if args.output_text_file:
        export_marker_table(adata, args.output_text_file)

    if args.output_plot:
        sc.set_figure_params(format=args.output_plot_format)
        sc.settings.verbosity = 0
        sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=args.show_n_genes, key=args.key,
                                               show=False, save=True)
        save_output_plot('stacked_violin', args.output_plot)

    logging.info('Done')
    return 0


def recarray_to_array(rec):
    dtyp = rec.dtype[0]
    return rec.view(dtyp)


def export_marker_table(adata, filename):
    result = adata.uns['rank_genes_groups']
    n_genes = len(result['names'])
    n_group = len(result['names'][0])
    groups = np.tile(np.arange(n_group), n_genes)
    names = recarray_to_array(result['names'])
    scores = recarray_to_array(result['scores'])
    # FIXME wilcoxon test doesn't produce logfoldchanges?
    logfc = recarray_to_array(result['logfoldchanges'])
    pvals = recarray_to_array(result['pvals'])
    pvals_adj = recarray_to_array(result['pvals_adj'])
    marker_tbl = pd.DataFrame({'names':names, 'groups':groups, 'scores':scores, 'logfc':logfc,
                               'pvals':pvals, 'pvals_adj':pvals_adj})
    marker_tbl.sort_values(by=['groups', 'pvals_adj']).to_csv(filename, index=None)


if __name__ == '__main__':
    main()
