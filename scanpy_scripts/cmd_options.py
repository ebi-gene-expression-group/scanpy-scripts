"""
Provide cmd options
"""

import click
from .click_utils import (
    CommaSeparatedText,
    Dictionary,
    valid_limit,
    valid_parameter_limits,
    mutually_exclusive_with,
    required_by,
)

COMMON_OPTIONS = {
    'input': [
        click.argument(
            'input_obj',
            metavar='<input_obj>',
            type=click.Path(exists=True, dir_okay=False),
        ),
        click.option(
            '--input-format', '-f',
            type=click.Choice(['anndata', 'loom']),
            default='anndata',
            show_default=True,
            help='Input object format.',
        ),
    ],

    'output': [
        click.argument(
            'output_obj',
            metavar='<output_obj>',
            type=click.Path(dir_okay=False, writable=True),
        ),
        click.option(
            '--output-format', '-F',
            type=click.Choice(['anndata', 'loom', 'zarr']),
            default='anndata',
            show_default=True,
            help='Output object format.',
        ),
        click.option(
            '--zarr-chunk-size', '-z',
            type=click.INT,
            default=1000,
            show_default=True,
            help='Chunk size for writing output in zarr format.',
        ),
        click.option(
            '--export-mtx', '-X',
            type=click.Path(dir_okay=True, writable=True),
            default=None,
            show_default=True,
            help='When specified, using it as prefix for exporting mtx files. '
            'If not empty and not ending with "/" or "_", a "_" will be '
            'appended.',
        ),
        click.option(
            '--show-obj',
            type=click.Choice(['stdout', 'stderr']),
            default=None,
            show_default=True,
            help='Print output object summary info to specified stream.',
        ),
    ],

    'plot': [
        click.argument(
            'output_fig',
            metavar='<output_fig>',
            type=click.Path(dir_okay=False, writable=True),
        ),
        click.option(
            '--fig-size',
            type=CommaSeparatedText(click.INT, length=2),
            default="7,7",
            show_default=True,
            help='Figure size.',
        ),
        click.option(
            '--fig-dpi',
            type=click.INT,
            default=80,
            show_default=True,
            help='Figure DPI.',
        ),
        click.option(
            '--fig-fontsize',
            type=click.INT,
            default=15,
            show_default=True,
            help='Figure font size.',
        ),
    ],

    'use_pc': [
        click.option(
            '--n-pcs', '-n',
            type=click.INT,
            default=None,
            show_default=True,
            help='Use this many PCs. Use `.X` if --n-pcs is 0 when --use-rep is '
            'None.',
        ),

        click.option(
            '--use-rep', '-u',
            type=click.STRING,
            default=None,
            show_default=True,
            help='Use the indicated representation. If None, the representation is '
            'chosen automatically: for `.n_vars` < 50, `.X` is used, otherwise '
            '`X_pca` is used. If `X_pca` is not present, it\'s computed with '
            'default parameters.'
        ),
    ],

    'knn_graph': [
        click.option(
            '--use-graph',
            type=click.STRING,
            default='neighbors',
            show_default=True,
            help='Slot name under `.uns` that contains the KNN graph of which '
            'sparse adjacency matrix is used for clustering.',
        ),
        click.option(
            '--directed/--undirected', 'directed',
            default=True,
            show_default=True,
            help='Interpret the adjacency matrix as directed graph.',
        ),
        click.option(
            '--use-weights',
            is_flag=True,
            default=False,
            show_default=True,
            help='Use weights from KNN graph.',
        ),
    ],

    'n_comps': click.option(
        '--n-comps',
        type=click.INT,
        default=50,
        show_default=True,
        help='Number of components to compute',
    ),

    'key_added': click.option(
        '--key-added',
        type=CommaSeparatedText(simplify=True),
        default=None,
        show_default=True,
        help='Key under which to add the computed results',
    ),

    'random_state': click.option(
        '--random-state', '-S',
        type=click.INT,
        default=0,
        show_default=True,
        help='Seed for random number generator.',
    ),

    'use_raw': click.option(
        '--use-raw/--no-raw', 'use_raw',
        default=True,
        show_default=True,
        help='Use expression values in `.raw` if present.',
    ),

    'init_pos': click.option(
        '--init-pos',
        type=click.STRING,
        default=None,
        help='Use precomputed coordinates for initialization. Can be any key '
        'of `.obsm` or "paga" if .uns["paga"] is present',
    ),

    'zero_center': click.option(
        '--no-zero-center', 'zero_center',
        is_flag=True,
        flag_value=False,
        default=True,
        help='When set, omit zero-centering variables to allow efficient '
        'handling of sparse input.',
    ),

    'n_jobs': click.option(
        '--n-jobs', '-J',
        type=click.INT,
        default=None,
        show_default=True,
        help='Number of jobs for parallel computation.',
    ),

    'restrict_to': click.option(
        '--restrict-to',
        type=(click.STRING, CommaSeparatedText()),
        default=(None, None),
        show_default=True,
        help='Restrict the clustering to the categories within the key for '
        'sample annotation, in the form of "obs_key list_of_categories".',
    ),

    'export_embedding': click.option(
        '--export-embedding', '-E',
        type=click.Path(dir_okay=False, writable=True),
        default=None,
        show_default=True,
        help='Export embeddings in a tab-separated text table.',
    ),
}

READ_CMD_OPTIONS = [
    click.option(
        '--input-10x-h5', '-i',
        type=click.Path(exists=True, dir_okay=False),
        callback=mutually_exclusive_with('--input-10x-mtx'),
        help='Input 10x data in Cell-Ranger hdf5 format.',
    ),
    click.option(
        '--input-10x-mtx', '-x',
        type=click.Path(exists=True, file_okay=False),
        callback=mutually_exclusive_with('--input-10x-h5'),
        help='Path of input folder containing 10x data in mtx format.',
    ),
    *COMMON_OPTIONS['output'],
    click.option(
        '--genome', '-g',
        callback=required_by('--input-10x-h5'),
        default='hg19',
        show_default=True,
        help='Name of the genome group in hdf5 file, required by '
        '"--input-10x-h5".',
    ),
    click.option(
        '--var-names', '-v',
        type=click.Choice(['gene_symbols', 'gene_ids']),
        callback=required_by('--input-10x-mtx'),
        default='gene_symbols',
        show_default=True,
        help='Attribute to be used as the index of the variable table, '
        'required by "--input-10x-mtx".',
    ),
    click.option(
        '--extra-obs',
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        show_default=True,
        help='Extra cell metadata table, must be tab-separated with a header '
        'row and an index column, and with matched dimension.',
    ),
    click.option(
        '--extra-var',
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        show_default=True,
        help='Extra gene metadata table, must be tab-separated with a header '
        'row and an index column, and with matched dimension.',
    ),
]

FILTER_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    click.option(
        '--gene-name', '-g',
        type=click.STRING,
        default='index',
        show_default=True,
        help='Name of the variable that contains gene names, '
        'used for flagging mitochondria genes.',
    ),
    click.option(
        '--list-attr', '-l',
        is_flag=True,
        default=False,
        help='When set, list attributes that can be filtered on.',
    ),
    click.option(
        '--param', '-p',
        type=(click.STRING, click.FLOAT, click.FLOAT),
        multiple=True,
        callback=valid_parameter_limits,
        help='Numerical parameters used to filter the data, '
        'in the format of "-p name min max". '
        'Multiple -p entries allowed.',
    ),
    click.option(
        '--category', '-c',
        type=(click.STRING, CommaSeparatedText()),
        multiple=True,
        help='Categorical attributes used to filter the data, '
        'in the format of "-c <name> <values>", '
        'where entries with attribute <name> with value in <values> are kept. '
        'Multiple -c entries allowed.',
    ),
    click.option(
        '--subset', '-s',
        type=(click.STRING, click.File()),
        multiple=True,
        help='Similar to --category in the format of "-s <name> <file>", '
        'but the <file> to be a one-column table that provides the values. '
        'Multiple -s entries allowed.',
    ),
]

NORM_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    click.option(
        '--save-raw', '-r',
        type=click.Choice(['yes', 'no', 'counts']),
        default='yes',
        show_default=True,
        help='Save raw data existing raw data.',
    ),
    click.option(
        '--normalize-to', '-t', 'target_sum',
        type=float,
        default=10_000,
        show_default=True,
        help='Normalize per cell nUMI to this number.',
    ),
    click.option(
        '--fraction',
        type=float,
        default=0.9,
        show_default=True,
        help='Only use genes that make up less than this fraction of the total '
        'count in every cell. So only these genes will sum up to the number '
        'specified by --normalize-to.',
    ),
]

HVG_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    click.option(
        '--mean-limits', '-m',
        type=(click.FLOAT, click.FLOAT),
        callback=valid_limit,
        default=(0.0125, 3),
        show_default=True,
        help='Cutoffs for the mean of expression'
        'in the format of "-m min max".',
    ),
    click.option(
        '--disp-limits', '-d',
        type=(click.FLOAT, click.FLOAT),
        callback=valid_limit,
        default=(0.5, float('inf')),
        show_default=True,
        help='Cutoffs for the dispersion of expression'
        'in the format of "-d min max".',
    ),
    click.option(
        '--n-bins', '-b',
        type=click.INT,
        default=20,
        show_default=True,
        help='Number of bins for binning the mean gene expression.',
    ),
    click.option(
        '--n-top-genes', '-t',
        type=click.INT,
        default=2000,
        show_default=True,
        help='Number of highly-variable genes to keep.',
    ),
    click.option(
        '--flavor', '-v',
        type=click.Choice(['seurat', 'cellranger']),
        default='seurat',
        show_default=True,
        help='Choose the flavor for computing normalized dispersion.',
    ),
    click.option(
        '--subset', '-s',
        is_flag=True,
        default=False,
        help='When set, inplace subset to highly-variable genes, otherwise '
        'only flag highly-variable genes.',
    ),
    click.option(
        '--by-batch', '-B',
        type=(click.STRING, click.INT),
        default=(None, None),
        show_default=True,
        help='Find highly variable genes within each batch defined by <TEXT> '
        'then pool and keep those found in at least <INTEGER> batches.',
    ),
]

SCALE_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    COMMON_OPTIONS['zero_center'],
    click.option(
        '--max-value', '-m',
        type=click.FLOAT,
        default=None,
        show_default=True,
        help='When specified, clip to this value after scaling, otherwise do '
        'not clip',
    ),
]

REGRESS_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    COMMON_OPTIONS['n_jobs'],
    click.option(
        '--keys', '-k',
        type=CommaSeparatedText(simplify=True),
        default=None,
        show_default=True,
        help='Key(s) for observation annotation on which to regress.',
    ),
]

PCA_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    COMMON_OPTIONS['zero_center'],
    COMMON_OPTIONS['random_state'],
    COMMON_OPTIONS['export_embedding'],
    COMMON_OPTIONS['n_comps'],
    click.option(
        '--svd-solver', '-V',
        type=click.Choice(['auto', 'arpack', 'randomized']),
        default='auto',
        show_default=True,
        help='SVD solver to use.'
    ),
    click.option(
        '--use-all', '-a', 'use_highly_variable',
        is_flag=True,
        flag_value=False,
        default=True,
        help='When set, use all genes for PCA, otherwise use '
        'highly-variable genes by default.'
    ),
    click.option(
        '--chunked', '-K',
        is_flag=True,
        default=False,
        help='When set, perform an incremental PCA on segments of '
        '--chunk-size, which automatically zero centers and ignore settings of '
        '--random-state and --svd-solver.',
    ),
    click.option(
        '--chunk-size', '-Z',
        type=click.INT,
        callback=required_by('--chunked'),
        default=None,
        show_default=True,
        help='Number of observations to include in each chunk, required by '
        '--chunked.',
    ),
]

NEIGHBOR_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    *COMMON_OPTIONS['use_pc'],
    COMMON_OPTIONS['key_added'],
    COMMON_OPTIONS['random_state'],
    click.option(
        '--n-neighbors', '-k',
        type=CommaSeparatedText(click.INT, simplify=True),
        default=15,
        show_default=True,
        help='The size of local neighborhood (in terms of number of '
        'neighboring data points) used for manifold approximation. Larger '
        'values result in more global views of the manifold, while smaller '
        'values result in more local data being preserved. In general values '
        'should be in the range 2 to 100.  If --knn is set, number of nearest '
        'neighbors to be searched, othwise a Gaussian kernel width is set to '
        'the distance of the --n-neighbors neighbor.',
    ),
    click.option(
        '--no-knn', 'knn',
        is_flag=True,
        flag_value=False,
        default=True,
        show_default=True,
        help='When NOT set, use a hard threshold to restrict the number of '
        'neighbors to --n-neighbors. Otherwise, use a Gaussian kernel to '
        'assign low weights to neighbors more distant than the --n-neighbors '
        'nearest neighbor',
    ),
    click.option(
        '--method', '-m',
        type=click.Choice(['umap', 'gauss']),
        default='umap',
        show_default=True,
        help='Use umap or gauss with adaptive width for computing '
        'connectivities.'
    ),
]

UMAP_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    COMMON_OPTIONS['knn_graph'][0], # --use-graph
    COMMON_OPTIONS['init_pos'],
    COMMON_OPTIONS['random_state'],
    COMMON_OPTIONS['key_added'],
    COMMON_OPTIONS['export_embedding'],
    click.option(
        '--init-pos',
        type=click.STRING,
        default='spectral',
        show_default=True,
        help='How to initialize the low dimensional embedding.',
    ),
    click.option(
        '--min-dist',
        type=click.FLOAT,
        default=0.5,
        show_default=True,
        help='The effective minimum distance between embedded points. Smaller '
        'values will result in a more clustered embedding, while larger values '
        'will results in a more even dispersal of points.',
    ),
    click.option(
        '--spread',
        type=click.FLOAT,
        default=1.0,
        show_default=True,
        help='The effective scale of embedded points, which determines the '
        'scale at which embedded points will be spread out.',
    ),
    click.option(
        '--n-components',
        type=click.INT,
        default=2,
        show_default=True,
        help='The number of dimensions of the embedding.',
    ),
    click.option(
        '--maxiter',
        type=click.INT,
        default=None,
        show_default=True,
        help='The number of iterations of the optimization.',
    ),
    click.option(
        '--alpha',
        type=click.FLOAT,
        default=1.0,
        show_default=True,
        help='The initial learning rate for the embedding optimization.',
    ),
    click.option(
        '--gamma',
        type=click.FLOAT,
        default=1.0,
        show_default=True,
        help='Weighting applied to negative samples in low dimensional '
        'embedding optimization.',
    ),
    click.option(
        '--negative-sample-rate',
        type=click.INT,
        default=5,
        show_default=True,
        help='The number of negative edge samples to use per positive edge '
        'sample in optimizing the low dimensional embedding.',
    ),
]

TSNE_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    *COMMON_OPTIONS['use_pc'],
    COMMON_OPTIONS['random_state'],
    COMMON_OPTIONS['key_added'],
    COMMON_OPTIONS['n_jobs'],
    COMMON_OPTIONS['export_embedding'],
    click.option(
        '--perplexity',
        type=click.FLOAT,
        default=30,
        show_default=True,
        help='The perplexity is related to the number of nearest neighbors '
        'that is used in other manifold learning algorithms. Larger datasets '
        'usually require a larger perplexity. Consider selecting a value '
        'between 5 and 50. The choice is not extremely critical since t-SNE '
        'is quite insensitive to this parameter.',
    ),
    click.option(
        '--early-exaggeration',
        type=click.FLOAT,
        default=12,
        show_default=True,
        help='Controls how tight natural clusters in the original space are in '
        'the embedded space and how much space will be between them. For '
        'larger values, the space between natural clusters will be larger in '
        'the embedded space. Again, the choice of this parameter is not very '
        'critical. If the cost function increases during initial optimization, '
        'the early exaggeration factor or the learning rate might be too high.',
    ),
    click.option(
        '--learning-rate',
        type=click.FLOAT,
        default=1000,
        show_default=True,
        help='Note that the R-package "Rtsne" uses a default of 200. The '
        'learning rate can be a critical parameter. It should be between 100 '
        'and 1000. If the cost function increases during initial optimization, '
        'the early exaggeration factor or the learning rate might be too high. '
        'If the cost function gets stuck in a bad local minimum increasing the '
        'learning rate helps sometimes.',
    ),
    click.option(
        '--no-fast-tsne', 'use_fast_tsne',
        is_flag=True,
        flag_value=False,
        default=True,
        show_default=True,
        help='Use the MulticoreTSNE package by D. Ulyanov if it is installed.',
    ),
]

FDG_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    COMMON_OPTIONS['knn_graph'][0], # --use-graph
    COMMON_OPTIONS['init_pos'],
    COMMON_OPTIONS['random_state'],
    COMMON_OPTIONS['key_added'],
    COMMON_OPTIONS['export_embedding'],
    click.option(
        '--layout',
        type=click.Choice(['fa', 'fr', 'grid_fr', 'kk', 'lgl', 'drl', 'rt']),
        default='fa',
        show_default=True,
        help='Name of any valid igraph layout, including "fa" (ForceAtlas2), '
        '"fr" (Fruchterman Reingold), "grid_fr" (Grid Fruchterman Reingold, '
        'faster than "fr"), "kk" (Kamadi Kawai, slower than "fr"), "lgl" '
        '(Large Graph Layout, very fast), "drl" (Distributed Recursive Layout, '
        'pretty fast) and "rt" (Reingold Tilford tree layout).',
    ),
    click.option(
        '--init-pos',
        type=click.STRING,
        default=None,
        show_default=True,
        help='How to initialize the low dimensional embedding. Can be "paga", '
        'or any valid key of `.obsm`.',
    ),
]

LOUVAIN_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    *COMMON_OPTIONS['knn_graph'],
    COMMON_OPTIONS['restrict_to'],
    COMMON_OPTIONS['random_state'],
    COMMON_OPTIONS['key_added'],
    click.option(
        '--flavor',
        type=click.Choice(['vtraag', 'igraph']),
        default='vtraag',
        show_default=True,
        help='Choose between two packages for computing the clustering. '
        '"vtraag" is much powerful, and the default.',
    ),
    click.option(
        '--resolution', '-r',
        type=CommaSeparatedText(click.FLOAT, simplify=True),
        default=1,
        show_default=True,
        help='For the default flavor "vtraag", you can provide a resolution. '
        'Higher resolution means finding more and smaller clusters.',
    ),
]

LEIDEN_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    *COMMON_OPTIONS['knn_graph'],
    COMMON_OPTIONS['restrict_to'],
    COMMON_OPTIONS['random_state'],
    COMMON_OPTIONS['key_added'],
    click.option(
        '--resolution', '-r',
        type=CommaSeparatedText(click.FLOAT, simplify=True),
        default=1,
        show_default=True,
        help='A parameter value controlling the coarseness of the clustering. '
        'Higher values lead to more clusters. Set to "None" if overriding '
        '--partition_type to one that doesn\'t accept `resolution_parameter`.',
    ),
    click.option(
        '--n-iterations',
        type=click.INT,
        default=-1,
        show_default=True,
        help='How many iterations of the Leiden clustering algorithm to '
        'perform. -1 has the algorithm run until it reaches its optimal '
        'clustering.',
    ),
]

DIFFEXP_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    COMMON_OPTIONS['use_raw'],
    click.option(
        '--groupby', '-g',
        type=click.STRING,
        required=True,
        help='The key of the observations grouping to consider.',
    ),
    click.option(
        '--groups',
        type=CommaSeparatedText(simplify=True),
        default='all',
        show_default=True,
        help='Subset of groups to which comparison shall be restricted.',
    ),
    click.option(
        '--reference',
        type=click.STRING,
        default='rest',
        show_default=True,
        help='If "rest", compare each group to the union of the rest of the '
        'groups. If a group identifier, compare with respect to this group.',
    ),
    click.option(
        '--n-genes', '-n',
        type=click.INT,
        default=None,
        show_default=True,
        help='The number of genes that appear in the retured tables. By '
        'default return all available genes depending on the value of '
        '--use-raw.'
    ),
    click.option(
        '--method',
        type=click.Choice(
            ['logreg', 't-test', 'wilcoxon', 't-test_overestim_var']),
        default='t-test_overestim_var',
        show_default=True,
        help='Method of performing differential expression analysis.',
    ),
    click.option(
        '--corr-method',
        type=click.Choice(['benjamini-hochberg', 'bonferroni']),
        default='benjamini-hochberg',
        show_default=True,
        help='P-value correction method. Used only for "t-test", '
        '"t-test_overestim_var" and "wilcoxon".',
    ),
    click.option(
        '--rankby-abs',
        is_flag=True,
        default=False,
        show_default=True,
        help='Rank genes by the absolute value of the score, not by the score. '
        'The returned scores are never the absolute values.',
    ),
    click.option(
        '--filter-params',
        type=Dictionary(keys=[
            'min_in_group_fraction',
            'max_out_group_fraction',
            'min_fold_change',
        ]),
        default=None,
        show_default=True,
        help='Parameters for filtering DE results, valid parameters are: '
        '"min_in_group_fraction" (float), "max_out_group_fraction" (float), '
        '"min_fold_change" (float).',
    ),
    click.option(
        '--logreg-param',
        type=Dictionary(),
        default=None,
        show_default=True,
        help='Parameters passed to `sklearn.linear_model.LogisticRegression`.',
    ),
    click.option(
        '--save',
        type=click.Path(dir_okay=False, writable=True),
        default=None,
        show_default=True,
        help='Tab-separated table to store results of differential expression '
        'analysis.',
    ),
]

PAGA_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    COMMON_OPTIONS['knn_graph'][0], # --use-graph
    COMMON_OPTIONS['key_added'],
    click.option(
        '--groups',
        type=click.STRING,
        required=True,
        help='Key for categorical in `.obs`. You can pass your predefined '
        'groups by choosing any categorical annotation of observations.',
    ),
    click.option(
        '--model',
        type=click.Choice(['v1.2', 'v1.0']),
        default='v1.2',
        show_default=True,
        help='The PAGA connectivity model.',
    ),
]

DIFFMAP_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    COMMON_OPTIONS['knn_graph'][0], # --use-graph
    COMMON_OPTIONS['key_added'],
    COMMON_OPTIONS['export_embedding'],
    COMMON_OPTIONS['n_comps'],
]

DPT_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    COMMON_OPTIONS['knn_graph'][0], # --use-graph
    COMMON_OPTIONS['key_added'],
    click.option(
        '--root',
        type=(click.STRING, click.STRING),
        default=(None, None),
        show_default=True,
        help='Specify a categorical annotaion of observations (`.obs`) and a '
        'value representing the root cells.',
    ),
    click.option(
        '--n-dcs',
        type=click.INT,
        default=10,
        show_default=True,
        help='The number of diffusion components to use.',
    ),
    click.option(
        '--n-branchings',
        type=click.INT,
        default=0,
        show_default=True,
        help='Number of branchings to detect.',
    ),
    click.option(
        '--min-group-size',
        type=click.FLOAT,
        default=0.01,
        show_default=True,
        help='During recursive splitting of branches for --n-branchings > 1, '
        'do not consider branches/groups that contain fewer than this fraction '
        'of the total number of data points.',
    ),
]

PLOT_EMBED_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['plot'],
    click.option(
        '--basis',
        type=click.STRING,
        default='umap',
        show_default=True,
        help='Name of the embedding to plot, must be a key of `.obsm` without '
        'the prefix "X_".',
    ),
    click.option(
        '--color',
        type=CommaSeparatedText(simplify=True),
        default=None,
        show_default=True,
        help='Keys for annotations of observations/cells or variables/genes.',
    ),
    click.option(
        '--use-raw/--no-raw',
        default=None,
        show_default=True,
        help='Use `.raw` attribute for coloring with gene expression. If '
        '`None`, uses `.raw` if present.',
    ),
    click.option(
        '--layers',
        type=click.STRING,
        default=None,
        help='Name of the AnnData object layer to be plotted. By default '
        '`.raw.X` is plotted. If --no-raw, then `.X` is plotted. If --layer '
        'is set to a valid layer name, then the layer is plotted. --layer '
        'takes precedence over --use-raw/--no-raw.'
    ),
    click.option(
        '--legend-loc',
        type=click.Choice(['right margin', 'on data']),
        default='right margin',
        show_default=True,
        help='Location of legend, either "on data", "right margin" or valid '
        'keywords for `matplotlib.legend`.',
    ),
    click.option(
        '--legend-fontsize',
        type=click.INT,
        default=15,
        show_default=True,
        help='Legend font size.',
    ),
]
