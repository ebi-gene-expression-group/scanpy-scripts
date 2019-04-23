"""
Provide cmd options
"""

import click
from .click_utils import (
    CommaSeparatedText,
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
            type=click.Path(),
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
            '--show-obj',
            type=click.Choice(['stdout', 'stderr']),
            default=None,
            show_default=True,
            help='Print output object summary info to specified stream.',
        ),
    ],

    'copy': click.option(
        '--copy',
        is_flag=True,
        default=False,
        show_default=True,
        help='Return a copy instead of writing to supplied input',
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

    'zero_center': click.option(
        '--no-zero-center', 'zero_center',
        is_flag=True,
        flag_value=False,
        default=True,
        help='When set, omit zero-centering variables to allow efficient '
        'handling of sparse input.',
    ),

    'n_pcs': click.option(
        '--n-pcs', '-n',
        type=click.INT,
        default=None,
        show_default=True,
        help='Use this many PCs. Use `.X` if --n-pcs is 0 when --use-rep is '
        'None.',
    ),

    'use_rep': click.option(
        '--use-rep', '-u',
        type=click.STRING,
        default=None,
        show_default=True,
        help='Use the indicated representation. If None, the representation is '
        'chosen automatically: for `.n_vars` < 50, `.X` is used, otherwise '
        '`X_pca` is used. If `X_pca` is not present, it\'s computed with '
        'default parameters.'
    ),

    'n_jobs': click.option(
        '--n-jobs', '-J',
        type=click.INT,
        default=None,
        show_default=True,
        help='Number of jobs',
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
        type=(str, CommaSeparatedText()),
        multiple=True,
        help='Categorical attributes used to filter the data, '
        'in the format of "-c <name> <values>", '
        'where entries with attribute <name> with value in <values> are kept. '
        'Multiple -c entries allowed.',
    ),
    click.option(
        '--subset', '-s',
        type=(str, click.File()),
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
        '--normalize-to', '-t', 'counts_per_cell_after',
        type=float,
        default=10_000,
        show_default=True,
        help='Normalize per cell nUMI to this number.',
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
        multiple=False,
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
    click.option(
        '--keys', '-k',
        type=CommaSeparatedText(simplify=True),
        default=None,
        show_default=True,
        help='Key(s) for observation annotation on which to regress.',
    ),
    click.option(
        '--n-jobs', '-j',
        type=click.INT,
        default=None,
        show_default=True,
        help='Number of jobs for parallel computation.',
    ),
]

PCA_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    click.option(
        '--n-comps', '-n',
        type=click.INT,
        default=50,
        show_default=True,
        help='Number of principal components to compute',
    ),
    click.option(
        '--svd-solver', '-V',
        type=click.Choice(['auto', 'arpack', 'randomized']),
        default='auto',
        show_default=True,
        help='SVD solver to use.'
    ),
    COMMON_OPTIONS['zero_center'],
    COMMON_OPTIONS['random_state'],
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
    COMMON_OPTIONS['n_pcs'],
    COMMON_OPTIONS['use_rep'],
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
    COMMON_OPTIONS['key_added'],
    COMMON_OPTIONS['random_state'],
    COMMON_OPTIONS['copy'],
]

UMAP_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
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
    click.option(
        '--init-pos',
        type=click.STRING,
        default='spectral',
        show_default=True,
        help='How to initialize the low dimensional embedding.',
    ),
]

TSNE_CMD_OPTIONS = [
    *COMMON_OPTIONS['input'],
    *COMMON_OPTIONS['output'],
    COMMON_OPTIONS['n_pcs'],
    COMMON_OPTIONS['use_rep'],
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
    COMMON_OPTIONS['random_state'],
    click.option(
        '--no-fast-tsne', 'use_fast_tsne',
        is_flag=True,
        flag_value=False,
        default=True,
        show_default=True,
        help='Use the MulticoreTSNE package by D. Ulyanov if it is installed.',
    ),
    COMMON_OPTIONS['n_jobs'],
]
