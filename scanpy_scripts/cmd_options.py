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
            '--loom-write-obsm-varm', '-b',
            is_flag=True,
            default=False,
            show_default=True,
            help='Write obsm and varm to the Loom file?',
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

    'save': [
        click.option(
            '--save-raw', '-r',
            is_flag=True,
            default=False,
            show_default=True,
            help='Save adata to adata.raw before processing.',
        ),
        click.option(
            '--save-layer', '-y',
            type=click.STRING,
            default=None,
            show_default=True,
            help='Save adata.X to the specified layer before processing.',
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

    'frame_title': [
        click.option(
            '--frameon/--frameoff', 'frameon',
            default=True,
            show_default=True,
            help='Draw a frame around the plot',
        ),
        click.option(
            '--title',
            type=CommaSeparatedText(simplify=True),
            default=None,
            show_default=True,
            help='Provide title for the plot or panels.',
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
            '--neighbors-key',
            type=click.STRING,
            default=None,
            show_default=False,
            help='If not specified, look in .uns[‘neighbors’] for neighbors '
            'settings and .obsp[‘connectivities’], .obsp[‘distances’] for connectivities and '
            'distances respectively (default storage places for pp.neighbors). If specified, '
            'look in .uns[neighbors_key] for neighbors settings and '
            '.obsp[.uns[neighbors_key][‘connectivities_key’]], '
            '.obsp[.uns[neighbors_key][‘distances_key’]] for connectivities and distances '
            'respectively.'
        ),
        click.option(
            '--obsp',
            type=click.STRING,
            default=None,
            show_default=True,
            help='Use .obsp[obsp] as adjacency. You can’t specify both obsp and '
            'neighbors_key at the same time.'
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

    'neighbor_metric': click.option(
        '--metric', '-t',
        type=click.Choice(['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan', 'braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']),
        default='euclidean',
        show_default=True,
        help='A known metric’s name.'
    ),

    'layer':click.option(
        '--layer',
        type=CommaSeparatedText(simplify=True),
        default=None,
        show_default=True,
        help='Name of the AnnData object layer that wants to be plotted. By '
        'default adata.raw.X is plotted. If use_raw=False is set, then adata.X '
        'is plotted. If layer is set to a valid layer name, then the layer is '
        'plotted. layer takes precedence over use_raw.',
    ),

    'n_comps': click.option(
        '--n-comps',
        type=click.INT,
        default=None,
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
        default=None,
        show_default=True,
        help='Use expression values in `.raw` if present.',
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

    'export_cluster': click.option(
        '--export-cluster',
        type=click.Path(dir_okay=False, writable=True),
        default=None,
        show_default=True,
        help='Export embeddings in a tab-separated text table.',
    ),

    'var_names': click.option(
            '--var-names',
            type=(CommaSeparatedText()),
            show_default=True,
            help='var_names should be a valid subset of adata.var_names.',
    ),

    'gene_symbols': click.option(
            '--gene-symbols',
            type=CommaSeparatedText(simplify=True),
            default=None,
            show_default=True,
            help='Column name in .var DataFrame that stores gene symbols. By '
            'default this is assumed to be the index column of the .var '
            'DataFrame. Setting this option allows alternative names to be '
            'used.',
    ),

    'diffexp_plot': [
        click.option(
            '--rgg',
            is_flag=True,
            default=False,
            show_default=True,
            help='When set, use the rank_genes_groups_ form of the function, '
            'where gene lists are automatically selected.',
        ),
        click.option(
             '--groupby',
            type=CommaSeparatedText(simplify=True),
            default=None,
            show_default=True,
            help='The key of the observation grouping to consider.',
        ),
        click.option(
            '--log',
            is_flag=True,
            default=False,
            show_default=True,
            help='Plot on logarithmic axis.',
        ),
        click.option(
            '--num-categories',
            type=click.INT,
            default=7,
            show_default=True,
            help='Only used if groupby observation is not categorical. This value '
            'determines the number of groups into which the groupby observation '
            'should be subdivided.',
        ),
        click.option(
            '--dendrogram',
            is_flag=True,
            default=False,
            show_default=False,
            help='If True, a dendrogram based on the hierarchical clustering '
            'between the groupby categories is added. The dendrogram information is '
            'computed using scanpy.tl.dendrogram(). If tl.dendrogram has not been '
            'called previously the function is called with default parameters.',
        ),
        click.option(
            '--standard-scale',
            type=click.Choice(['var', 'obs']),
            default=None,
            show_default=True,
            help='Whether or not to standardize that dimension between 0 and 1, '
            'meaning for each variable or group, subtract the minimum and divide '
            'each by its maximum.'
        ),

    ],

    'sviol': [
        click.option(
            '--no-stripplot', 'stripplot',
            is_flag=True,
            default=True,
            show_default=True,
            help='When set, do not add a stripplot on top of the violin plot.',
        ),
        click.option(
            '--no-jitter', 'jitter',
            is_flag=True,
            default=True,
            show_default=True,
            help='Suppress jitter in the stripplot (only when stripplot is True)'
        ),
        click.option(
            '--size',
            type=click.INT,
            default=1,
            show_default=True,
            help='Size of the jitter points.'
        ),
        click.option(
            '--order',
            type=CommaSeparatedText(),
            default=None,
            show_default=True,
            help='Order in which to show the categories.'
        ),
        click.option(
            '--scale',
            type=click.Choice(['area', 'count', 'width']),
            default='width',
            show_default=True,
            help='The method used to scale the width of each violin. If ‘area’, '
            'each violin will have the same area. If ‘count’, the width of the '
            'violins will be scaled by the number of observations in that bin. If '
            '‘width’, each violin will have the same width.'
        ),
        click.option(
            '--row-palette',
            type=CommaSeparatedText(simplify=True),
            default='muted',
            show_default=True,
            help='The row palette determines the colors to use in each of the '
            'stacked violin plots. The value should be a valid seaborn palette name '
            'or a valic matplotlib colormap (see '
            'https://seaborn.pydata.org/generated/seaborn.color_palette.html). '
            'Alternatively, a single color name or hex value can be passed. E.g. '
            '‘red’ or ‘#cc33ff’.'
        ),
    ],


    'dot': [
        click.option(
            '--expression-cutoff',
            type=click.FLOAT,
            default=0,
            show_default=True,
            help='Expression cutoff that is used for binarizing the gene expression '
            'and determining the fraction of cells expressing given genes. A gene is '
            'expressed only if the expression value is greater than this threshold.'
        ),
        click.option(
            '--mean-only-expressed',
            is_flag=True,
            default=False,
            show_default=True,
            help='If True, gene expression is averaged only over the cells '
            'expressing the given genes.',
        ),
        click.option(
            '--color-map',
            type=CommaSeparatedText(simplify=True),
            default='Reds',
            show_default=True,
            help='String denoting matplotlib color map.',
        ),
        click.option(
            '--dot-max',
            type=click.FLOAT,
            default=None,
            show_default=True,
            help='If none, the maximum dot size is set to the maximum fraction '
            'value found (e.g. 0.6). If given, the value should be a number between '
            '0 and 1. All fractions larger than dot_max are clipped to this value.'
        ),
        click.option(
            '--dot-min',
            type=click.FLOAT,
            default=None,
            show_default=True,
            help='If none, the minimum dot size is set to 0. If given, the value '
            'should be a number between 0 and 1. All fractions smaller than dot_min '
            'are clipped to this value.'
        ),
        click.option(
            '--smallest-dot',
            type=click.FLOAT,
            default=0,
            show_default=True,
            help='If none, the smallest dot has size 0. All expression levels with '
            'dot_min are potted with smallest_dot dot size.'
        ),
    ],

    'heat': [
        click.option(
            '--show-gene-labels',
            is_flag=True,
            default=None,
            show_default=True,
            help='By default gene labels are shown when there are 50 or less '
            'genes. Otherwise the labels are removed.'
        ),
    ],

    'swap_axes': click.option(
        '--swap-axes',
        is_flag=True,
        default=False,
        show_default=True,
        help='By default, the x axis contains var_names (e.g. genes) and the y '
        'axis the groupby categories. By setting swap_axes then x are the '
        'groupby categories and y the var_names. When swapping axes '
        'var_group_positions are no longer used.',
    ),

    'rank_genes_groups_plots': [
        click.option(
            '--groups',
            type=CommaSeparatedText(),
            default=None,
            show_default=True,
            help='The groups for which to show the gene ranking.'
        ),
        click.option(
            '--n-genes', '-n',
            type=click.INT,
            default=10,
            show_default=True,
            help='Number of genes to show.'
        ),
    ],

    'root': click.option(
        '--root',
        type=click.INT,
        default=0,
        show_default=True,
        help='If choosing a tree layout, this is the index of the root node.',
    ),

    'plot_embed': [
        click.option(
            '--use-raw/--no-raw',
            default=None,
            show_default=True,
            help='Use `.raw` attribute for coloring with gene expression. If '
            '`None`, uses `.raw` if present.',
        ),
        click.option(
            '--groups',
            type=click.STRING,
            default=None,
            help='Key for categorical in `.obs`. You can pass your predefined '
            'groups by choosing any categorical annotation of observations.',
        ),
    ],

    'batch_key': click.option(
        '--batch-key', 'key',
        type=click.STRING,
        required=True,
        help='The name of the column in adata.obs that differentiates among '
        'experiments/batches.'
    ),
    
    'batch_layer': click.option(
        '--layer', '-l',
        type=click.STRING,
        default=None,
        show_default=True,
        help="Layer to batch correct. By default corrects the contents of .X."
    ),

    'scrublet': [
        click.option(
            '--sim-doublet-ratio',
            type=click.FLOAT,
            default=2.0,
            show_default=True,
            help='Number of doublets to simulate relative to the number of '
            'observed transcriptomes.',
        ),
        click.option(
            '--synthetic-doublet-umi-subsampling',
            type=click.FLOAT,
            default=1.0,
            show_default=True,
            help='Where input_obj_sim not suplied, rate for sampling UMIs when '
            'creating synthetic doublets. If 1.0, each doublet is created by '
            'simply adding the UMI counts from two randomly sampled observed '
            'transcriptomes.  For values less than 1, the UMI counts are added '
            'and then randomly sampled at the specified rate.'
        ),
    ],
}

COMMON_OPTIONS['opt_output'] = [
    click.option(
        '--output-obj',
        type=click.Path(dir_okay=False, writable=True),
        help='Optionally output an object to the specified path.',
    ),
    *COMMON_OPTIONS['output'][1:], 
]

CMD_OPTIONS = {
    'read': [
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
    ],

    'filter': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        COMMON_OPTIONS['save'][0], # --save-raw
        click.option(
            '--gene-name', '-g',
            type=click.STRING,
            default='index',
            show_default=True,
            help='Name of the variable that contains gene names, used for flagging '
            'mitochondria genes when column "mito" is absent from `.var`.',
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
            'If <values> is preceded by "!", entries with value in <values> are '
            'removed. Multiple -c entries allowed.',
        ),
        click.option(
            '--subset', '-s',
            type=(click.STRING, click.File()),
            multiple=True,
            help='Similar to --category in the format of "-s <name> <file>", '
            'but the <file> to be a one-column table that provides the values. '
            'Multiple -s entries allowed.',
        ),
        click.option(
            '--force-recalc',
            is_flag=True,
            default=False,
            help='When set, re-calculate `pct_counts_<qc_variable>` and '
            '`pct_counts_in_top_<n>_genes` even if they exist.',
        ),
    ],

    'norm': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        *COMMON_OPTIONS['save'],
        COMMON_OPTIONS['key_added'],
        click.option(
            '--no-log-transform', 'log_transform',
            is_flag=True,
            default=True,
            show_default=True,
            help='When set, do not apply (natural) log transform following normalisation.',
        ),
        click.option(
            '--normalize-to', '-t', 'target_sum',
            type=float,
            default=10_000,
            show_default=True,
            help='Normalize per cell nUMI to this number.',
        ),
        click.option(
            '--exclude-highly-expressed', '-e', 'exclude_highly_expressed',
            is_flag=True,
            default=False,
            show_default=True,
            help='Exclude (very) highly expressed genes for the computation of '
            'the normalization factor (size factor) for each cell. A gene is considered '
            'highly expressed, if it has more than max_fraction of the total counts in at '
            'least one cell. The not-excluded genes will sum up to the number '
            'specified by --normalize-to.'
        ),
        click.option(
            '--max-fraction', '-m', 'max_fraction',
            type=float,
            default=0.05,
            show_default=True,
            help='If exclude_highly_expressed=True, consider cells as highly '
            'expressed that have more counts than max_fraction of the original total counts '
            'in at least one cell.'
        ),
        click.option(
            '--layers', '-l',
            type=CommaSeparatedText(simplify=True),
            default=None,
            show_default=True,
            help="List of layers to normalize. Set to 'all' to normalize all layers."
        ),
        click.option(
            '--layer-norm', '-n', 'layer_norm',
            type=click.Choice(['after', 'X']),
            default=None,
            show_default=True,
            help="Specifies how to normalize layers: 1) If None, after "
            "normalization, for each layer in layers each cell has a total count equal to "
            "the median of the counts_per_cell before normalization of the layer. 2) If "
            "'after', for each layer in layers each cell has a total count equal to "
            "target_sum. 3) If 'X', for each layer in layers each cell has a total count "
            "equal to the median of total counts for observations (cells) of adata.X before "
            "normalization.'"
        ),
    ],

    'hvg': [
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
            '--span',
            type=click.FLOAT,
            default=0.3,
            show_default=True,
            help="The fraction of the data (cells) used when estimating the "
            "variance in the loess model fit if flavor='seurat_v3'."    
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
            default=None,
            show_default=True,
            help='Number of highly-variable genes to keep.',
        ),
        click.option(
            '--flavor', '-v',
            type=click.Choice(['seurat', 'cell_ranger', 'seurat_v3']),
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
            '--batch-key', 'batch_key',
            type=click.STRING,
            default=None,
            help="If specified, highly-variable genes are selected within each "
            "batch separately and merged. This simple process avoids the selection of "
            "batch-specific genes and acts as a lightweight batch correction method. For all "
            "flavors, genes are first sorted by how many batches they are a HVG. For "
            "dispersion-based flavors ties are broken by normalized dispersion. If flavor = "
            "'seurat_v3', ties are broken by the median (across batches) rank based on "
            "within-batch normalized variance."
        ),
    ],

    'scale': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        *COMMON_OPTIONS['save'],
        COMMON_OPTIONS['zero_center'],
        click.option(
            '--max-value', '-m',
            type=click.FLOAT,
            default=None,
            show_default=True,
            help='When specified, clip to this value after scaling, otherwise do '
            'not clip',
        ),
        click.option(
            '--layer', '-l',
            type=CommaSeparatedText(simplify=True),
            default=None,
            help="If provided, which element of layers to scale."
        ),
    ],

    'regress': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        *COMMON_OPTIONS['save'],
        COMMON_OPTIONS['n_jobs'],
        click.option(
            '--keys', '-k',
            type=CommaSeparatedText(simplify=True),
            default=None,
            show_default=True,
            help='Key(s) for observation annotation on which to regress.',
        ),
    ],

    'pca': [
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
    ],

    'neighbor': [
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
            type=click.Choice(['umap', 'gauss', 'rapids']),
            default='umap',
            show_default=True,
            help='Use umap or gauss with adaptive width for computing '
            'connectivities. Use rapids for the RAPIDS implementation of UMAP '
            '(experimental, GPU only).'
        ),
        COMMON_OPTIONS['neighbor_metric'],
    ],

    'umap': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        COMMON_OPTIONS['knn_graph'][0], # --neighbors-key
        COMMON_OPTIONS['random_state'],
        COMMON_OPTIONS['key_added'],
        COMMON_OPTIONS['export_embedding'],
        click.option(
            '--init-pos',
            type=click.STRING,
            default='spectral',
            show_default=True,
            help='How to initialize the low dimensional embedding. Can be '
            '"spectral", "paga" or "random", or any key of `.obsm`.',
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
        click.option(
            '--method',
            type=click.Choice(['umap', 'rapids']),
            default='umap',
            show_default=True,
            help='Use the original ‘umap’ implementation, or ‘rapids’ '
            '(experimental, GPU only).'
        ),
    ],

    'tsne': [
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
            help='When NOT set, use the MulticoreTSNE package by D. Ulyanov if '
            'installed.',
        ),
    ],

    'fdg': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        COMMON_OPTIONS['random_state'],
        COMMON_OPTIONS['export_embedding'],
        COMMON_OPTIONS['root'],
        click.option(
            '--init-pos',
            type=click.STRING,
            default=None,
            help='Use precomputed coordinates for initialization. Can be any key '
            'of `.obsm` or "paga" if .uns["paga"] is present',
        ),
        click.option(
            '--layout',
            type=click.Choice(['fa', 'fr', 'grid_fr', 'kk', 'lgl', 'drl', 'rt', 'rt_circular']),
            default='fa',
            show_default=True,
            help='Name of any valid igraph layout, including "fa" (ForceAtlas2), '
            '"fr" (Fruchterman Reingold), "grid_fr" (Grid Fruchterman Reingold, '
            'faster than "fr"), "kk" (Kamadi Kawai, slower than "fr"), "lgl" '
            '(Large Graph Layout, very fast), "drl" (Distributed Recursive Layout, '
            'pretty fast) and "rt" (Reingold Tilford tree layout).',
        ),
        click.option(
            '--key-added-ext',
            type=click.STRING,
            default=None,
            show_default=True,
            help="By default, append 'layout'"
        ),
        click.option(
            '--init-pos',
            type=click.STRING,
            default=None,
            show_default=True,
            help='How to initialize the low dimensional embedding. Can be "paga", '
            'or any valid key of `.obsm`.',
        ),
        COMMON_OPTIONS['knn_graph'][0], # --neighbors-key
        COMMON_OPTIONS['knn_graph'][1], # --obsp
    ],

    'louvain': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        COMMON_OPTIONS['export_cluster'],
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
    ],

    'leiden': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        COMMON_OPTIONS['export_cluster'],
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
    ],

    'diffexp': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        COMMON_OPTIONS['use_raw'],
        COMMON_OPTIONS['key_added'],
        click.option(
            '--layer', '-l',
            type=click.STRING,
            default=None,
            help="Key from adata.layers whose value will be used to perform tests on."
        ),
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
            '--pts',
            is_flag=True,
            default=False,
            show_default=True,
            help='Compute the fraction of cells expressing the genes.'
        ),
        click.option(
            '--tie-correct',
            is_flag=True,
            default=False,
            show_default=True,
            help="Use tie correction for 'wilcoxon' scores. Used only for "
            "'wilcoxon'."
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
    ],

    'paga': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        COMMON_OPTIONS['knn_graph'][0], # --neighbors-key
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
        click.option(
            '--use-rna-velocity',
            is_flag=True,
            default=False,
            show_default=True,
            help='Use RNA velocity to orient edges in the abstracted graph and '
            'estimate transitions. Requires that adata.uns contains a directed single-cell '
            'graph with key velocity_graph. This feature might be subject to change in the '
            'future.',
        ),
    ],

    'diffmap': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        COMMON_OPTIONS['knn_graph'][0], # --neighbors-key
        COMMON_OPTIONS['key_added'],
        COMMON_OPTIONS['export_embedding'],
        COMMON_OPTIONS['n_comps'],
    ],

    'dpt': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        COMMON_OPTIONS['knn_graph'][0], # --neighbors-key
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
        click.option(
            '--disallow-kendall-tau-shift', 'allow_kendall_tau_shift',
            is_flag=True,
            default=True,
            show_default=True,
            help='By default: If a very small branch is detected upon '
            'splitting, shift away from maximum correlation in Kendall tau criterion of '
            '[Haghverdi16] to stabilize the splitting. Use flag to disable this.'
        ),
    ],

    'combat': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        COMMON_OPTIONS['batch_key'],
        COMMON_OPTIONS['batch_layer'],
        click.option(
            '--key-added',
            type=click.STRING,
            default=None,
            show_default=True,
            help="Key under which to add the computed results. By default a new "
            "layer will be created called 'combat', 'combat_{layer}' or "
            "'combat_layer_{key_added}' where those parameters were specified. A value of 'X' "
            "causes batch-corrected values to overwrite the original content of .X."
        ),
        click.option(
            '--covariates',
            type=(CommaSeparatedText()),
            default=None,
            show_default=True,
            help="Comma-separated list of additional covariates besides the "
            "batch variable such as adjustment variables or biological condition. This "
            "parameter refers to the design matrix X in Equation 2.1 in [Johnson07] and to "
            "the mod argument in the original combat function in the sva R package.  Note "
            "that not including covariates may introduce bias or lead to the removal of "
            "biological signal in unbalanced designs."
        ),
        
    ],

    'harmony': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        COMMON_OPTIONS['batch_key'],
        click.option(
            '--basis',
            type=click.STRING,
            default='X_pca',
            show_default=True,
            help="The name of the field in adata.obsm where the PCA table is "
            "stored. Defaults to 'X_pca', which is the default for sc.tl.pca()."
        ),
        click.option(
            '--adjusted-basis',
            type=click.STRING,
            default='X_pca_harmony',
            show_default=True,
            help='The name of the field in adata.obsm where the adjusted PCA '
            'table will be stored after running this function.'
        ),
        click.option(
            '--theta',
            type=click.FLOAT,
            default=2,
            show_default=True,
            help='Diversity clustering penalty parameter. theta=0 does not encourage any '
            'diversity. Larger values of theta result in more diverse clusters.'
        ),
        click.option(
            '--lambda', 'lamb',
            type=click.FLOAT,
            default=1,
            show_default=True,
            help='Ridge regression penalty parameter. Lambda must be strictly '
            'positive.  Smaller values result in more aggressive correction.'
        ),
        click.option(
            '--sigma',
            type=click.FLOAT,
            default=0.1,
            show_default=True,
            help='Width of soft kmeans clusters. Sigma scales the distance from '
            'a cell to cluster centroids. Larger values of sigma result in cells assigned to '
            'more clusters. Smaller values of sigma make soft kmeans cluster approach hard '
            'clustering.'
        ),
        click.option(
            '--n-clust', 'nclust',
            type=click.INT,
            default=None,
            show_default=False,
            help='Number of clusters in model. nclust=1 equivalent to simple '
            'linear regression.'
        ),
        click.option(
            '--tau',
            type=click.INT,
            default=0,
            show_default=True,
            help='Protection against overclustering small datasets with large ones. '
            'tau is the expected number of cells per cluster.'
        ),
        click.option(
            '--block-size',
            type=click.FLOAT,
            default=0.05,
            show_default=True,
            help='What proportion of cells to update during clustering. Between '
            '0 to 1, default 0.05. Larger values may be faster but less accurate.'
        ),
        click.option(
            '--max-iter-cluster', 'max_iter_kmeans',
            type=click.INT,
            default=20,
            show_default=True,
            help='Maximum number of rounds to run clustering at each round of '
            'Harmony.'
        ),
        click.option(
            '--max-iter-harmony',
            type=click.INT,
            default=10,
            show_default=True,
            help='Maximum number of rounds to run Harmony. One round of Harmony '
            'involves one clustering and one correction step.'
        ),
        click.option(
            '--epsilon-cluster',
            type=click.FLOAT,
            default=1e-5,
            show_default=True,
            help='Convergence tolerance for clustering round of Harmony Set to '
            '-Inf to never stop early.'
        ),
        click.option(
            '--epsilon-harmony',
            type=click.FLOAT,
            default=1e-5,
            show_default=True,
            help='Convergence tolerance for clustering round of Harmony Set to '
            '-Inf to never stop early.'
        ),
        COMMON_OPTIONS['random_state'],
    ],

    'mnn': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        *COMMON_OPTIONS['save'],
        COMMON_OPTIONS['batch_key'],
        COMMON_OPTIONS['batch_layer'],
        click.option(
            '--key-added',
            type=click.STRING,
            default=None,
            show_default=True,
            help="Key under which to add the computed results. By default a new "
            "layer will be created called 'mnn', 'mnn_{layer}' or "
            "'mnn_layer_{key_added}' where those parameters were specified. A value of 'X' "
            "causes batch-corrected values to overwrite the original content of .X."
        ),
        click.option(
            '--var-subset',
            type=(click.STRING, CommaSeparatedText()),
            multiple=True,
            help="The subset of vars (list of str) to be used when performing "
            "MNN correction in the format of '--var-subset <name> <values>'. Typically, use "
            "the highly variable genes (HVGs) like '--var-subset highly_variable True'. When "
            "unset, uses all vars."
        ),
        click.option(
            '--n-neighbors', '-k',
            type=CommaSeparatedText(click.INT, simplify=True),
            default=20,
            show_default=True,
            help='Number of mutual nearest neighbors.'
        ),
        click.option(
            '--sigma',
            type=click.FLOAT,
            default=1.0,
            show_default=True,
            help='The bandwidth of the Gaussian smoothing kernel used to '
            'compute the correction vectors.'
        ),
        click.option(
            '--no-cos_norm_in', 'cos_norm_in',
            is_flag=True,
            default=True,
            help='Default behaviour is to perform cosine normalization on the '
            'input data prior to calculating distances between cells. Use this '
            'flag to disable that behaviour.'
        ),
        click.option(
            '--no-cos_norm_out', 'cos_norm_out',
            is_flag=True,
            default=True,
            help='Default behaviour is to perform cosine normalization prior to '
            'computing corrected expression values. Use this flag to disable that '
            'behaviour.'
        ),
        click.option(
            '--svd-dim',
            type=click.INT,
            default=None,
            show_default=True,
            help='The number of dimensions to use for summarizing biological '
            'substructure within each batch. If not set, biological components '
            'will not be removed from the correction vectors.'
        ),
        click.option(
            '--no-var-adj',
            is_flag=True,
            default=True,
            help='Default behaviour is to adjust variance of the correction '
            'vectors. Use this flag to disable that behaviour. Note this step takes most '
            'computing time.'
        ),
        click.option(
            '--compute-angle',
            is_flag=True,
            default=False,
            help='When set, compute the angle between each cell’s correction '
            'vector and the biological subspace of the reference batch.'
        ),
        click.option(
            '--svd-mode',
            type=click.Choice(['svd', 'rsvd', 'irlb']),
            default='rsvd',
            show_default=True,
            help="'svd' computes SVD using a non-randomized SVD-via-ID "
            "algorithm, while 'rsvd' uses a randomized version. 'irlb' performs truncated "
            "SVD by implicitly restarted Lanczos bidiagonalization (forked from "
            "https://github.com/airysen/irlbpy)."
        ),
    ],

    'bbknn': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        COMMON_OPTIONS['key_added'],
        COMMON_OPTIONS['batch_key'],
        click.option(
            '--use-rep', '-u',
            type=click.STRING,
            default='X_pca',
            show_default=True,
            help='The dimensionality reduction in .obsm to use for neighbour '
            'detection.'
        ),
        COMMON_OPTIONS['use_pc'][0], # --n-pcs
        click.option(
            '--no-approx', 'approx',
            is_flag=True,
            default=True,
            help='Default behaviour is to use annoy’s approximate neighbour '
            'finding. This results in a quicker run time for large datasets while also '
            'potentially increasing the degree of batch correction. Use this flag to disable '
            'that behaviour.',
        ),
        COMMON_OPTIONS['neighbor_metric'],
        click.option(
            '--neighbors-within-batch',
            type=click.INT,
            default=3,
            show_default=True,
            help='How many top neighbours to report for each batch; total '
            'number of neighbours will be this number times the number of batches.'
        ),
        click.option(
            '--trim',
            type=click.INT,
            default=None,
            show_default=True,
            help='Trim the neighbours of each cell to these many top '
            'connectivities. May help with population independence and improve the tidiness '
            'of clustering. The lower the value the more independent the individual '
            'populations, at the cost of more conserved batch effect. If None, sets the '
            'parameter value automatically to 10 times the total number of neighbours for '
            'each cell. Set to 0 to skip.'
        ),
        click.option(
            '--annoy-n-trees',
            type=click.INT,
            default=10,
            show_default=True,
            help='Only used when approx=True. The number of trees to construct '
            'in the annoy forest. More trees give higher precision when querying, at the '
            'cost of increased run time and resource intensity.'
        ),
        click.option(
            '--no-use-faiss', 'use_faiss',
            is_flag=True,
            default=True,
            help='Default behaviour If approx=False and the metric is '
            '“euclidean”, is to use the faiss package to compute nearest neighbours if '
            'installed. This improves performance at a minor cost to numerical precision as '
            'faiss operates on float32. Use this flag to disable that behaviour.'
        ),
        click.option(
            '--set-op-mix-ratio',
            type=click.FLOAT,
            default=1,
            show_default=True,
            help='UMAP connectivity computation parameter, float between 0 and '
            '1, controlling the blend between a connectivity matrix formed exclusively from '
            'mutual nearest neighbour pairs (0) and a union of all observed neighbour '
            'relationships with the mutual pairs emphasised (1).'
        ),
        click.option(
            '--local-connectivity',
            type=click.INT,
            default=1,
            show_default=True,
            help='UMAP connectivity computation parameter, how many nearest '
            'neighbors of each cell are assumed to be fully connected (and given a '
            'connectivity value of 1)'
        ),
    ],
    
    'scrublet': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        click.option(
            '--input-obj-sim', 'adata_sim',
            type=click.Path(exists=True, dir_okay=False),
            default=None,
            help='(Advanced use case) Optional annData object generated by '
            'sc.external.pp.scrublet_simulate_doublets(), with same number of  '
            'vars as adata. This should have been built from input_obj after '
            'filtering genes and cells and selcting highly-variable genes.'
        ),
        click.option(
            '--threshold',
            type=click.FLOAT,
            default=None,
            show_default=True,
            help='Doublet score threshold for calling a transcriptome a '
            'doublet. If not set, this is set automatically by looking for the '
            'minimum between the two modes of the doublet_scores_sim_ histogram. '
            'It is best practice to check the threshold visually using the '
            'doublet_scores_sim_ histogram and/or based on co-localization of '
            'predicted doublets in a 2-D embedding.'
        ),
        *COMMON_OPTIONS['scrublet'],
        click.option(
            '--expected-doublet-rate',
            type=click.FLOAT,
            default=0.05,
            show_default=True,
            help='Where input_obj_sim not suplied, the estimated doublet rate '
            'for the experiment.'
        ),
        click.option(
            '--stdev-doublet-rate',
            type=click.FLOAT,
            default=0.02,
            show_default=True,
            help='Where input_obje_sim not suplied, uncertainty in the expected '
            'doublet rate.'
        ),
         click.option(
            '--knn-dist-metric', '-t',
            type=click.Choice(['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan', 'braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']),
            default='euclidean',
            show_default=True,
            help='A known metric’s name.'
        ),
        click.option(
            '--no-normalize-variance', 'normalize_variance',
            is_flag=True,
            default=True,
            help='Default is to normalize the data such that each gene has a '
            'variance of 1. sklearn.decomposition.TruncatedSVD will be used for '
            'dimensionality reduction, if --no-mean-center is set. Use this flag '
            'to disable that behaviour.'
        ),
        click.option(
            '--log-transform',
            is_flag=True,
            default=False,
            show_default=True,
            help='Whether to use :func:~scanpy.pp.log1p to log-transform the '
            'data prior to PCA.'
        ),
        click.option(
            '--no-mean-center', 'mean_center',
            is_flag=True,
            default=True,
            help='If True, center the data such that each gene has a mean of 0. '
            'sklearn.decomposition.PCA will be used for dimensionality '
            'reduction.'
        ),
        click.option(
            '--n-pcs', 'n_prin_comps',
            type=click.INT,
            default=30,
            show_default=True,
            help='Number of principal components used to embed the '
            'transcriptomes prior to k-nearest-neighbor graph construction.'
        ),
        click.option(
            '--no-approx', 'use_approx_neighbors',
            is_flag=True,
            default=True,
            help='Default behaviour is to use the approximate nearest neighbor '
            'method (annoy) for the KNN classifier. Use this flag to disable '
            'that behaviour.'
        ),
        click.option(
            '--get-doublet-neighbor-parents',
            is_flag=True,
            default=False,
            show_default=True,
            help='If set, return (in .uns) the parent transcriptomes that '
            'generated the doublet neighbors of each observed transcriptome. '
            'This information can be used to infer the cell states that '
            'generated a given doublet state.'
        ),
        click.option(
            '--n-neighbors', '-k',
            type=CommaSeparatedText(click.INT, simplify=True),
            default=None,
            show_default=True,
            help='Number of neighbors used to construct the KNN graph of '
            'observed transcriptomes and simulated doublets. If not set, this is '
            'automatically set to np.round(0.5 * np.sqrt(n_obs)).'
        ),
        click.option(
            '--filter', 'filter',
            is_flag=True,
            default=False,
            help='By default, the output object is annotated but not filtered '
            'according to the scrublet status. Setting this flag will cause '
            'predicted multiplet elements to be removed.'
        ),
        click.option(
            '--no-verbose', 'verbose',
            is_flag=True,
            default=True,
            help='Default behaviour is to print progress updates. Use this flag '
            'to disable that.'
        ),
        click.option(
            '--export-table',
            type=click.Path(dir_okay=False, writable=True),
            default=None,
            show_default=True,
            help='Export a table of double scores and calls to the specified file.',
        ),
        COMMON_OPTIONS['random_state'],
    ],
    
    'plot_scrublet': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['plot'],
         click.option(
            '--scale-hist-obs', '-b',
            type=click.Choice(['linear', 'log', 'symlog', 'logit']),
            default='log',
            show_default=True,
            help='Set y axis scale transformation in matplotlib for the plot of observed transcriptomes.'
        ),
         click.option(
            '--scale-hist-sim', '-s',
            type=click.Choice(['linear', 'log', 'symlog', 'logit']),
            default='linear',
            show_default=True,
            help='Set y axis scale transformation in matplotlib for the plot of observed transcriptomes.'
        ),
     ],

    'scrublet_simulate_doublets': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['output'],
        *COMMON_OPTIONS['scrublet'],
        click.option(
            '--layer', '-l',
            type=click.STRING,
            default=None,
            help="Layer of adata where raw values are stored, or ‘X’ if values "
            "are in .X."
        ),
    ],

    'embed': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['plot'],
        *COMMON_OPTIONS['frame_title'],
        COMMON_OPTIONS['layer'],
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
        click.option(
            '--size',
            type=click.FLOAT,
            default=None,
            show_default=True,
            help='Point size. Automatically computed if not specified.',
        ),
        COMMON_OPTIONS['gene_symbols'],
        click.option(
            '--edges',
            is_flag=True,
            default=False,
            show_default=True,
            help='Show edges.',
        ),
        click.option(
            '--edges-width',
            type=click.FLOAT,
            default=0.1,
            show_default=True,
            help='Width of edges.',
        ),
        click.option(
            '--edges-color',
            type=click.STRING,
            default=None,
            show_default=True,
            help='Color of edges. See draw_networkx_edges().',
        ),
        COMMON_OPTIONS['knn_graph'][0], # --neighbors-key
        click.option(
            '--no-sort-order', 'sort_order',
            is_flag=True,
            default=True,
            show_default=True,
            help='Disable default behaviour: for continuous annotations used as '
            'color parameter, plot data points with higher values on top of others.',
        ),
        *COMMON_OPTIONS['plot_embed'],
        click.option(
            '--components',
            type=click.STRING,
            default=None,
            show_default=True,
            help="For instance, ['1,2', '2,3']. To plot all available components use 'all'.",
        ),
        click.option(
            '--projection',
            type=click.Choice(['2d', '3d']),
            default='2d',
            show_default=True,
            help="Projection of plot."
        ),
         
    ],

    'plot_paga': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['plot'],
        *COMMON_OPTIONS['frame_title'],
        *COMMON_OPTIONS['plot_embed'],
        COMMON_OPTIONS['random_state'],
        click.option(
            '--use-key',
            type=click.STRING,
            default='paga',
            show_default=True,
            help='The key in `.uns` that contains trajectory information.',
        ),
        click.option(
            '--layout',
            type=click.Choice(['fa', 'fr', 'grid_fr', 'kk', 'lgl', 'drl', 'rt']),
            default='fr',
            show_default=True,
            help='Plotting layout that computes positions.',
        ),
        click.option(
            '--init-pos',
            type=click.STRING,
            default=None,
            show_default=True,
            help='Plotting layout that computes positions.',
        ),
        click.option(
            '--threshold',
            type=click.FLOAT,
            default=0.01,
            show_default=True,
            help='Do not draw edges for weights below this threshold. Set to 0 to '
            'include all edges.',
        ),
        COMMON_OPTIONS['root'],
        click.option(
            '--root',
            type=click.INT,
            default=0,
            show_default=True,
            help='If choosing a tree layout, this is the index of the root node.',
        ),
        click.option(
            '--transitions',
            type=click.STRING,
            default=None,
            show_default=True,
            help='Key for `.uns["paga"]` that specifies the matrix, e.g. '
            '`transition_confidence`, that stores the arrows.',
        ),
        click.option(
            '--single-component',
            is_flag=True,
            default=False,
            show_default=True,
            help='Restrict to largest connected component',
        ),
        click.option(
            '--solid-edges',
            type=click.Choice(['connectivities', 'connectivities_tree']),
            default='connectivities',
            show_default=True,
            help='Key for `.uns["paga"]` that specifies the matrix that stores the '
            'edges to be drawn solid black.',
        ),
        click.option(
            '--basis',
            type=click.STRING,
            default=None,
            show_default=True,
            help='Name of the embedding to plot, must be a key of `.obsm` without '
            'the prefix "X_".',
        ),
        click.option(
            '--color',
            type=CommaSeparatedText(simplify=True),
            default=None,
            show_default=True,
            help='Key(s) for annotation of observations/cells or variables/genes. Comma-separated if more than one',
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
            '--size',
            type=click.FLOAT,
            default=None,
            show_default=True,
            help='Point size. Automatically computed if not specified.',
        ),
        click.option(
            '--node-size-scale',
            type=click.FLOAT,
            default=1.0,
            show_default=True,
            help='Increase of decrease the size of the nodes.',
        ),
        click.option(
            '--fontsize',
            type=click.INT,
            default=None,
            show_default=True,
            help='Font size for node labels.',
        ),
        click.option(
            '--edge-width-scale',
            type=click.FLOAT,
            default=1.0,
            show_default=True,
            help='Increase of decrease the width of the edges.',
        ),
        click.option(
            '--arrowsize',
            type=click.INT,
            default=30,
            show_default=True,
            help='For directed graphs, specify the length and width of the arrowhead.',
        ),
        *COMMON_OPTIONS['opt_output'],
    ],

    'sviol': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['plot'],
        COMMON_OPTIONS['use_raw'],
        COMMON_OPTIONS['var_names'],
        *COMMON_OPTIONS['rank_genes_groups_plots'],
        COMMON_OPTIONS['layer'],
        *COMMON_OPTIONS['diffexp_plot'],
        COMMON_OPTIONS['gene_symbols'],
        *COMMON_OPTIONS['sviol'],
        COMMON_OPTIONS['swap_axes'],
    ],

    'dot': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['plot'],
        COMMON_OPTIONS['use_raw'],
        COMMON_OPTIONS['var_names'],
        *COMMON_OPTIONS['rank_genes_groups_plots'],
        COMMON_OPTIONS['layer'],
        *COMMON_OPTIONS['diffexp_plot'],
        COMMON_OPTIONS['gene_symbols'],
        *COMMON_OPTIONS['dot'],
    ],

    'matrix': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['plot'],
        COMMON_OPTIONS['use_raw'],
        COMMON_OPTIONS['var_names'],
        *COMMON_OPTIONS['rank_genes_groups_plots'],
        COMMON_OPTIONS['layer'],
        *COMMON_OPTIONS['diffexp_plot'],
        COMMON_OPTIONS['gene_symbols'],
    ],

    'heat': [
        *COMMON_OPTIONS['input'],
        *COMMON_OPTIONS['plot'],
        COMMON_OPTIONS['use_raw'],
        COMMON_OPTIONS['var_names'],
        *COMMON_OPTIONS['rank_genes_groups_plots'],
        COMMON_OPTIONS['layer'],
        *COMMON_OPTIONS['diffexp_plot'],
        COMMON_OPTIONS['gene_symbols'],
        *COMMON_OPTIONS['heat'],
        COMMON_OPTIONS['swap_axes'],
    ],

}
