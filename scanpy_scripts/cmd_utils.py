"""
Provide helper functions for constructing sub-commands
"""

import click
import pandas as pd
import scanpy as sc
import scanpy.external as sce

from .cmd_options import CMD_OPTIONS
from .lib._paga import plot_paga
from .lib._scrublet import plot_scrublet
from .obj_utils import _save_matrix


def make_subcmd(cmd_name, func, cmd_desc, arg_desc, opt_set=None):
    """
    Factory function that returns a sub-command function
    """
    opt_set = opt_set if opt_set else cmd_name
    options = CMD_OPTIONS[opt_set]
    option_spec = [click.command(cmd_name)]
    option_spec.extend(options)

    def add_docstring(cmd_desc, arg_desc):
        def docstring_dec(obj):
            obj.__doc__ = obj.__doc__.format(cmd_desc=cmd_desc, arg_desc=arg_desc)
            return obj

        return docstring_dec

    @add_options(option_spec)
    @add_docstring(cmd_desc, arg_desc)
    def cmd(
        input_obj=None,
        output_obj=None,
        input_format=None,
        output_format=None,
        zarr_chunk_size=None,
        loom_write_obsm_varm=False,
        export_mtx=None,
        mtx_compression=None,
        show_obj=None,
        **kwargs,
    ):
        """{cmd_desc}\n\n\b\n{arg_desc}"""
        if input_obj:
            adata = _read_obj(input_obj, input_format=input_format)
            func(adata, **kwargs)
        else:
            adata = func(**kwargs)

        if output_obj:
            _write_obj(
                adata,
                output_obj,
                output_format=output_format,
                chunk_size=zarr_chunk_size,
                write_obsm_varm=loom_write_obsm_varm,
                export_mtx=export_mtx,
                mtx_compression=mtx_compression,
                show_obj=show_obj,
            )
        return 0

    return cmd


def add_options(options):
    """
    Returns a decorator to group multiple click decorators
    """

    def _add_options(func):
        for option in reversed(options):
            func = option(func)
        return func

    return _add_options


def _fix_booleans(df):
    for var in df.columns:
        if (
            df[var].dtype.kind == "O"
            and df[var].dtype.name == "object"
            and set(pd.Categorical(df[var])).issubset(set(["True", "False", "nan"]))
        ):
            d = {"False": True, "False": False, "nan": False}
            df[var] = df[var].map(d).astype(bool)
    return df


def _read_obj(input_obj, input_format="anndata", **kwargs):
    if input_format == "anndata":
        adata = sc.read_h5ad(input_obj, **kwargs)
    elif input_format == "loom":
        adata = sc.read_loom(input_obj, **kwargs)
    else:
        raise NotImplementedError("Unsupported input format: {}".format(input_format))
    adata.var = _fix_booleans(adata.var)
    adata.obs = _fix_booleans(adata.obs)

    return adata


def _write_obj(
    adata,
    output_obj,
    output_format="anndata",
    chunk_size=None,
    export_mtx=None,
    mtx_compression=None,
    show_obj=None,
    write_obsm_varm=False,
    **kwargs,
):
    if output_format == "anndata":
        adata.write(output_obj, compression="gzip")
    elif output_format == "loom":
        adata.write_loom(output_obj, write_obsm_varm=write_obsm_varm)
    elif output_format == "zarr":
        adata.write_zarr(output_obj, chunk_size=chunk_size, **kwargs)
    else:
        raise NotImplementedError("Unsupported output format: {}".format(output_format))
    if export_mtx:
        compression = None
        if mtx_compression is not None:
            compression = {"method": mtx_compression}

        write_mtx(adata, fname_prefix=export_mtx, compression=compression, **kwargs)
    if show_obj:
        click.echo(adata, err=show_obj == "stderr")
    return 0


def write_mtx(
    adata,
    fname_prefix="",
    var=None,
    obs=None,
    use_raw=False,
    use_layer=None,
    compression=None,
):
    """Export AnnData object to mtx formt
    * Parameters
        + adata : AnnData
        An AnnData object
        + fname_prefix : str
        Prefix of the exported files. If not empty and not ending with '/' or '_',
        a '_' will be appended. Full names will be <fname_prefix>matrix.mtx,
        <fname_prefix>genes.tsv, <fname_prefix>barcodes.tsv
        + var : list
        A list of column names to be exported to gene table
        + obs : list
        A list of column names to be exported to barcode/cell table
        + use_raw : bool
        Take data the matrix from .raw.X?
        + use_layer: str
        Specify a layer to use instead of .X (non-raw only)
        + compression: None, str or dict
        Compression parameter for Pandas' to_csv(). For compression, a dict
        with a 'method' key, e.g. {'method': 'gzip', 'compresslevel': 1,
        'mtime': 1}

    >>> import os
    >>> from pathlib import Path
    >>> adata = sc.datasets.pbmc3k()
    >>> # Test uncompressed write
    >>> Path("uncompressed").mkdir(parents=True, exist_ok=True)
    >>> write_mtx(adata, fname_prefix = 'uncompressed/', use_raw = False, use_layer = None, var = ['gene_name'])
    >>> sorted(os.listdir('uncompressed'))
    ['barcodes.tsv', 'genes.tsv', 'matrix.mtx']
    >>> # Test that the matrix is the same when we read it back
    >>> test_readable = sc.read_10x_mtx('uncompressed')
    >>> if any(test_readable.obs_names != adata.obs_names) or any(test_readable.var_names != adata.var_names) or (test_readable.X[1].sum() - adata.X[1].sum()) > 1e-5:
    ...   print("Re-read matrix is different to the one we stored, something is wrong with the writing")
    >>> # Test compressed write
    >>> Path("compressed").mkdir(parents=True, exist_ok=True)
    >>> write_mtx(adata, fname_prefix = 'compressed/', use_raw = False, use_layer = None, var = ['gene_name'], compression = {'method': 'gzip'})
    >>> sorted(os.listdir('compressed'))
    ['barcodes.tsv.gz', 'genes.tsv.gz', 'matrix.mtx.gz']
    """
    if fname_prefix and not (fname_prefix.endswith("/") or fname_prefix.endswith("_")):
        fname_prefix = fname_prefix + "_"
    if var is None:
        var = []
    if obs is None:
        obs = []

    import scipy.sparse as sp

    if use_raw:
        var_source = adata.raw.var
        mat = sp.coo_matrix(adata.raw.X)
    else:
        var_source = adata.var
        if use_layer is not None:
            mat = sp.coo_matrix(adata.layers[use_layer])
        else:
            mat = sp.coo_matrix(adata.X)

    obs = list(set(obs) & set(adata.obs.columns))
    var = list(set(var) & set(var_source.columns))

    n_obs, n_var = mat.shape
    n_entry = len(mat.data)

    # Define the header lines as a Pandas DataFrame so we can use the same compression
    header = pd.DataFrame(
        ["%%MatrixMarket matrix coordinate real general", f"{n_var} {n_obs} {n_entry}"]
    )
    df = pd.DataFrame({"col": mat.col + 1, "row": mat.row + 1, "data": mat.data})

    # Define outputs
    mtx_fname = fname_prefix + "matrix.mtx"
    gene_fname = fname_prefix + "genes.tsv"
    barcode_fname = fname_prefix + "barcodes.tsv"

    # Write matrix with Pandas CSV and use its compression where requested
    if (
        compression is not None
        and type(compression) is dict
        and "method" in compression
    ):
        compressed_exts = {"zip": "zip", "gzip": "gz", "bz2": "bz2", "zstd": "zst"}
        ext = compressed_exts.get(compression["method"], "None")

        if ext is None:
            errmsg = "Invalid compression method"
            raise Exception(errmsg)

        mtx_fname += f".{ext}"
        gene_fname += f".{ext}"
        barcode_fname += f".{ext}"
    else:
        compression = None

    header.to_csv(mtx_fname, header=False, index=False, compression=compression)
    df.to_csv(
        mtx_fname, sep=" ", header=False, index=False, compression=compression, mode="a"
    )

    # Now write the obs and var, also with compression if appropriate
    obs_df = adata.obs[obs].reset_index(level=0)
    obs_df.to_csv(
        barcode_fname, sep="\t", header=False, index=False, compression=compression
    )
    var_df = var_source[var].reset_index(level=0)
    if not var:
        var_df["gene"] = var_df["index"]
    var_df.to_csv(
        gene_fname, sep="\t", header=False, index=False, compression=compression
    )


def make_plot_function(func_name, kind=None):
    """Make plot function that handles common plotting parameters"""

    # Provide a function translation

    plot_funcs = {
        "embedding": sc.pl.embedding,
        "scatter": sc.pl.scatter,
        "sviol": sc.pl.stacked_violin,
        "rgg_sviol": sc.pl.rank_genes_groups_stacked_violin,
        "dot": sc.pl.dotplot,
        "rgg_dot": sc.pl.rank_genes_groups_dotplot,
        "matrix": sc.pl.matrixplot,
        "rgg_matrix": sc.pl.rank_genes_groups_matrixplot,
        "heat": sc.pl.heatmap,
        "rgg_heat": sc.pl.rank_genes_groups_heatmap,
    }

    def plot_function(
        adata,
        output_fig=None,
        fig_size=None,
        fig_dpi=300,
        fig_fontsize=15,
        **kwargs,
    ):
        sc.settings.set_figure_params(dpi=fig_dpi, fontsize=fig_fontsize)
        if fig_size:
            from matplotlib import rcParams

            rcParams.update({"figure.figsize": fig_size})

        # Choose the function to run

        is_rgg = False

        if func_name in plot_funcs:
            if "rgg" in kwargs:
                if kwargs["rgg"] == True:
                    is_rgg = True
                    func = plot_funcs["rgg_" + func_name]
                    kwargs.pop("var_names", None)
                else:
                    func = plot_funcs[func_name]
                    kwargs.pop("groups", None)
                    kwargs.pop("n_genes", None)

                kwargs.pop("rgg")
            else:
                func = plot_funcs[func_name]
        else:
            func = globals()[func_name]

        # Generate the output file name

        figname = False
        showfig = True
        if output_fig:
            import os

            import matplotlib.pyplot as plt

            sc.settings.figdir = os.path.dirname(output_fig) or "."

            figname = os.path.basename(output_fig)
            showfig = False

        # Run the selected function

        func(adata, save=figname, show=showfig, **kwargs)

        # Rename output to the spefied file name. We need to work out what
        # prefix the function will have used for its output files.

        if output_fig:
            prefix = ""
            if func_name == "scatter" or func_name == "embedding":
                prefix = kwargs.get("basis", func.__name__)
            elif kind:
                prefix = kind
            elif func_name in plot_funcs:
                prefix = plot_funcs[func_name].__name__.split(".")[-1]
                if func_name in [
                    "sviol",
                    "rgg_sviol",
                    "dot",
                    "rgg_dot",
                    "matrix",
                    "rgg_matrix",
                ]:
                    prefix = prefix + "_"

            os.rename(os.path.join(sc.settings.figdir, prefix + figname), output_fig)
            plt.close()

    return plot_function


# Wrap matrix-processing functions in logic to back up .X or specified input
# layers prior to processing


def make_matrix_function(func):
    def matrix_function(
        adata,
        save_raw=True,
        save_layer=None,
        **kwargs,
    ):

        # For the subset of matrix functions that allow layer specification,
        # pass that as the thing to save.

        layer = None
        if "layer" in kwargs:
            layer = kwargs["layer"]

        _save_matrix(adata, save_raw, save_layer=save_layer, layer=layer)
        func(adata, **kwargs)
        return adata

    return matrix_function
