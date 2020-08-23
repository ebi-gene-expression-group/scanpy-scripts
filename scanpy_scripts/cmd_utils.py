"""
Provide helper functions for constructing sub-commands
"""

import click
import pandas as pd
import scanpy as sc
from .cmd_options import CMD_OPTIONS
from .lib._paga import plot_paga

def make_subcmd(cmd_name, func, cmd_desc, arg_desc, opt_set = None):
    """
    Factory function that returns a sub-command function
    """
    opt_set = opt_set if opt_set else cmd_name
    options = CMD_OPTIONS[opt_set]
    option_spec = [click.command(cmd_name)]
    option_spec.extend(options)

    def add_docstring(cmd_desc, arg_desc):
        def docstring_dec(obj):
            obj.__doc__ = obj.__doc__.format(
                cmd_desc=cmd_desc, arg_desc=arg_desc)
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
            show_obj=None,
            **kwargs
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


def _read_obj(input_obj, input_format='anndata', **kwargs):
    if input_format == 'anndata':
        adata = sc.read(input_obj, **kwargs)
    elif input_format == 'loom':
        adata = sc.read_loom(input_obj, **kwargs)
    else:
        raise NotImplementedError(
            'Unsupported input format: {}'.format(input_format))
    return adata


def _write_obj(
        adata,
        output_obj,
        output_format='anndata',
        chunk_size=None,
        export_mtx=None,
        show_obj=None,
        write_obsm_varm=False,
        **kwargs
):
    if output_format == 'anndata':
        adata.write(output_obj, compression='gzip')
    elif output_format == 'loom':
        adata.write_loom(output_obj, write_obsm_varm=write_obsm_varm )
    elif output_format == 'zarr':
        adata.write_zarr(output_obj, chunk_size=chunk_size, **kwargs)
    else:
        raise NotImplementedError(
            'Unsupported output format: {}'.format(output_format))
    if export_mtx:
        write_mtx(adata, fname_prefix=export_mtx, **kwargs)
    if show_obj:
        click.echo(adata, err=show_obj == 'stderr')
    return 0


def write_mtx(adata, fname_prefix='', var=None, obs=None, use_raw=False):
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
    """
    if fname_prefix and not (fname_prefix.endswith('/') or fname_prefix.endswith('_')):
        fname_prefix = fname_prefix + '_'
    if var is None:
        var = []
    if obs is None:
        obs = []
    if use_raw:
        adata = adata.raw
    obs = list(set(obs) & set(adata.obs.columns))
    var = list(set(var) & set(adata.var.columns))

    import scipy.sparse as sp
    mat = sp.coo_matrix(adata.X)
    n_obs, n_var = mat.shape
    n_entry = len(mat.data)
    header = '%%MatrixMarket matrix coordinate real general\n%\n{} {} {}\n'.format(
        n_var, n_obs, n_entry)
    df = pd.DataFrame({'col': mat.col + 1, 'row': mat.row + 1, 'data': mat.data})
    mtx_fname = fname_prefix + 'matrix.mtx'
    gene_fname = fname_prefix + 'genes.tsv'
    barcode_fname = fname_prefix + 'barcodes.tsv'
    with open(mtx_fname, 'a') as fh:
        fh.write(header)
        df.to_csv(fh, sep=' ', header=False, index=False)

    obs_df = adata.obs[obs].reset_index(level=0)
    obs_df.to_csv(barcode_fname, sep='\t', header=False, index=False)
    var_df = adata.var[var].reset_index(level=0)
    if not var:
        var_df['gene'] = var_df['index']
    var_df.to_csv(gene_fname, sep='\t', header=False, index=False)


def make_plot_function(func_name, kind=None):
    """Make plot function that handles common plotting parameters
    """

    # Provide a function translation

    plot_funcs = {
        'embedding': sc.pl.embedding,
        'scatter': sc.pl.scatter,
        'sviol': sc.pl.stacked_violin,
        'rgg_sviol': sc.pl.rank_genes_groups_stacked_violin,
        'dot': sc.pl.dotplot,
        'rgg_dot': sc.pl.rank_genes_groups_dotplot,
        'matrix': sc.pl.matrixplot,
        'rgg_matrix': sc.pl.rank_genes_groups_matrixplot,
        'heat': sc.pl.heatmap,
        'rgg_heat': sc.pl.rank_genes_groups_heatmap,
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
            rcParams.update({'figure.figsize': fig_size})

        # Choose the function to run

        is_rgg=False

        if func_name in plot_funcs:
            if 'rgg' in kwargs:
                if kwargs['rgg'] == True:
                    is_rgg=True
                    func=plot_funcs[ 'rgg_' + func_name ]
                    kwargs.pop('var_names', None)
                else:
                    func = plot_funcs[ func_name ]
                    kwargs.pop('groups', None)
                    kwargs.pop('n_genes', None)
        
                kwargs.pop('rgg')
            else:
                func = plot_funcs[ func_name ]
        else:
            func = globals()[func_name]

        # Generate the output file name

        figname = False
        showfig = True
        if output_fig:
            import os
            import matplotlib.pyplot as plt
            sc.settings.figdir = os.path.dirname(output_fig) or '.'

            figname = os.path.basename(output_fig)
            showfig = False

        # Run the selected function

        func(
            adata,
            save=figname,
            show=showfig,
            **kwargs)

        # Rename output to the spefied file name. We need to work out what
        # prefix the function will have used for its output files.

        if output_fig:
            prefix=''
            if func_name == 'scatter' or func_name == 'embedding':
                prefix =  kwargs.get('basis', func.__name__)
            elif kind:
                prefix = kind 
            elif func_name in plot_funcs:
                prefix = plot_funcs[ func_name ].__name__.split('.')[-1]
                if func_name in [ 'sviol', 'rgg_sviol', 'dot', 'rgg_dot', 'matrix', 'rgg_matrix' ]:
                    prefix = prefix + '_'

            os.rename(
                os.path.join(sc.settings.figdir, prefix + figname), output_fig)
            plt.close()

    return plot_function
