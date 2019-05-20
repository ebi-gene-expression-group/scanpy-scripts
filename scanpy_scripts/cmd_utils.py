"""
Provide helper functions for constructing sub-commands
"""

import click
import pandas as pd
import scanpy as sc


def make_subcmd(cmd_name, options, func, cmd_desc, arg_desc):
    """
    Factory function that returns a sub-command function
    """
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

        _write_obj(
            adata,
            output_obj,
            output_format=output_format,
            chunk_size=zarr_chunk_size,
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


def _read_obj(input_obj, input_format='anndata'):
    if input_format == 'anndata':
        adata = sc.read(input_obj)
    elif input_format == 'loom':
        adata = sc.read_loom(input_obj)
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
        **kwargs
):
    if output_format == 'anndata':
        adata.write(output_obj, compression='gzip')
    elif output_format == 'loom':
        adata.write_loom(output_obj, **kwargs)
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


def write_cluster(adata, names, cluster_fn, sep='\t'):
    """Export cell clustering as a text table
    """
    if not isinstance(names, (list, tuple)):
        names = [names]
    for name in names:
        if name not in adata.obs.keys():
            raise KeyError(f'{name} is not a valid `.uns` key')
    adata[names].to_csv(cluster_fn, sep=sep, header=True, index=True)


def write_embedding(adata, name, embed_fn, n_comp=None, sep='\t'):
    """Export cell embeddings as a txt table
    """
    if name not in adata.obsm.keys():
        raise KeyError(f'{name} is not a valid `.obsm` key')
    mat = adata.obsm[name].copy()
    if n_comp is not None and mat.shape[1] >= n_comp:
        mat = mat[:, 0:n_comp]
    pd.DataFrame(mat, index=adata.obs_names).to_csv(
        embed_fn, sep=sep, header=False, index=True)


def plot_embeddings(
        adata,
        basis,
        save_fig=None,
        fig_size=None,
        fig_dpi=300,
        fig_fontsize=15,
        **kwargs,
):
    """Make scatter plot of cell embeddings
    """
    name = 'X_' + basis
    if name not in adata.obsm.keys():
        raise KeyError(f'{name} is not a valid `.obsm` key')

    sc.settings.set_figure_params(dpi=fig_dpi, fontsize=fig_fontsize)
    if fig_size:
        from matplotlib import rcParams
        rcParams.update({'figure.figsize': fig_size})
    if save_fig:
        import os
        import matplotlib.pyplot as plt
        sc.settings.figdir = os.path.dirname(save_fig)

        figname = os.path.basename(save_fig)
        sc.pl.scatter(adata, basis, save=figname, show=False, **kwargs)
        os.rename(os.path.join(sc.settings.figdir, 'scatter' + figname), save_fig)
        plt.close()
    else:
        sc.pl.scatter(adata, basis, save=False, show=True, **kwargs)
