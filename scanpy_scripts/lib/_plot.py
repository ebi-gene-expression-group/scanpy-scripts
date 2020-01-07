"""
Plotting related functions
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LinearSegmentedColormap
import scanpy as sc
from scanpy.plotting._tools.scatterplots import plot_scatter

from ._diffexp import extract_de_table


def expression_colormap(background_level=0.01):
    """Returns a nice color map for highlighting gene expression
    """
    background_nbin = int(100 * background_level)
    reds = plt.cm.Reds(np.linspace(0, 1, 100-background_nbin))
    greys = plt.cm.Greys_r(np.linspace(0.7, 0.8, background_nbin))
    palette = np.vstack([greys, reds])
    return LinearSegmentedColormap.from_list('expression', palette)


def _is_numeric(x):
    return x.dtype.kind in ('i', 'f')


def cross_table(adata, x, y, normalise=None, highlight=None, subset=None):
    """Make a cross table comparing two categorical annotations
    """
    x_attr = adata.obs[x]
    y_attr = adata.obs[y]
    if subset is not None:
        x_attr = x_attr[subset]
        y_attr = y_attr[subset]
    assert not _is_numeric(x_attr.values), f'Can not operate on numerical {x}'
    assert not _is_numeric(y_attr.values), f'Can not operate on numerical {y}'
    crs_tbl = pd.crosstab(x_attr, y_attr)
    if normalise == 'x':
        x_sizes = x_attr.groupby(x_attr).size().values
        crs_tbl = (crs_tbl.T / x_sizes * 100).round(2).T
    elif normalise == 'y':
        y_sizes = x_attr.groupby(y_attr).size().values
        crs_tbl = (crs_tbl / y_sizes * 100).round(2)
    elif normalise == 'xy':
        x_sizes = x_attr.groupby(x_attr).size().values
        crs_tbl = (crs_tbl.T / x_sizes * 100).T
        y_sizes = crs_tbl.sum(axis=0)
        crs_tbl = (crs_tbl / y_sizes * 100).round(2)
    elif normalise == 'yx':
        y_sizes = x_attr.groupby(y_attr).size().values
        crs_tbl = (crs_tbl / y_sizes * 100).round(2)
        x_sizes = crs_tbl.sum(axis=1)
        crs_tbl = (crs_tbl.T / x_sizes * 100).round(2).T
    if highlight is not None:
        return crs_tbl.style.background_gradient(cmap='viridis', axis=highlight)
    return crs_tbl


def set_figsize(dim):
    if len(dim) == 2 and (isinstance(dim[0], (int, float))
            and isinstance(dim[1], (int, float))):
        rcParams.update({'figure.figsize': dim})
    else:
        raise ValueError(f'Invalid {dim} value, must be an iterable of '
                          'length two in the form of (width, height).')


def plot_df_heatmap(
        df,
        cmap='viridis',
        title=None,
        figsize=(7, 7),
        rotation=90,
        save=None,
        **kwargs,
):
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(df, cmap=cmap, aspect='auto', **kwargs)
    if 0 < rotation < 90:
        horizontalalignment = 'right'
    else:
        horizontalalignment = 'center'
    plt.xticks(
        range(len(df.columns)),
        df.columns,
        rotation=rotation,
        horizontalalignment=horizontalalignment,
    )
    plt.yticks(range(len(df.index)), df.index)
    if title:
        fig.suptitle(title)
    fig.colorbar(im)
    if save:
        plt.savefig(fname=save, bbox_inches='tight', pad_inches=0.1)


def plot_qc(adata, groupby=None, groupby_only=False):
    qc_metrics = ('n_counts', 'n_genes', 'percent_mito', 'percent_ribo', 'percent_hb')
    for qmt in qc_metrics:
        if qmt not in adata.obs.columns:
            raise ValueError(f'{qmt} not found.')
    if groupby:
        ax = sc.pl.violin(adata, keys=qc_metrics, groupby=groupby, rotation=45, show=False)
    if not groupby_only:
        ax = sc.pl.violin(adata, keys=qc_metrics, multi_panel=True)
        old_figsize = rcParams.get('figure.figsize')
        rcParams.update({'figure.figsize': (4,3)})
        sc.pl.scatter(adata, x='n_counts', y='n_genes', color=groupby, alpha=0.5)
        sc.pl.scatter(adata, x='n_counts', y='n_genes', color='percent_mito', alpha=0.5)
        sc.pl.scatter(adata, x='n_counts', y='percent_mito', color=groupby, alpha=0.5)
        sc.pl.scatter(adata, x='n_counts', y='n_genes', color='percent_ribo', alpha=0.5)
        sc.pl.scatter(adata, x='n_counts', y='percent_ribo', color=groupby, alpha=0.5)
        sc.pl.scatter(adata, x='n_counts', y='n_genes', color='percent_hb', alpha=0.5)
        sc.pl.scatter(adata, x='n_counts', y='percent_hb', color=groupby, alpha=0.5)
        rcParams.update({'figure.figsize': old_figsize})


def plot_metric_by_rank(
        adata,
        subject='cell',
        metric='n_counts',
        kind='rank',
        nbins=50,
        decreasing=True,
        ax=None,
        title=None,
        hpos=None,
        vpos=None,
        logx=True,
        logy=True,
        swap_axis=False,
        **kwargs,
):
    """Plot metric by rank from top to bottom"""
    kwargs['c'] = kwargs.get('c', 'black')
    metric_names = {
        'n_counts': 'nUMI', 'n_genes': 'nGene', 'n_cells': 'nCell', 'mt_prop': 'fracMT'}
    if subject not in ('cell', 'gene'):
        print('`subject` must be "cell" or "gene".')
        return
    if metric not in metric_names:
        print('`metric` must be "n_counts", "n_genes", "n_cells" or "mt_prop".')
        return

    if 'mt_prop' not in adata.obs.columns:
        k_mt = adata.var_names.str.startswith('MT-')
        if sum(k_mt) > 0:
            adata.obs['mt_prop'] = np.squeeze(
                np.asarray(adata.X[:, k_mt].sum(axis=1))) / adata.obs['n_counts']

    if not ax:
        fig, ax = plt.subplots()

    if kind == 'rank':
        order_modifier = -1 if decreasing else 1

        if subject == 'cell':
            if swap_axis:
                x = 1 + np.arange(adata.shape[0])
                k = np.argsort(order_modifier * adata.obs[metric])
                y = adata.obs[metric].values[k]
            else:
                y = 1 + np.arange(adata.shape[0])
                k = np.argsort(order_modifier * adata.obs[metric])
                x = adata.obs[metric].values[k]
        else:
            if swap_axis:
                x = 1 + np.arange(adata.shape[1])
                k = np.argsort(order_modifier * adata.var[metric])
                y = adata.var[metric].values[k]
            else:
                y = 1 + np.arange(adata.shape[1])
                k = np.argsort(order_modifier * adata.var[metric])
                x = adata.var[metric].values[k]

        if kwargs['c'] is not None and not isinstance(kwargs['c'], str):
            kwargs_c = kwargs['c'][k]
            del kwargs['c']
            ax.scatter(x, y, c=kwargs_c, **kwargs)
        else:
            ax.plot(x, y, **kwargs)

        if logy:
            ax.set_yscale('log')
        if swap_axis:
            ax.set_xlabel('Rank')
        else:
            ax.set_ylabel('Rank')

    elif kind == 'hist':
        if subject == 'cell':
            value = adata.obs[metric]
            logbins = np.logspace(np.log10(np.min(value[value>0])), np.log10(np.max(adata.obs[metric])), nbins)
            h = ax.hist(value[value>0], logbins)
            y = h[0]
            x = h[1]
        else:
            value = adata.var[metric]
            logbins = np.logspace(np.log10(np.min(value[value>0])), np.log10(np.max(adata.var[metric])), nbins)
            h = ax.hist(value[value>0], logbins)
            y = h[0]
            x = h[1]

        ax.set_ylabel('Count')

    if hpos is not None:
        ax.hlines(hpos, xmin=x.min(), xmax=x.max(), linewidth=1, colors='red')
    if vpos is not None:
        ax.vlines(vpos, ymin=y.min(), ymax=y.max(), linewidth=1, colors='green')
    if logx:
        ax.set_xscale('log')
    if swap_axis:
        ax.set_ylabel(metric_names[metric])
    else:
        ax.set_xlabel(metric_names[metric])
    if title:
        ax.set_title(title)
    ax.grid(linestyle='--')
    if 'label' in kwargs:
        ax.legend(loc='upper right' if decreasing else 'upper left')

    return ax


def plot_embedding(
        adata, basis, groupby, color=None, annot=True, highlight=None, size=None,
        save=None, savedpi=300, **kwargs):
    if f'X_{basis}' not in adata.obsm.keys():
        raise KeyError(f'"X_{basis}" not found in `adata.obsm`.')
    if isinstance(groupby, (list, tuple)):
        groupby = groupby[0]
    if groupby not in adata.obs.columns:
        raise KeyError(f'"{groupby}" not found in `adata.obs`.')
    if adata.obs[groupby].dtype.name != 'category':
        raise ValueError(f'"{groupby}" is not categorical.')
    groups = adata.obs[groupby].copy()
    categories = list(adata.obs[groupby].cat.categories)
    rename_dict = {ct: f'{i}: {ct}' for i, ct in enumerate(categories)}
    restore_dict = {f'{i}: {ct}': ct for i, ct in enumerate(categories)}

    size_ratio = 1.2

    ad = adata
    marker_size = size
    kwargs['show'] = False
    kwargs['save'] = False
    kwargs['frameon'] = True
    kwargs['legend_loc'] = 'right margin'

    color = groupby if color is None else color
    offset = 0 if 'diffmap' in basis else -1
    xi, yi = 1, 2
    if 'components' in kwargs:
        xi, yi = components
    xi += offset
    yi += offset

    if highlight:
        k = adata.obs[groupby].isin(highlight).values
        ad = adata[~k, :]
        marker_size = size / size_ratio if size else None
    if annot:
        adata.obs[groupby].cat.rename_categories(rename_dict, inplace=True)
        if highlight:
            ad.obs[groupby].cat.rename_categories(rename_dict, inplace=True)
    else:
        kwargs['frameon'] = False
        kwargs['title'] = ''
        kwargs['legend_loc'] = None

    fig, ax = plt.subplots()
    try:
        plot_scatter(ad, basis, color=color, ax=ax, size=marker_size, **kwargs)

        if highlight:
            embed = adata.obsm[f'X_{basis}']
            for i, ct in enumerate(categories):
                if ct not in highlight:
                    continue
                k_hl = groups == ct
                ax.scatter(embed[k_hl, xi], embed[k_hl, yi], marker='D', c=adata.uns[f'{groupby}_colors'][i], s=marker_size)
                ax.scatter([], [], marker='D', c=adata.uns[f'{groupby}_colors'][i], label=adata.obs[groupby].cat.categories[i])
            if annot:
                ax.legend(frameon=False, loc='center left', bbox_to_anchor=(1, 0.5),
                          ncol=(1 if len(categories) <= 14 else 2 if len(categories) <= 30 else 3))
    finally:
        if annot:
            adata.obs[groupby].cat.rename_categories(restore_dict, inplace=True)
    if annot:
        centroids = pseudo_bulk(adata, groupby, use_rep=f'X_{basis}', FUN=np.median).T
        texts = [ax.text(x=row[xi], y=row[yi], s=f'{i:d}', fontsize=8, fontweight='bold') for i, row in centroids.reset_index(drop=True).iterrows()]
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, text_from_points=False, autoalign=False)
    if save:
        plt.savefig(fname=save, dpi=savedpi, bbox_inches='tight', pad_inches=0.1)


def plot_diffexp(
    adata,
    basis='umap',
    key='rank_genes_groups',
    top_n=4,
    extra_genes=None,
    figsize1=(4,4),
    figsize2=(2.5, 2.5),
    dotsize=None,
    **kwargs
):
    grouping = adata.uns[key]['params']['groupby']
    de_tbl = extract_de_table(adata.uns[key])
    de_tbl = de_tbl.loc[de_tbl.genes.astype(str) != 'nan', :]
    de_genes = list(de_tbl.groupby('cluster').head(top_n)['genes'].values)
    de_clusters = list(de_tbl.groupby('cluster').head(top_n)['cluster'].astype(str).values)
    if extra_genes:
        de_genes.extend(extra_genes)
        de_clusters.extend(['known'] * len(extra_genes))

    rcParams.update({'figure.figsize': figsize1})
    sc.pl.rank_genes_groups(adata, key=key, show=False)

    sc.pl.dotplot(adata, var_names=de_genes, groupby=grouping, show=False)

    expr_cmap = expression_colormap(0.01)
    rcParams.update({'figure.figsize':figsize2})
    plot_scatter(
        adata,
        basis=basis,
        color=de_genes,
        color_map=expr_cmap,
        use_raw=True,
        size=dotsize,
        title=[f'{c}, {g}' for c, g in zip(de_clusters, de_genes)],
        show=False,
        **kwargs,
    )

    rcParams.update({'figure.figsize':figsize1})
    plot_embedding(
        adata, basis=basis, groupby=grouping, size=dotsize, show=False)

    return de_tbl
