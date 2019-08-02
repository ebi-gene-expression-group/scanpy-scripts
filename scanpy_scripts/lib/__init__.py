"""
Provides exported functions
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scanpy as sc

from ._read import read_10x
from ._filter import filter_anndata
from ._norm import normalize
from ._hvg import hvg
from ._pca import pca
from ._neighbors import neighbors
from ._umap import umap
from ._fdg import fdg
from ._tsne import tsne
from ._louvain import louvain
from ._leiden import leiden
from ._diffexp import diffexp, diffexp_paired, extract_de_table
from ._diffmap import diffmap
from ._dpt import dpt
from ._paga import paga
from ..cmd_utils import _read_obj as read_obj
from ..cmd_utils import _write_obj as write_obj
from ..cmd_utils import switch_layer


def expression_colormap():
    """Returns a nice color map for highlighting gene expression
    """
    reds = plt.cm.Reds(np.linspace(0, 1, 118))
    greys = plt.cm.Greys_r(np.linspace(0.7, 0.8, 10))
    palette = np.vstack([greys, reds])
    return LinearSegmentedColormap.from_list('expression', palette)


def cross_table(adata, x, y, normalise=None, highlight=False, subset=None):
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
    if highlight:
        return crs_tbl.style.background_gradient(cmap='viridis', axis=0)
    return crs_tbl


def _is_numeric(x):
    return x.dtype.kind in ('i', 'f')


def run_harmony(
        adata,
        batch,
        theta=2.0,
        use_rep='X_pca',
        key_added='hm',
):
    if not isinstance(batch, (tuple, list)):
        batch = [batch]
    if not isinstance(theta, (tuple, list)):
        theta = [theta]
    for b in batch:
        if b not in adata.obs.columns:
            raise KeyError(f'{b} is not a valid obs annotation.')
    if use_rep not in adata.obsm.keys():
        raise KeyError(f'{use_rep} is not a valid embedding.')
    meta = adata.obs[batch].reset_index()
    embed = adata.obsm[use_rep]

    # ===========
    from rpy2.robjects.packages import importr
    harmony = importr('harmony')
    from rpy2.robjects import numpy2ri, pandas2ri
    numpy2ri.activate()
    pandas2ri.activate()
    hm_embed = harmony.HarmonyMatrix(
            embed, meta, batch, theta, do_pca=False, verbose=False)
    pandas2ri.deactivate()
    numpy2ri.deactivate()
    hm_embed = numpy2ri.ri2py(hm_embed)
    # ===========
    if key_added:
        adata.obsm[f'{use_rep}_{key_added}'] = hm_embed
    else:
        adata.obsm[use_rep] = hm_embed


def split_by_group(adata, groupby):
    if groupby not in adata.obs.columns:
        raise KeyError(f'{groupby} is not a valid obs annotation.')
    groups = adata.obs[groupby].unique().to_list()
    adata_dict = {}
    for grp in groups:
        adata_dict[grp] = adata[adata.obs[groupby] == grp, :]
    return adata_dict


def regroup(adata, groupby, regroups):
    if groupby not in adata.obs.columns:
        raise KeyError(f'{groupby} is not a valid obs annotation.')
    groups = adata.obs[groupby].astype(str)
    new_groups = groups.copy()
    for new_grp, old_grps in regroups.items():
        for grp in old_grps:
            new_groups[groups == grp] = new_grp
    return new_groups.astype('category')


def subsample(adata, fraction, groupby=None, min_n=0, max_n=10000, **kwargs):
    if groupby:
        if groupby not in adata.obs.columns:
            raise KeyError(f'{groupby} is not a valid obs annotation.')
        groups = adata.obs[groupby].unique()
        n_obs_per_group = {}
        for grp in groups:
            k = adata.obs[groupby] == grp
            grp_size = sum(k)
            n_obs_per_group[grp] = int(min(max_n,
                                           max(np.ceil(grp_size * fraction),
                                               min(min_n, grp_size))))
        sampled_groups = [sc.pp.subsample(adata[adata.obs[groupby] == grp, :],
                                          n_obs=n_obs_per_group[grp],
                                          copy=True,
                                          **kwargs) for grp in groups]
        sampled_obs_names = np.concatenate(
                [sg.obs_names.values for sg in sampled_groups])
        subsampled = adata[adata.obs_names.isin(sampled_obs_names), :]
    else:
        subsampled = sc.pp.subsample(adata, fraction, **kwargs, copy=True)
    return subsampled


def pseudo_bulk(
        adata, groupby, use_rep='X', highly_variable=False, FUN=np.mean):
    """Make pseudo bulk data from grouped sc data
    """
    if adata.obs[groupby].dtype.name == 'category':
        group_attr = adata.obs[groupby].values
        groups = adata.obs[groupby].cat.categories.values
    else:
        group_attr = adata.obs[groupby].astype(str).values
        groups = np.unique(group_attr)
    n_level = len(groups)
    if highly_variable:
        if isinstance(highly_variable, (list, tuple)):
            k_hv = adata.var_names.isin(highly_variable)
        else:
            k_hv = adata.var['highly_variable'].values
    if use_rep == 'X':
        x = adata.X
        features = adata.var_names.values
        if highly_variable:
            x = x[:, k_hv]
            features = features[k_hv]
    elif use_rep in adata.layers.keys():
        x = adata.layers[use_rep]
        features = adata.var_names.values
        if highly_variable:
            x = x[:, k_hv]
            features = features[k_hv]
    elif use_rep in adata.obsm.keys():
        x = adata.obsm[use_rep]
        features = np.arange(x.shape[1])
    elif (isinstance(use_rep, np.ndarray) and
            use_rep.shape[0] == adata.shape[0]):
        x = use_rep
        features = np.arange(x.shape[1])
    else:
        raise KeyError(f'{use_rep} invalid.')
    summarised = np.zeros((n_level, x.shape[1]))
    for i, grp in enumerate(groups):
        k_grp = group_attr == grp
        summarised[i] = FUN(x[k_grp, :], axis=0, keepdims=True)
    return pd.DataFrame(summarised.T, columns=groups, index=features)


def LR_annotate(
        adata,
        train_label,
        train_x,
        use_rep='X',
        subset=None,
        penalty='l1',
        sparcity=0.2,
        return_label=False,
):
    """Annotate cells using logistic regression
    """
    if subset is not None:
        train_label = train_label[subset]
        train_x = train_x[subset, :]
    from sklearn.linear_model import LogisticRegression
    lr = LogisticRegression(penalty=penalty, C=sparcity)
    lr.fit(train_x, train_label)
    if use_rep == 'X':
        predict_x = adata.X
    elif use_rep in adata.layers.keys():
        predict_x = adata.layers[use_rep]
    elif use_rep in adata.obsm.keys():
        predict_x = adata.obsm[use_rep]
    if return_label:
        return lr.predict(predict_x)
    else:
        return lr.predict_proba(predict_x)


def annotate(adata, groupby, annot, annotation_matrix):
    """Annotate a clustering based on a matrix of values

    + groupby: the clustering to be annotated
    + annot : the annotationo to be aligned/transferred
    + annotation_matrix: a matrix with shape (n,m) where n equals the number
                         of clusters and m equals the number of annotation
                         categories.
    """
    if groupby not in adata.obs.columns:
        raise KeyError(f'"{groupby}" not found.')
    if annot not in adata.obs.columns:
        raise KeyError(f'"{annot}" not found.')
    groupby_annot = f'{groupby}_annot'
    clustering = adata.obs[groupby]
    annotation = adata.obs[annot]
    k_nan = annotation == 'nan'
    cell_types = np.unique(annotation[~k_nan].values)
    n_celltype = cell_types.shape[0]
    n_cluster = np.unique(clustering).shape[0]
    n_cell = clustering.shape[0]

    def dedup(names):
        names_copy = names.copy()
        count = {}
        for i, name in enumerate(names):
            n = count.get(name, 0) + 1
            count[name] = n
            if n > 1 or name in names[(i+1):]:
                names_copy[i] = f'{name}{n}'
            else:
                names_copy[i] = name
        return names_copy

    if (annotation_matrix.shape[0] == n_cluster and
            annotation_matrix.shape[1] == n_celltype):
        k_max = np.argmax(annotation_matrix, axis=1)
        cluster_annot = cell_types[k_max]
        adata.obs[groupby_annot] = adata.obs[groupby].cat.rename_categories(
                dedup(cluster_annot))
    else:
        raise ValueError(f'[annotation_matrix] not in compatible shape')


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


def plot_qc(adata, groupby=None):
    qc_metrics = ('n_counts', 'n_genes', 'percent_mito')
    for qmt in qc_metrics:
        if qmt not in adata.obs.columns:
            raise ValueError(f'{qmt} not found.')
    if groupby:
        sc.pl.violin(adata, keys=qc_metrics, groupby=groupby, rotation=45)
    sc.pl.violin(adata, keys=qc_metrics, multi_panel=True, rotation=45)
    sc.pl.scatter(
        adata, x='n_counts', y='n_genes', color='percent_mito', alpha=0.5)
    sc.pl.scatter(adata, x='n_counts', y='n_genes', color=groupby, alpha=0.5)
    sc.pl.scatter(
        adata, x='n_counts', y='percent_mito', color=groupby, alpha=0.5)


def plot_metric_by_rank(
        adata,
        subject='cell',
        metric='n_counts',
        decreasing=True,
        ax=None,
        title=None,
        hpos=None,
        vpos=None,
        logx=True,
        logy=True,
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
        adata.obs['mt_prop'] = np.squeeze(
            np.asarray(adata.X[:, k_mt].sum(axis=1))) / adata.obs['n_counts']

    order_modifier = -1 if decreasing else 1

    if subject == 'cell':
        x = 1 + np.arange(adata.shape[0])
        k = np.argsort(order_modifier * adata.obs[metric])
        y = adata.obs[metric].values[k]
    else:
        x = 1 + np.arange(adata.shape[1])
        k = np.argsort(order_modifier * adata.var[metric])
        y = adata.var[metric].values[k]

    if not ax:
        fig, ax = plt.subplots()

    if kwargs['c'] is not None and not isinstance(kwargs['c'], str):
        kwargs_c = kwargs['c'][k]
        del kwargs['c']
        for c in np.unique(kwargs_c):
            k = kwargs_c == c
            ax.plot(x[k], y[k], c=c, **kwargs)
    else:
        ax.plot(x, y, **kwargs)

    if hpos is not None:
        ax.hlines(hpos, xmin=x.min(), xmax=x.max(), linewidth=1, colors='red')
    if vpos is not None:
        ax.vlines(vpos, ymin=y.min(), ymax=y.max(), linewidth=1, colors='green')
    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')
    ax.set_xlabel('Rank')
    ax.set_ylabel(metric_names[metric])
    if title:
        ax.set_title(title)
    ax.grid(linestyle='--')
    if 'label' in kwargs:
        ax.legend(loc='upper right' if decreasing else 'upper left')

    return ax


def plot_embedding(
        adata, basis, groupby, annot=True, highlight=None, size=None,
        save=None, savedpi=300, **kwargs):
    if f'X_{basis}' not in adata.obsm.keys():
        raise KeyError(f'"X_{basis}" not found in `adata.obsm`.')
    if isinstance(groupby, (list, tuple)):
        groupby = groupby[0]
    if groupby not in adata.obs.columns:
        raise KeyError(f'"{groupby}" not found in `adata.obs`.')
    if adata.obs[groupby].dtype.name != 'category':
        raise ValueError(f'"{groupby}" is not categorical.')
    from scanpy.plotting._tools.scatterplots import plot_scatter
    groups = adata.obs[groupby].copy()
    categories = list(adata.obs[groupby].cat.categories)
    rename_dict = {ct: f'{i+1}: {ct}' for i, ct in enumerate(categories)}
    restore_dict = {f'{i+1}: {ct}': ct for i, ct in enumerate(categories)}

    size_ratio = 1.2

    ad = adata
    marker_size = size
    kwargs['show'] = False
    kwargs['save'] = False
    kwargs['frameon'] = True
    kwargs['legend_loc'] = 'right margin'

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
        plot_scatter(ad, basis, color=groupby, ax=ax, size=marker_size, **kwargs)

        if highlight:
            embed = adata.obsm[f'X_{basis}']
            for i, ct in enumerate(categories):
                if ct not in highlight:
                    continue
                k_hl = groups == ct
                ax.scatter(embed[k_hl, 0], embed[k_hl, 1], marker='D', c=adata.uns[f'{groupby}_colors'][i], s=size/8)
                ax.scatter([], [], marker='D', c=adata.uns[f'{groupby}_colors'][i], label=adata.obs[groupby].cat.categories[i])
            if annot:
                ax.legend(frameon=False, loc='center left', bbox_to_anchor=(1, 0.5),
                          ncol=(1 if len(categories) <= 14 else 2 if len(categories) <= 30 else 3))
    finally:
        if annot:
            adata.obs[groupby].cat.rename_categories(restore_dict, inplace=True)
    if annot:
        centroids = pseudo_bulk(adata, groupby, use_rep=f'X_{basis}', FUN=np.median).T
        texts = [ax.text(x=row[0], y=row[1], s=f'{i+1:d}', fontsize=8, fontweight='bold') for i, row in centroids.reset_index(drop=True).iterrows()]
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, text_from_points=False, autoalign=False)
    if save:
        plt.savefig(fname=save, dpi=savedpi, bbox_inches='tight', pad_inches=0.1)


def simple_default_pipeline(
        adata,
        qc_only=False,
        batch=None,
        filter_params={'min_genes': 200, 'min_cells': 3, 'max_counts': 25000, 'max_mito': 20},
        norm_params={'target_sum': 1e4, 'fraction': 0.9},
        combat_args={'key': None},
        hvg_params={'flavor': 'seurat', 'by_batch': None},
        scale_params={'max_value': 10},
        pca_params={'n_comps': 50, 'svd_solver': 'arpack', 'use_highly_variable': True},
        harmony_params={'batch': None, 'theta': 2.0},
        nb_params={'n_neighbors': 15, 'n_pcs': 40},
        umap_params={},
        tsne_params={},
        diffmap_params={'n_comps': 15},
        leiden_params={
            'resolution': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]},
        fdg_params={'layout': 'fa'},
):
    """
    Scanpy pipeline
    """
    if qc_only:
        adata.var['mito'] = adata.var_names.str.startswith('MT-')
        qc_tbls = sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mito'], percent_top=None)
        adata.obs['n_counts'] = qc_tbls[0]['total_counts'].values
        adata.obs['n_genes'] = qc_tbls[0]['n_genes_by_counts'].values
        adata.obs['percent_mito'] = qc_tbls[0]['pct_counts_mito'].values
        adata.var['n_cells'] = qc_tbls[1]['n_cells_by_counts'].values
        plot_qc(adata, batch)
    else:
        if filter_params is not None and isinstance(filter_params, dict):
            if 'min_genes' in filter_params:
                sc.pp.filter_cells(adata, min_genes=filter_params['min_genes'])
            if 'min_cells' in filter_params:
                sc.pp.filter_genes(adata, min_cells=filter_params['min_cells'])
            if 'min_counts' in filter_params:
                k = adata.obs['n_counts'] >= filter_params['min_counts']
                adata._inplace_subset_obs(k)
            if 'max_counts' in filter_params:
                k = adata.obs['n_counts'] <= filter_params['max_counts']
                adata._inplace_subset_obs(k)
            if 'max_mito' in filter_params:
                k = adata.obs['percent_mito'] <= filter_params['max_mito']
                adata._inplace_subset_obs(k)
            if 'counts' not in adata.layers.keys():
                adata.layers['counts'] = adata.X
        if norm_params is not None and isinstance(norm_params, dict):
            if 'counts' in adata.layers.keys():
                adata.X = adata.layers['counts']
            sc.pp.normalize_total(adata, **norm_params)
            sc.pp.log1p(adata)
            adata.raw = adata
        if (combat_args is not None and (
                isinstance(combat_args, dict) and
                combat_args.get('key', None) and
                combat_args['key'] in adata.obs.keys())):
            adata.layers['X'] = adata.X
            adata.X = adata.raw.X
            sc.pp.combat(adata, **combat_args)
        if hvg_params is not None and isinstance(hvg_params, dict):
            hvg(adata, **hvg_params)
        if scale_params is not None and isinstance(scale_params, dict):
            sc.pp.scale(adata, **scale_params)
        if pca_params is not None and isinstance(pca_params, dict):
            pca(adata, **pca_params)
        if (harmony_params is not None and (
                isinstance(harmony_params, dict) and
                harmony_params.get('batch', None))):
            run_harmony(adata, **harmony_params)
        if nb_params is not None and isinstance(nb_params, dict):
            neighbors(adata, **nb_params)
        if umap_params is not None and isinstance(umap_params, dict):
            umap(adata, **umap_params)
        if tsne_params is not None and isinstance(tsne_params, dict):
            tsne(adata, **tsne_params)
        if diffmap_params is not None and isinstance(diffmap_params, dict):
            diffmap(adata, **diffmap_params)
        if leiden_params is not None and isinstance(leiden_params, dict):
            leiden(adata, **leiden_params)
        if fdg_params is not None and isinstance(fdg_params, dict):
            fdg(adata, **fdg_params)
    return adata
