"""
scanpy filter
"""

import logging
import re
import click
import numpy as np
import scanpy as sc


def filter_anndata(
        adata,
        gene_name='index',
        list_attr=False,
        param=None,
        category=None,
        subset=None,
):
    """
    Wrapper function for sc.pp.filter_cells() and sc.pp.filter_genes(), mainly
    for supporting arbitrary filtering
    """
    param = [] if param is None else param
    category = [] if category is None else category
    subset = [] if subset is None else subset

    logging.debug('--gene-name=%s', gene_name)
    logging.debug('--param=%s', param)
    logging.debug('--category=%s', category)
    logging.debug('--subset=%s', subset)

    if gene_name:
        try:
            gene_names = getattr(adata.var, gene_name)
            adata.var['mito'] = gene_names.str.startswith('MT-')
        except AttributeError:
            logging.warning(
                'Specified gene column [%s] not found, skip calculating '
                'expression of mitochondria genes', gene_name)

    attributes = _get_attributes(adata)
    if list_attr:
        click.echo(_repr_obj(attributes))
        return 0

    conditions, qc_vars, pct_top = _get_filter_conditions(
        adata, attributes, param, category, subset)

    sc.pp.filter_cells(adata, min_genes=0)
    sc.pp.filter_cells(adata, min_counts=0)
    sc.pp.filter_genes(adata, min_cells=0)
    sc.pp.filter_genes(adata, min_counts=0)

    if qc_vars or pct_top:
        if not pct_top:
            pct_top = [50]
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=qc_vars, percent_top=pct_top, inplace=True)

    k_cell = np.ones(len(adata.obs)).astype(bool)
    for cond in conditions['cell']['numerical']:
        name, vmin, vmax = cond
        if name.startswith('cell:'):
            name = name[5:]
        attr = getattr(adata.obs, name)
        k_cell = k_cell & (attr >= vmin) & (attr <= vmax)

    for cond in conditions['cell']['categorical']:
        name, values = cond
        attr = getattr(adata.obs, name)
        k_cell = k_cell & attr.isin(values)

    k_gene = np.ones(len(adata.var)).astype(bool)
    for cond in conditions['gene']['numerical']:
        name, vmin, vmax = cond
        if name.startswith('gene:'):
            name = name[5:]
        attr = getattr(adata.var, name)
        k_gene = k_gene & (attr >= vmin) & (attr <= vmax)

    for cond in conditions['gene']['categorical']:
        name, values = cond
        attr = getattr(adata.var, name)
        k_gene = k_gene & attr.isin(values)

    adata = adata[k_cell, :]
    adata = adata[:, k_gene]

    return adata


def _get_attributes(adata):
    attributes = {
        'cell': {
            'numerical': [],
            'categorical': ['cell:index'],
            'bool': [],
        },
        'gene': {
            'numerical': [],
            'categorical': ['gene:index'],
            'bool': [],
        },
    }

    for attr, dtype in adata.obs.dtypes.to_dict().items():
        typ = dtype.kind
        if typ == 'O':
            attributes['cell']['categorical'].append(attr)
        elif typ in ('i', 'f'):
            attributes['cell']['numerical'].append(attr)
        elif typ == 'b':
            attributes['cell']['bool'].append(attr)

    for attr, dtype in adata.var.dtypes.to_dict().items():
        typ = dtype.kind
        if typ == 'O':
            attributes['gene']['categorical'].append(attr)
        elif typ in ('i', 'f'):
            attributes['gene']['numerical'].append(attr)
        elif typ == 'b':
            attributes['gene']['bool'].append(attr)

    attributes['cell']['numerical'].extend([
        'n_genes',
        'cell:n_counts',
        'pct_counts_in_top_<n>_genes',
    ])

    for attr in attributes['gene']['bool']:
        attributes['cell']['numerical'].append('pct_counts_' + attr)

    attributes['gene']['numerical'].extend([
        'n_cells',
        'gene:n_counts',
        'mean_counts',
        'pct_dropout_by_counts',
    ])
    logging.debug(attributes)
    return attributes


def _get_filter_conditions(adata, attributes, param, category, subset):
    conditions = {
        'cell': {
            'numerical': [],
            'categorical': [],
            'bool': [],
        },
        'gene': {
            'numerical': [],
            'categorical': [],
            'bool': [],
        },
    }
    percent_top_pattern = re.compile(r'^pct_counts_in_top_(?P<n>\d+)_genes$')
    pct_top = []
    qc_vars_pattern = re.compile(r'^pct_counts_(?P<qc_var>\S+)$')
    qc_vars = []
    for name, vmin, vmax in param:
        if name in ('n_counts', 'index'):
            logging.error(
                'Ambiguous parameter name [%s] given, choose from '
                '"gene:{%s}" or "cell:{%s}"', name, name, name)
            raise ValueError
        pt_match = percent_top_pattern.match(name)
        if pt_match:
            if name not in adata.obs.columns:
                pct_top.append(int(pt_match['n']))
            conditions['cell']['numerical'].append([name, vmin*100, vmax*100])
            continue
        qv_match = qc_vars_pattern.match(name)
        if qv_match and qv_match['qc_var'] in attributes['gene']['bool']:
            if name not in adata.obs.columns:
                qc_vars.append(qv_match['qc_var'])
            conditions['cell']['numerical'].append([name, vmin*100, vmax*100])
            continue
        if name in attributes['cell']['numerical']:
            conditions['cell']['numerical'].append([name, vmin, vmax])
        elif name in attributes['gene']['numerical']:
            conditions['gene']['numerical'].append([name, vmin, vmax])
        else:
            logging.warning('Parameter [%s] undefined, '
                            'dropped from filtering', name)
    for ct in category + subset:
        name = ct[0]
        if not isinstance(ct[1], list):
            fh = ct[1]
            values = fh.read().rstrip().split('\n')
            fh.close()
            ct[1] = values
        if name in attributes['cell']['categorical']:
            conditions['cell']['categorical'].append(ct)
        elif name in attributes['gene']['categorical']:
            conditions['gene']['categorical'].append(ct)
        else:
            logging.warning('Attribute [%s] undefined, '
                            'dropped from filtering', name)
    logging.debug((conditions, qc_vars, pct_top))
    return conditions, qc_vars, sorted(pct_top)


def _repr_obj(obj, padding='  ', level=0):
    if isinstance(obj, dict):
        obj_str = '\n'.join(['\n'.join([
            padding * level + k + ':', _repr_obj(v, level=level+1)
        ]) for k, v in obj.items()])
    elif isinstance(obj, (tuple, list, set)):
        obj_str = '\n'.join([_repr_obj(elm, level=level) for elm in obj])
    else:
        obj_str = padding * level + repr(obj)
    return obj_str
