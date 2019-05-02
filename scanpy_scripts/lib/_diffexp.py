"""
scanpy diffexp
"""

import pandas as pd
import scanpy as sc


def diffexp(
        adata,
        use_raw=True,
        n_genes=None,
        key_added=None,
        logreg_param=None,
        filter_params=None,
        save=None,
        **kwargs,
):
    """
    Wrapper function for sc.tl.rank_genes_groups.
    """
    if adata.raw is None:
        use_raw = False

    if n_genes is None:
        n_genes = adata.raw.shape[1] if use_raw else adata.shape[1]

    if logreg_param and isinstance(logreg_param, dict):
        for key, val in logreg_param:
            kwargs[key] = val

    sc.tl.rank_genes_groups(
        adata, use_raw=use_raw, n_genes=n_genes, key_added=key_added, **kwargs)

    key_added = key_added if key_added else 'rank_genes_groups'

    de_tbl = extract_de_table(adata.uns[key_added])

    if isinstance(filter_params, dict):
        sc.tl.filter_rank_genes_groups(
            adata,
            key=key_added,
            key_added=key_added + '_filtered',
            use_raw=use_raw,
            **filter_params,
        )

        de_tbl = extract_de_table(adata.uns[key_added + '_filtered'])
        de_tbl = de_tbl.loc[de_tbl.genes.astype(str) != 'nan', :]

    if save:
        filtered_de_tbl.to_csv(save, sep='\t', header=True, index=False)

    return de_tbl


def diffexp_paired(adata, groupby, pair, **kwargs):
    """
    Restrict DE to between a pair of clusters, return both up and down genes
    """
    test, ref = pair
    de_key = f'de.{test}-{ref}'
    up_de = diffexp(
        adata,
        key_added=de_key,
        groupby=groupby,
        groups=[test],
        reference=ref,
        **kwargs,
    )
    ref, test = pair
    de_key = f'de.{test}-{ref}'
    down_de = diffexp(
        adata,
        key_added=de_key,
        groupby=groupby,
        groups=[test],
        reference=ref,
        **kwargs,
    )
    return up_de, down_de


def extract_de_table(de_dict):
    """
    Extract DE table from adata.uns
    """
    if de_dict['params']['method'] == 'logreg':
        requested_fields = ('scores',)
    else:
        requested_fields = ('scores', 'logfoldchanges', 'pvals', 'pvals_adj',)
    gene_df = _recarray_to_dataframe(de_dict['names'], 'genes')[
        ['cluster', 'rank', 'genes']]
    gene_df['ref'] = de_dict['params']['reference']
    gene_df = gene_df[['cluster', 'ref', 'rank', 'genes']]
    de_df = pd.DataFrame({
        field: _recarray_to_dataframe(de_dict[field], field)[field]
        for field in requested_fields if field in de_dict
    })
    return gene_df.merge(de_df, left_index=True, right_index=True)


def _recarray_to_dataframe(array, field_name):
    return pd.DataFrame(array).reset_index().rename(
        columns={'index': 'rank'}).melt(
            id_vars='rank', var_name='cluster', value_name=field_name)
