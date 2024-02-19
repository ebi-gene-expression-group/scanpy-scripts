"""
scanpy diffexp
"""

import logging
import math

import pandas as pd
import scanpy as sc


def diffexp(
    adata,
    use_raw=None,
    n_genes=None,
    key_added="rank_genes_groups",
    layer=None,
    logreg_param=None,
    filter_params=None,
    save=None,
    groupby=None,
    groups=None,
    **kwargs,
):
    """
    Wrapper function for sc.tl.rank_genes_groups.

    Test that we can load a single group.
    >>> import os
    >>> from pathlib import Path
    >>> adata = sc.datasets.krumsiek11()
    >>> tbl = diffexp(adata, groupby='cell_type', groups='Mo', reference='progenitor')
    >>> # get the size of the data frame
    >>> tbl.shape
    (11, 8)
    """
    if adata.raw is None:
        use_raw = False

    if n_genes is None:
        n_genes = adata.raw.shape[1] if use_raw else adata.shape[1]

    if logreg_param and isinstance(logreg_param, dict):
        for key, val in logreg_param:
            kwargs[key] = val

    key_added = key_added if key_added else "rank_genes_groups"
    diff_key = (key_added + f"_{layer}") if layer else key_added

    if groups == "all":

        # Avoid divisions by zeros for singlet groups. See
        # https://github.com/theislab/scanpy/pull/1490#issuecomment-726031442.

        groups_to_test = list(
            adata.obs[groupby].value_counts().loc[lambda x: x > 1].index
        )

        if len(groups_to_test) < len(adata.obs[groupby].cat.categories):
            groups = groups_to_test
            logging.warning(
                "Singlet groups removed before passing to rank_genes_groups()"
            )

    # avoid issue when groups is a single group as a string simplified by click
    # https://github.com/ebi-gene-expression-group/scanpy-scripts/issues/123
    if groups != "all" and isinstance(groups, str):
        groups = [groups]

    sc.tl.rank_genes_groups(
        adata,
        use_raw=use_raw,
        n_genes=n_genes,
        key_added=diff_key,
        groupby=groupby,
        groups=groups,
        **kwargs,
    )

    de_tbl = extract_de_table(adata.uns[diff_key])

    if isinstance(filter_params, dict):
        key_filtered = diff_key + "_filtered"
        sc.tl.filter_rank_genes_groups(
            adata,
            key=diff_key,
            key_added=key_filtered,
            use_raw=use_raw,
            **filter_params,
        )

        # there are non strings on recarray object at this point, in
        # adata.uns['rank_genes_groups_filtered']['names']
        # for instance:
        # adata.uns['rank_genes_groups_filtered']['names'][0]
        # (nan, nan, 'NKG7', nan, nan, 'PPBP')
        # this now upsets h5py > 3.0
        de_tbl = extract_de_table(adata.uns[key_filtered])
        de_tbl = de_tbl.loc[de_tbl.genes.astype(str) != "nan", :]

        # change nan for strings in adata.uns['rank_genes_groups_filtered']['names']
        # TODO on scanpy updates, check if this is not done within scanpy so that we can remove this
        for row in range(0, len(adata.uns[key_filtered]["names"])):
            for col in range(0, len(adata.uns[key_filtered]["names"][row])):
                element = adata.uns[key_filtered]["names"][row][col]
                if isinstance(element, float) and math.isnan(element):
                    adata.uns[key_filtered]["names"][row][col] = "nan"

    if save:
        de_tbl.to_csv(save, sep="\t", header=True, index=False)

    return de_tbl


def diffexp_paired(adata, groupby, pair, **kwargs):
    """
    Restrict DE to between a pair of clusters, return both up and down genes
    """
    test, ref = pair
    de_key = f"de.{test}-{ref}"
    up_de = diffexp(
        adata,
        key_added=de_key,
        groupby=groupby,
        groups=[test],
        reference=ref,
        **kwargs,
    )
    ref, test = pair
    de_key = f"de.{test}-{ref}"
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
    if de_dict["params"]["method"] == "logreg":
        requested_fields = ("scores",)
    else:
        requested_fields = (
            "scores",
            "logfoldchanges",
            "pvals",
            "pvals_adj",
        )
    gene_df = _recarray_to_dataframe(de_dict["names"], "genes")[
        ["cluster", "rank", "genes"]
    ]
    gene_df["ref"] = de_dict["params"]["reference"]
    gene_df = gene_df[["cluster", "ref", "rank", "genes"]]
    de_df = pd.DataFrame(
        {
            field: _recarray_to_dataframe(de_dict[field], field)[field]
            for field in requested_fields
            if field in de_dict
        }
    )
    return gene_df.merge(de_df, left_index=True, right_index=True)


def _recarray_to_dataframe(array, field_name):
    return (
        pd.DataFrame(array)
        .reset_index()
        .rename(columns={"index": "rank"})
        .melt(id_vars="rank", var_name="cluster", value_name=field_name)
    )
