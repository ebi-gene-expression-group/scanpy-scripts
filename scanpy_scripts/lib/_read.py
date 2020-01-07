"""
Provides read_10x()
"""

import os
import pandas as pd
import scanpy as sc


def read_10x(
        input_10x_h5,
        input_10x_mtx,
        genome='hg19',
        var_names='gene_symbols',
        extra_obs=None,
        extra_var=None
):
    """
    Wrapper function for sc.read_10x_h5() and sc.read_10x_mtx(), mainly to
    support adding extra metadata
    """
    if input_10x_h5 is not None:
        adata = sc.read_10x_h5(input_10x_h5, genome=genome)
    elif input_10x_mtx is not None:
        adata = sc.read_10x_mtx(input_10x_mtx, var_names=var_names)

    if extra_obs:
        obs_tbl = pd.read_csv(extra_obs, sep='\t', header=0, index_col=0)
        adata.obs = adata.obs.merge(
            obs_tbl,
            how='left',
            left_index=True,
            right_index=True,
            suffixes=(False, False),
        )

    if extra_var:
        var_tbl = pd.read_csv(extra_var, sep='\t', header=0, index_col=0)
        adata.var = adata.var.merge(
            var_tbl,
            how='left',
            left_index=True,
            right_index=True,
            suffixes=(False, False),
        )
    return adata


def read_10x_atac(
        input_10x_mtx,
        extra_obs=None,
        extra_var=None
):
    """
    Wrapper function for sc.read_10x_mtx() to read 10x ATAC data
    """
    matrix_mtx = os.path.join(input_10x_mtx, 'matrix.mtx')
    barcodes_tsv = os.path.join(input_10x_mtx, 'barcodes.tsv')
    peaks_bed = os.path.join(input_10x_mtx, 'peaks.bed')
    if not os.path.exists(peaks_bed):
        raise FileNotFoundError

    import tempfile
    tmp_dir = tempfile.mkdtemp()
    genes_tsv = os.path.join(tmp_dir, 'genes.tsv')

    os.symlink(os.path.abspath(matrix_mtx), os.path.join(tmp_dir, 'matrix.mtx'))
    os.symlink(os.path.abspath(barcodes_tsv), os.path.join(tmp_dir, 'barcodes.tsv'))
    with open(peaks_bed) as f, open(genes_tsv, 'w') as fw:
        for line in f:
            fields = line.rstrip().split('\t')
            peak_id = '_'.join(fields)
            print(f'{peak_id}\t{peak_id}', file=fw)

    adata = sc.read_10x_mtx(tmp_dir)
    adata.var.rename({'gene_id': 'peak'}, inplace=True)

    os.remove(os.path.join(tmp_dir, 'matrix.mtx'))
    os.remove(os.path.join(tmp_dir, 'barcodes.tsv'))
    os.remove(genes_tsv)
    os.removedirs(tmp_dir)

    if extra_obs:
        obs_tbl = pd.read_csv(extra_obs, sep='\t', header=0, index_col=0)
        adata.obs = adata.obs.merge(
            obs_tbl,
            how='left',
            left_index=True,
            right_index=True,
            suffixes=(False, False),
        )

    if extra_var:
        var_tbl = pd.read_csv(extra_var, sep='\t', header=0, index_col=0)
        adata.var = adata.var.merge(
            var_tbl,
            how='left',
            left_index=True,
            right_index=True,
            suffixes=(False, False),
        )
    return adata
