# scanpy-cli
Scripts for using scanpy from the command line

In order to wrap scanpy's internal workflow in any given workflow language, it's important to have scripts to call each of those steps. These scripts are being written here, and will improve in completeness as time progresses. 

## Install

```
pip3 install scanpy-cli
```

## Test installation

There is an example script included:

```
scanpy-cli-tests.sh
```

This downloads [a well-known test 10X dataset]('https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) and executes all of the scripts described below.

## Commands

Currently wrapper scripts are described below. Each script has usage instructions available via --help, consult function documentation in scanpy for further details.

### read-10x: read 10X data and create the AnnData object (calls `sc.read()`)

```
scanpy read-10x -d <10x data directory> -o <raw data object in .h5ad format>
```

### filter-cells: filter out poor-quality cells (calls `sc.pp.filter_cells()`)

```
scanpy filter-cells -i <raw data object in .h5ad format> -p n_genes,n_counts -l <min_genes>,<min_counts> -o <cell-filtered object in .h5ad format>
``` 

### filter-genes: filter out poorly covered genes (calls `sc.pp.filter_genes()`)

```
scanpy filter-genes -i <cell-filtered object in .h5ad format> -p n_cells,n_counts -l <min_celss>,<min_counts> -o <cell-and-gene-filtered object in .h5ad format>
``` 

### normalise-data: normalise the expression values (calls `sc.pp.normalize_per_cell()`)

```
scanpy normalise-data.py -i <cell-and-gene-filtered object in .h5ad format> -s <scale_factor> -o <object with normalised expression values in .h5ad format> [--save-raw]
```

### find-variable-genes: find variable genes (calls `sc.pp.filter_genes_dispersion()`)

```
scanpy find-variable-genes -i <object with normalised expression values in .h5ad format> --flavor <method to normalise dispersion> -p mean,disp -l <min_mean>,<min_disp> -j <high_mean>,<high_disp> -b <n bins> -n <n top genes> -o <object with variable genes in .h5ad format>
```

### scale-data: scale expression values (calls `sc.pp.scale()`, optionally `sc.pp.log1p()` and `sc.pp.regress_out()`)

```
scanpy scale-data -i <object with variable genes in .h5ad format>  -V <variables to regress> -x <scale max> -o <object with scaled expression values in .h5ad format>
```

### run-pca: run principal components analysis (calls `sc.tl.pca()`, `sc.pl.pca()`)

```
scanpy run-pca -i <object with scaled expression values in .h5ad format> -n <number of pcs to compute> -o <output object with PCs in .h5ad format> --output-embeddings-file <pca embedding in csv format> -output-loadings-file <pca loadings file in csv format> --output-stdev-file <pca stdev file in text format> --output-var-ratio-file <pca proportion of explained variance file in text format> -P <PCA plot file> --color <variable to color cells by>
```

### neighbours: compute neighbourhood graph (calls `sc.pp.neighbors()`)

```
scanpy neighbours -i <object with PCs .h5ad format> -N <number of neighbors to consider> -n <number of PCs to use> -m <method to compute connectivity> -M <distance metric> -o <output object with neighbourhood graph in .h5ad format>
```

### find-cluster: find clusters (calls `sc.tl.louvain()`)

```
scanpy find-cluster -i <object with neighbourhood graph in .h5ad format> --flavor <method to compute clustering> -o <output object with clusters in .h5ad format> --output-text-file <output cluster assignment table in csv table>
```

### run-umap: run UMAP analysis (calls `sc.tl.umap()`, `sc.pl.umap()`)

```
scanpy run-umap -i <object with clusters in .h5ad format> -n <number of dimensions> -o <output object with umap embeddings in .h5ad format> --output-embeddings-file <umap embeddings in csv format> -P <umap plot file> --color <variable to color cells by>
```

### run-tsne: run tSNE analysis (calls `sc.tl.tsne()`, `sc.pl.tsne()`)

```
scanpy run-tsne -i <object with clusters in .h5ad format> -n <number of dimensions> -o <output object with tsne embeddings in .h5ad format> --output-embeddings-file <tsne embeddings in csv format> -P <tsne plot file> --color <variable to color cells by>
```

### find-markers: find marker genes for each group/cluster of cells (calls `sc.tl.rank_genes_groups()`)

```
scanpy find-clusters -i <object with clusters in .h5ad format> -g <groupby> -n <number of genes to test for each group> -m <method of testing> --reference <reference group to compare agains> -o <output object in .h5ad format> --output-text-file <table of top tested candidate marker genes in csv format> -P <plot of candidate gene expression across groups> --show-n-genes <number of genes to plot>
```
