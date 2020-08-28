# scanpy-scripts [![Anaconda-Server Badge](https://anaconda.org/bioconda/scanpy-scripts/badges/installer/conda.svg)](https://anaconda.org/bioconda/scanpy-scripts) 

A command-line interface for functions of the Scanpy suite, to facilitate flexible constrution of workflows, for example in Galaxy, Nextflow, Snakemake etc.

## Install

```bash
conda install scanpy-scripts
# or
pip3 install scanpy-scripts
```

## Test installation

There is an example script included:

```bash
scanpy-scripts-tests.bats
```

This requires the [bats](https://github.com/sstephenson/bats) testing framework to run. The script downloads [a well-known test 10X dataset]('https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) and executes all of the commands described below.

## Commands

Available commands are described below. Each has usage instructions available via --help, consult function documentation in scanpy for further details.

```
Usage: scanpy-cli [OPTIONS] COMMAND [ARGS]...

  Command line interface to [scanpy](https://github.com/theislab/scanpy)

Options:
  --debug              Print debug information
  --verbosity INTEGER  Set scanpy verbosity
  --version            Show the version and exit.
  --help               Show this message and exit.

Commands:
  read       Read 10x data and save in specified format.
  filter     Filter data based on specified conditions.
  norm       Normalise data per cell.
  hvg        Find highly variable genes.
  scale      Scale data per gene.
  regress    Regress-out observation variables.
  pca        Dimensionality reduction by PCA.
  neighbor   Compute a neighbourhood graph of observations.
  embed      Embed cells into two-dimensional space.
  cluster    Cluster cells into sub-populations.
  diffexp    Find markers for each clusters.
  paga       Trajectory inference by abstract graph analysis.
  dpt        Calculate diffusion pseudotime relative to the root cells.
  integrate  Integrate cells from different experimental batches.
  plot       Visualise data.
  ```
