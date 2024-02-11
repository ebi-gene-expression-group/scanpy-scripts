# scanpy-scripts [![Anaconda-Server Badge](https://anaconda.org/bioconda/scanpy-scripts/badges/installer/conda.svg)](https://anaconda.org/bioconda/scanpy-scripts) 

A command-line interface for functions of the Scanpy suite, to facilitate flexible constrution of workflows, for example in Galaxy, Nextflow, Snakemake etc.

## Install

The recommended way of installing scanpy-scripts is via conda:

```bash
conda install scanpy-scripts
```

pip installation is also possible, however the version of mnnpy is not patched as in the conda version, and so the `integrate` command will not work.

```bash
pip install scanpy-scripts
```

For development installation, we suggest following the github actions python-package.yml file.

Currently, tests run on python 3.9, so those are the recommended versions if not installing via conda. BKNN doesn't currently install on Python 3.10 due to a skip in Bioconda.

## How to update

Updating the stack that depends on this script is unfortunately complex, it gets all the way to Galaxy tools using this for production pipelines in at least two different institutions.

1. Open a branch here from develop
2. Change the conda dependencies versions desired in the test-env.yaml file (if any)
3. Make your changes in the branch (including adding any tests to the bats files for new areas or functionalities).
4. Iterate with the CI on the branch until all test pass with the new changes and dependencies.
5. Get reviews and merge to develop.
6. Open a bump branch in bioconda, making sure that it points to the develop branch of this repo rather than the pypi release and that all dependencies in bioconda reflect what is currently in the test-env.yaml file of this repo (which was used to test in the feature branch here in points 2 to 4).
7. Once the tests pass in bioconda, ask the bot to fetch the artifacts, download the linux artifact. Leave the bioconda branch in waiting, do not merge it (and probably add a message saying so)
8. Use the image .tar.gz inside the artifact to create a new local docker container on the linux machine where you intend to test the Galaxy tools.
9. Open a local feature branch on container-galaxy-sc repo to start trying the changes.
10. Change the scanpy-scripts macro2 requirements part to point to the newly created container.
11. Run planemo test, it should use the local container.

## Test installation

There is an example script included:

```bash
scanpy-scripts-tests.bats
```

This requires the [bats](https://github.com/sstephenson/bats) testing framework to run. The script downloads [a well-known test 10X dataset]('https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) and executes all of the commands described below.

## Commands

Available commands are described below. Each has usage instructions available via `--help`, consult function documentation in scanpy for further details.

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
  multiplet  Execute methods for multiplet removal.
  plot       Visualise data.
  ```
