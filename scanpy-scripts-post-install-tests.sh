#!/usr/bin/env bash

script_dir=$(readlink -f "$(dirname "${BASH_SOURCE[0]}" )")
script_name=$(basename "${BASH_SOURCE[0]}")

function usage {
    echo "usage: scanpy-scripts-post-install-tests.sh [action] [use_existing_outputs]"
    echo "  - action: what action to take 'test' or 'clean'"
    echo "  - use_existing_outputs: 'true' or 'false'"
    exit 1
}

action=${1:-'test'}
use_existing_outputs=${2:-'false'}

if [[ $action != 'test' ]] && [[ $action != 'clean' ]]; then
    echo "Invalid action"
    usage
fi

if [[ $use_existing_outputs != 'true' ]] && [[ $use_existing_outputs != 'false' ]]; then
    echo "Invalid value ($use_existing_outputs) for 'use_existing_outputs'"
    usage
fi

test_data_url='https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
test_working_dir=$script_dir/post_install_tests
export test_data_archive=$test_working_dir/$(basename "$test_data_url")

# Clean up if specified

if [[ $action == 'clean' ]]; then
    echo "Cleaning up $test_working_dir ..."
    rm -rf $test_working_dir
    exit 0
elif [[ $action != 'test' ]]; then
    echo "Invalid action '$action' supplied"
    exit 1
fi

# Initialise directories

output_dir=$test_working_dir/outputs
export data_dir=$test_working_dir/data

mkdir -p $test_working_dir
mkdir -p $output_dir
mkdir -p $data_dir

################################################################################
# Fetch test data
################################################################################

if [[ ! -e $test_data_archive ]]; then
    wget $test_data_url -P $test_working_dir

fi

################################################################################
# List tool outputs/ inputs
################################################################################

export raw_matrix="$data_dir/matrix.mtx"
export input_object="$output_dir/input.h5ad"
export filtered_cells_object="$output_dir/filtered_cells.h5ad"
export filtered_genes_object="$output_dir/filtered_genes.h5ad"
export normalised_object="$output_dir/normalised.h5ad"
export variable_genes_object="$output_dir/variable_genes.h5ad"
export variable_image_file="$output_dir/variable_genes.png"
export scaled_object="$output_dir/scaled.h5ad"
export pca_object="$output_dir/pca.h5ad"
export pca_image_file="$output_dir/pca.png"
export pca_embeddings_file="$output_dir/pca_embeddings.csv"
export pca_loadings_file="$output_dir/pca_loadings.csv"
export pca_stdev_file="$output_dir/pca_stdev.txt"
export pca_var_ratio_file="$output_dir/pca_var_ratio.txt"
export graph_object="$output_dir/graph.h5ad"
export cluster_object="$output_dir/cluster.h5ad"
export cluster_text_file="$output_dir/clusters.txt"
export umap_object="$output_dir/umap.h5ad"
export umap_embeddings_file="$output_dir/umap_embeddings.csv"
export umap_image_file="$output_dir/umap.png"
export tsne_object="$output_dir/tsne.h5ad"
export tsne_image_file="$output_dir/tsne.png"
export tsne_embeddings_file="$output_dir/tsne_embeddings.csv"
export marker_object="$output_dir/marker.h5ad"
export marker_image_file="$output_dir/markers.png"
export marker_text_file="$output_dir/markers.csv"

## Test parameters- would form config file in real workflow. DO NOT use these
## as default values without being sure what they mean.

# Filter cells.
export FC_parameters='n_genes'
export FC_min_genes=200
export FC_max_genes=2500

# Filter genes.
export FT_parameters='n_cells'
export FT_min_cells=3

# Normalisation.
export ND_scale_factor=10000
export ND_save_raw=True

# Find variable genes.
export FVG_flavor=seurat
export FVG_nbins=20
export FVG_parameters='mean,disp'
export FVG_low_mean='0.0125'
export FVG_high_mean=3
export FVG_low_disp='0.5'
export FVG_high_disp='Inf'

# Scale and center the data
export SD_vars_to_regress='n_counts'
export SD_zero_center='-z'
export SD_scale_max='10'

# Run PCA
export PCA_npcs=50
export PCA_svd_solver=arpack
export PCA_random_seed=0
export PCA_color=ENSG00000101439
export PCA_projection=2d
export PCA_frameon='--frameoff'

# Run compute graph
export CG_nneighbor=15
export CG_npcs=50
export CG_knn='--knn'
export CG_random_seed=0
export CG_method=umap

# Run find cluster
export FC_flavor=vtraag
export FC_resolution='1.0'
export FC_key_added=louvain
export FC_use_weight='--use-weights'
export FC_random_seed=0

# Run UMAP
export UMAP_min_dist='0.5'
export UMAP_spread='1.0'
export UMAP_use_raw='--no-raw'
export UMAP_ncomp=2
export UMAP_random_seed=0
export UMAP_alpha=1.0
export UMAP_gamma=1.0
export UMAP_initpos=spectral
export UMAP_color=louvain
export UMAP_projection=2d
export UMAP_frameon='--frameoff'

# Run t-SNE
export TSNE_perplexity=30
export TSNE_early_exaggeration=12
export TSNE_learning_rate=1000
export TSNE_random_seed=0
export TSNE_color=louvain
export TSNE_projection=2d
export TSNE_frameon='--frameoff'

# Run find marker
export FM_groupby=louvain
export FM_groups=all
export FM_reference=rest
export FM_n_genes=50
export FM_method='t-test_overestim_var'
export FM_show_n_genes=5
export FM_key='rank_genes_groups'

################################################################################
# Test individual scripts
################################################################################

# Make the script options available to the tests so we can skip tests e.g.
# where one of a chain has completed successfullly.

export use_existing_outputs

# Derive the tests file name from the script name

tests_file="${script_name%.*}".bats

# Execute the bats tests

./$tests_file
