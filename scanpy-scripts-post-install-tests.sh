#!/usr/bin/env bash

script_dir=$(realpath "$(dirname "${BASH_SOURCE[0]}" )")
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
export raw_scanpy_object="$output_dir/raw_scanpy.loom"
export filtered_cells_object="$output_dir/filtered_cells.loom"
export filtered_genes_object="$output_dir/filtered_genes.loom"
export normalised_per_cell_object="$output_dir/normalised_per_cell.loom"
export variable_genes_object="$output_dir/variable_genes.loom"
export variable_genes_list="$output_dir/variable_genes.txt"
export test_genes="$output_dir/random_genes.txt"
export scaled_seurat_object="$output_dir/scaled_seurat.rds"
export pca_seurat_object="$output_dir/pca_seurat.rds"
export pca_embeddings_file="$output_dir/pca_embeddings.csv"
export pca_loadings_file="$output_dir/pca_loadings.csv"
export pca_stdev_file="$output_dir/pca_stdev.txt"
export pca_image_file="$output_dir/pcatest.png"
export cluster_seurat_object="$output_dir/cluster_seurat.rds"
export cluster_text_file="$output_dir/clusters.txt"
export tsne_seurat_object="$output_dir/tsne_seurat.rds"
export tsne_embeddings_file="$output_dir/tsne_embeddings.csv"
export marker_text_file="$output_dir/markers.csv"

## Test parameters- would form config file in real workflow. DO NOT use these
## as default values without being sure what they mean.

# Normalisation.

export min_genes=200
export min_gene_counts=1000

# Normalisation. See Seurat ?NormalizeData

export transformation_method='LogNormalize'
export scale_factor=10000

# Find variable genes. See ?FindVariableGenes

export mean_function='ExpMean'
export dispersion_function='LogVMR'
export low_mean='0.0125'
export high_mean='3'
export low_disp='0.5'
export high_disp='Inf'

# Scale and center the data. See ?ScaleData
export vars_to_regress='nUMI'
export model_use='linear'
export use_umi='TRUE'
export do_scale='TRUE'
export do_center='TRUE'
export scale_max='10'
export block_size='1000'
export min_cells_to_block='1000'
export check_for_norm='TRUE'

# Run PCA
export pcs_compute=50
export use_imputed='FALSE'

# Plot PCA
export pca_dim_one=1
export pca_dim_two=2
export pt_size=1
export label_size=4
export do_label='FALSE'
export group_by='ident'
export pca_plot_title='Test PCA plot'
export pca_png_width=1000
export pca_png_height=1000

# Find clusters
export reduction_type='pca'
export dims_use='1,2,3,4,5,6,7,8,9,10'
export k_param=30
export resolution=0.8
export cluster_algorithm=1
export cluster_tmp_file_location='/tmp'

# t-SNE
export tsne_do_fast='TRUE'

# Marker detection
export logfc_threshold=0.25
export marker_min_pct=0.1
export marker_only_pos='FALSE'
export marker_test_use='wilcox'
export marker_max_cells_per_ident='Inf'
export marker_min_cells_gene=3
export marker_min_cells_group=3

################################################################################
# Test individual scripts
################################################################################

# Make the script options available to the tests so we can skip tests e.g.
# where one of a chain has completed successfullly.

export use_existing_outputs

# Derive the tests file name from the script name

tests_file="${script_name%.*}".bats

# Execute the bats tests

$tests_file
