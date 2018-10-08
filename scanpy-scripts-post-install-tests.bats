#!/usr/bin/env bats

# Extract the test data

@test "Extract .mtx matrix from archive" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$raw_matrix" ]; then
        skip "$raw_matrix exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $raw_matrix && tar -xvzf $test_data_archive --strip-components 2 -C $data_dir

    [ "$status" -eq 0 ]
    [ -f "$raw_matrix" ]
}

# Run scanpy-read-10x.py

@test "Scanpy AnnData object creation from 10x" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$raw_scanpy_object" ]; then
        skip "$raw_scanpy_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $raw_scanpy_object && bin/scanpy-read-10x.py -d $data_dir -o $raw_scanpy_object

    [ "$status" -eq 0 ]
    [ -f  "$raw_scanpy_object" ]
}


# Run scanpy-filter-cells.py

@test "Filter cells from a raw AnnData object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_cells_object" ]; then
        skip "$filtered_cells_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_cells_object && scanpy-filter-cells.py -i $raw_scanpy_object -s n_genes,n_counts -l $min_genes,$min_gene_counts -o $filtered_cells_object

    [ "$status" -eq 0 ]
    [ -f  "$filtered_cells_object" ]
}

# Run scanpy-filter-genes.py

@test "Filter genes from a cell-filtered AnnData object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_genes_object" ]; then
        skip "$filtered_genes_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_genes_object && scanpy-filter-genes.py -i $filtered_cells_object -o $filtered_genes_object

    [ "$status" -eq 0 ]
    [ -f  "$filtered_genes_object" ]
}

# Run scanpy-normalise-data.py

@test "Normalise expression values per cell" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$normalised_per_cell_object" ]; then
        skip "$normalised_per_cell_object exists and use_existing_outputs is set to 'true'"
    fi

    run scanpy-normalise-data.py -i $filtered_genes_object -t $transformation_method -s $scale_factor -o $normalised_per_cell_object

    [ "$status" -eq 0 ]
    [ -f  "$normalised_per_cell_object" ]
}

# Run scanpy-find-variable-genes.R

@test "Find variable genes" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$variable_genes_list" ]; then
        skip "$variable_genes_list exists and use_existing_outputs is set to 'true'"
    fi

    run scanpy-find-variable-genes.py -i $normalised_per_cell_object -m $mean_function -s mean,disp -l $min_mean,$min_disp -j $high_mean,$high_disp -o $variable_genes_object -v $variable_genes_list

    [ "$status" -eq 0 ]
    [ -f  "$variable_genes_list" ]
}

## # Get a random set of genes to use in testing argments to scale-data.R
##
## @test "Generate random gene list" {
##     if [ "$use_existing_outputs" = 'true' ] && [ -f "$test_genes" ]; then
##         skip "$test_genes exists and use_existing_outputs is set to 'true'"
##     fi
##
##     run seurat-get-random-genes.R $normalised_seurat_object $test_genes 10000
##
##     [ "$status" -eq 0 ]
##     [ -f  "$test_genes" ]
## }
##
## # Scale expression values
##
## @test "Scale expression values" {
##     if [ "$use_existing_outputs" = 'true' ] && [ -f "$scaled_seurat_object" ]; then
##         skip "$scaled_seurat_object exists and use_existing_outputs is set to 'true'"
##     fi
##
##     run seurat-scale-data.R -i $variable_genes_seurat_object -e $test_genes -v $vars_to_regress -m $model_use -u $use_umi -s $do_scale -c $do_center -x $scale_max -b $block_size -d $min_cells_to_block -a $assay_type -n $check_for_norm -o $scaled_seurat_object
##
##     [ "$status" -eq 0 ]
##     [ -f  "$scaled_seurat_object" ]
## }
##
## # Run PCA
##
## @test "Run principal component analysis" {
##     if [ "$use_existing_outputs" = 'true' ] && [ -f "$pca_seurat_object" ]; then
##         skip "$scaled_seurat_object exists and use_existing_outputs is set to 'true'"
##     fi
##
##     run seurat-run-pca.R -i $scaled_seurat_object -e $test_genes -p $pcs_compute -m $use_imputed -o $pca_seurat_object -b $pca_embeddings_file -l $pca_loadings_file -s $pca_stdev_file
##
##     [ "$status" -eq 0 ]
##     [ -f  "$pca_seurat_object" ]
## }
##
## # Plot the PCA
##
## @test "Plot dimension reduction" {
##     if [ "$use_existing_outputs" = 'true' ] && [ -f "$pca_image_file" ]; then
##         skip "$scaled_seurat_object exists and use_existing_outputs is set to 'true'"
##     fi
##
##     run seurat-dim-plot.r -i $pca_seurat_object -r pca -a $pca_dim_one -b $pca_dim_two -p $pt_size -l $label_size -d $do_label -f $group_by -t '$pca_plot_title' -w $pca_png_width -j $pca_png_height -o $pca_image_file
##
##     [ "$status" -eq 0 ]
##     [ -f  "$pca_image_file" ]
## }
##
## # Generate clusters
##
## @test "Generate cell clusters from expression values" {
##     if [ "$use_existing_outputs" = 'true' ] && [ -f "$cluster_text_file" ]; then
##         skip "$pca_image_file exists and use_existing_outputs is set to 'true'"
##     fi
##
##     run seurat-find-clusters.r -i $pca_seurat_object -e $test_genes -u $reduction_type -d $dims_use -k $k_param -r $resolution -a $cluster_algorithm -m $cluster_tmp_file_location -o $cluster_seurat_object -t $cluster_text_file
##
##     [ "$status" -eq 0 ]
##     [ -f  "$cluster_text_file" ]
## }
##
## # Run t-SNE
##
## @test "Run-tSNE analysis" {
##     if [ "$use_existing_outputs" = 'true' ] && [ -f "$tsne_seurat_object" ]; then
##         skip "$tsne_seurat_object exists and use_existing_outputs is set to 'true'"
##     fi
##
##     run seurat-run-tsne.r -i $pca_seurat_object -r $reduction_type -d $dims_use -e NULL -f $tsne_do_fast -o $tsne_seurat_object -b $tsne_embeddings_file
##
##     [ "$status" -eq 0 ]
##     [ -f  "$tsne_seurat_object" ]
## }
##
## # Run marker detection
##
## @test "Find markers for each cluster" {
##     if [ "$use_existing_outputs" = 'true' ] && [ -f "$marker_text_file" ]; then
##         skip "$marker_text_file exists and use_existing_outputs is set to 'true'"
##     fi
##
##     run seurat-find-markers.R -i $cluster_seurat_object -e NULL -l $logfc_threshold -m $marker_min_pct -p $marker_only_pos -t $marker_test_use -x $marker_max_cells_per_ident -c $marker_min_cells_gene -d $marker_min_cells_group -o $marker_text_file
##
##     [ "$status" -eq 0 ]
##     [ -f  "$marker_text_file" ]
## }