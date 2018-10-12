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
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$input_object" ]; then
        skip "$input_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $input_object && bin/scanpy-read-10x.py -d $data_dir -o $input_object

    [ "$status" -eq 0 ]
    [ -f  "$input_object" ]
}


# Run scanpy-filter-cells.py

@test "Filter cells from a raw AnnData object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_cells_object" ]; then
        skip "$filtered_cells_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_cells_object && scanpy-filter-cells.py -i $input_object -o $filtered_cells_object -p $FC_parameters -l $FC_min_genes -j $FC_max_genes

    [ "$status" -eq 0 ]
    [ -f  "$filtered_cells_object" ]
}

# Run scanpy-filter-genes.py

@test "Filter genes from a cell-filtered AnnData object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_genes_object" ]; then
        skip "$filtered_genes_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_genes_object && scanpy-filter-genes.py -i $filtered_cells_object -o $filtered_genes_object -p $FT_parameters -l $FT_min_cells

    [ "$status" -eq 0 ]
    [ -f  "$filtered_genes_object" ]
}

# Run scanpy-normalise-data.py

@test "Normalise expression values per cell" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$normalised_per_cell_object" ]; then
        skip "$normalised_per_cell_object exists and use_existing_outputs is set to 'true'"
    fi

    run scanpy-normalise-data.py -i $filtered_genes_object -s $ND_scale_factor -o $normalised_per_cell_object

    [ "$status" -eq 0 ]
    [ -f  "$normalised_per_cell_object" ]
}

# Run scanpy-find-variable-genes.R

@test "Find variable genes" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$variable_genes_list" ]; then
        skip "$variable_genes_list exists and use_existing_outputs is set to 'true'"
    fi

    run scanpy-find-variable-genes.py -i $normalised_per_cell_object --flavor $FVG_flavor -b $FVG_nbins -p $FVG_parameters -l $FVG_min_mean,$FVG_min_disp -j $FVG_high_mean,$FVG_high_disp -o $variable_genes_object -P $variable_image_file

    [ "$status" -eq 0 ]
    [ -f  "$variable_genes_object" ] && [ -f "$variable_image_file" ]
}


# Scale expression values

@test "Scale expression values" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$scaled_object" ]; then
        skip "$scaled_object exists and use_existing_outputs is set to 'true'"
    fi

    run scanpy-scale-data.py -i $variable_genes_object -v $SD_vars_to_regress -c $SD_do_center -x $SD_scale_max -o $scaled_object

    [ "$status" -eq 0 ]
    [ -f  "$scaled_object" ]
}

# Run PCA

@test "Run principal component analysis" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$pca_object" ]; then
        skip "$scaled_object exists and use_existing_outputs is set to 'true'"
    fi

    run scanpy-run-pca.py -i $scaled_object -n $PCA_npcs --svd-solver $PCA_svd_solver -s $PCA_random_seed -o $pca_object -P $pca_image_file

    [ "$status" -eq 0 ]
    [ -f  "$pca_object" ] && [ -f "$pca_image_file" ]
}


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

# Local Variables:
# mode: sh
# End:
