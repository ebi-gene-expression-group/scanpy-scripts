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

@test "Scanpy Loom object creation from 10x" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$input_object" ]; then
        skip "$input_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $input_object && bin/scanpy-read-10x.py -d $data_dir -o $input_object

    [ "$status" -eq 0 ]
    [ -f  "$input_object" ]
}


# Run scanpy-filter-cells.py

@test "Filter cells from a raw Loom object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_cells_object" ]; then
        skip "$filtered_cells_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_cells_object && bin/scanpy-filter-cells.py -i $input_object -o $filtered_cells_object -p $FC_parameters -l $FC_min_genes -j $FC_max_genes

    [ "$status" -eq 0 ]
    [ -f  "$filtered_cells_object" ]
}

# Run scanpy-filter-genes.py

@test "Filter genes from a cell-filtered Loom object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_genes_object" ]; then
        skip "$filtered_genes_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_genes_object && bin/scanpy-filter-genes.py -i $filtered_cells_object -o $filtered_genes_object -p $FT_parameters -l $FT_min_cells

    [ "$status" -eq 0 ]
    [ -f  "$filtered_genes_object" ]
}

# Run scanpy-normalise-data.py

@test "Normalise expression values per cell" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$normalised_object" ]; then
        skip "$normalised_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $normalised_object && bin/scanpy-normalise-data.py -i $filtered_genes_object -s $ND_scale_factor -o $normalised_object

    [ "$status" -eq 0 ]
    [ -f  "$normalised_object" ]
}

# Run scanpy-find-variable-genes.R

@test "Find variable genes" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$variable_genes_object" ]; then
        skip "$variable_genes_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $variable_genes_object $variable_image_file && bin/scanpy-find-variable-genes.py -i $normalised_object --flavor $FVG_flavor -b $FVG_nbins -p $FVG_parameters -l $FVG_low_mean,$FVG_low_disp -j $FVG_high_mean,$FVG_high_disp -o $variable_genes_object -P $variable_image_file

    [ "$status" -eq 0 ]
    [ -f  "$variable_genes_object" ] && [ -f "$variable_image_file" ]
}


# Scale expression values

@test "Scale expression values" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$scaled_object" ]; then
        skip "$scaled_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm-f $scaled_object && bin/scanpy-scale-data.py -i $variable_genes_object -V $SD_vars_to_regress -x $SD_scale_max $SD_zero_center -o $scaled_object

    [ "$status" -eq 0 ]
    [ -f  "$scaled_object" ]
}

# Run PCA

@test "Run principal component analysis" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$pca_object" ]; then
        skip "$pca_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $pca_object $pca_image_file && bin/scanpy-run-pca.py -i $scaled_object -n $PCA_npcs --svd-solver $PCA_svd_solver -s $PCA_random_seed -o $pca_object -P $pca_image_file

    [ "$status" -eq 0 ]
    [ -f  "$pca_object" ] && [ -f "$pca_image_file" ]
}

# Run compute graph

@test "Run compute neighbor graph" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$graph_object" ]; then
        skip "$scaled_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $graph_object $graph_image_file && bin/scanpy-compute-graph.py -i $pca_object -N $CG_nneighbor -n $CG_npcs $CG_knn --random-seed $CG_random_seed --method $CG_method -s $CG_random_seed -o $graph_object

    [ "$status" -eq 0 ]
    [ -f  "$graph_object" ]
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