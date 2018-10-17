#!/usr/bin/env bats

# Extract the test data

@test "Extract .mtx matrix from archive" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$raw_matrix" ]; then
        skip "$raw_matrix exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $raw_matrix && \
	tar -xvzf $test_data_archive --strip-components 2 -C $data_dir

    [ "$status" -eq 0 ]
    [ -f "$raw_matrix" ]
}

# Run scanpy-read-10x.py

@test "Scanpy object creation from 10x" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$input_object" ]; then
        skip "$input_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $input_object && \
	bin/scanpy-read-10x.py -d $data_dir \
			       -o $input_object

    [ "$status" -eq 0 ]
    [ -f  "$input_object" ]
}


# Run scanpy-filter-cells.py

@test "Filter cells from a raw object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_cells_object" ]; then
        skip "$filtered_cells_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_cells_object && \
	bin/scanpy-filter-cells.py -i $input_object \
				   -o $filtered_cells_object \
				   -p $FC_parameters \
				   -l $FC_min_genes \
				   -j $FC_max_genes

    [ "$status" -eq 0 ]
    [ -f  "$filtered_cells_object" ]
}

# Run scanpy-filter-genes.py

@test "Filter genes from a cell-filtered object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_genes_object" ]; then
        skip "$filtered_genes_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_genes_object && \
	bin/scanpy-filter-genes.py -i $filtered_cells_object \
				   -o $filtered_genes_object \
				   -p $FT_parameters \
				   -l $FT_min_cells

    [ "$status" -eq 0 ]
    [ -f  "$filtered_genes_object" ]
}

# Run scanpy-normalise-data.py

@test "Normalise expression values per cell" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$normalised_object" ]; then
        skip "$normalised_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $normalised_object && \
	bin/scanpy-normalise-data.py -i $filtered_genes_object \
				     -o $normalised_object \
				     -s $ND_scale_factor

    [ "$status" -eq 0 ]
    [ -f  "$normalised_object" ]
}

# Run scanpy-find-variable-genes.R

@test "Find variable genes" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$variable_genes_object" ]; then
        skip "$variable_genes_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $variable_genes_object $variable_image_file && \
	bin/scanpy-find-variable-genes.py -i $normalised_object \
					  -o $variable_genes_object \
					  --flavor $FVG_flavor \
					  -b $FVG_nbins \
					  -p $FVG_parameters \
					  -l $FVG_low_mean,$FVG_low_disp \
					  -j $FVG_high_mean,$FVG_high_disp \
					  -P $variable_image_file

    [ "$status" -eq 0 ]
    [ -f  "$variable_genes_object" ] && [ -f "$variable_image_file" ]
}


# Scale expression values

@test "Scale expression values" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$scaled_object" ]; then
        skip "$scaled_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $scaled_object && \
	bin/scanpy-scale-data.py -i $variable_genes_object \
				 -x $SD_scale_max \
				 -o $scaled_object \
				 -V $SD_vars_to_regress \
				 $SD_zero_center

    [ "$status" -eq 0 ]
    [ -f  "$scaled_object" ]
}

# Run PCA

@test "Run principal component analysis" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$pca_object" ]; then
        skip "$pca_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $pca_object $pca_image_file && \
	bin/scanpy-run-pca.py -i $scaled_object \
			      -o $pca_object \
			      --output-embeddings-file $pca_embeddings_file \
			      --output-loadings-file $pca_loadings_file \
			      --output-stdev-file $pca_stdev_file \
			      --output-var-ratio-file $pca_var_ratio_file \
			      -n $PCA_npcs \
			      --svd-solver $PCA_svd_solver \
			      -s $PCA_random_seed \
			      -P $pca_image_file \
			      --color $PCA_color \
			      --projection $PCA_projection \
			      $PCA_frameon

    [ "$status" -eq 0 ]
    [ -f  "$pca_object" ] && [ -f "$pca_image_file" ] && \
	[ -f "$pca_embeddings_file" ] && [ -f "$pca_loadings_file" ] && \
	[ -f "$pca_stdev_file" ] && [ -f "$pca_var_ratio_file" ]
}

# Run compute graph

@test "Run compute neighbor graph" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$graph_object" ]; then
        skip "$scaled_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $graph_object $graph_image_file && \
	bin/scanpy-neighbours.py -i $pca_object \
				 -o $graph_object \
				 -N $CG_nneighbor \
				 -n $CG_npcs \
				 -s $CG_random_seed \
				 --method $CG_method \
				 $CG_knn

    [ "$status" -eq 0 ]
    [ -f  "$graph_object" ]
}

# Run find cluster

@test "Run find cluster" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$cluster_object" ]; then
        skip "$cluster_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $cluster_object $cluster_text_file && \
	bin/scanpy-find-cluster.py -i $graph_object \
				   -o $cluster_object \
				   --output-text-file $cluster_text_file \
				   --flavor $FC_flavor \
				   --resolution $FC_resolution \
				   --key-added $FC_key_added \
				   -s $FC_random_seed \
				   $FC_use_weight

    [ "$status" -eq 0 ]
    [ -f  "$cluster_object" ] && [ -f "$cluster_text_file" ]
}

# # Run UMAP
# 
# @test "Run UMAP analysis" {
#     if [ "$use_existing_outputs" = 'true' ] && [ -f "$umap_object" ]; then
#         skip "$umap_object exists and use_existing_outputs is set to 'true'"
#     fi
# 
#     run rm -f $umap_object $umap_image_file && \
# 	bin/scanpy-run-umap.py -i $graph_object -o $umap_object \
# 			       -s $UMAP_random_seed \
# 			       -n $UMAP_ncomp \
# 			       --min-dist $UMAP_min_dist \
# 			       --spread $UMAP_spread \
# 			       --alpha $UMAP_alpha \
# 			       --gamma $UMAP_gamma \
# 			       --init-pos $UMAP_initpos \
# 			       -P $umap_image_file \
# 			       --color $UMAP_color \
# 			       --projection $UMAP_projection \
# 			       $UMAP_frameon
# 
#     [ "$status" -eq 0 ]
#     [ -f  "$umap_object" ] && [ -f "$umap_image_file" ]
# }
# 
# # Run TSNE
# 
# @test "Run TSNE analysis" {
#     if [ "$use_existing_outputs" = 'true' ] && [ -f "$tsne_object" ]; then
#         skip "$tsne_object exists and use_existing_outputs is set to 'true'"
#     fi
# 
#     run rm -f $tsne_object $tsne_image_file && \
# 	bin/scanpy-run-tsne.py -i $graph_object -o $tsne_object \
# 			       -s $TSNE_random_seed \
# 			       --perplexity $TSNE_perplexity \
# 			       --early-exaggeration $TSNE_early_exaggeration \
# 			       --learning-rate $TSNE_learning_rate \
# 			       -P $tsne_image_file \
# 			       --color $TSNE_color \
# 			       --projection $TSNE_projection \
# 			       $TSNE_frameon
# 
#     [ "$status" -eq 0 ]
#     [ -f  "$tsne_object" ] && [ -f "$tsne_image_file" ]
# }

# Local Variables:
# mode: sh
# End:
