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
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$input_h5ad" ]; then
        skip "$input_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $input_h5ad && \
	./scanpy-read-10x.py -d $data_dir/ \
			     -o $input_h5ad

    [ "$status" -eq 0 ]
    [ -f  "$input_h5ad" ]
}

@test "Scanpy object creation from 10x (loom)" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$input_loom" ]; then
        skip "$input_loom exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $input_loom && \
	./scanpy-read-10x.py -d $data_dir/ \
			     -o $input_loom

    [ "$status" -eq 0 ]
    [ -f  "$input_loom" ]
}


# Run scanpy-filter-cells.py

@test "Filter cells from a raw object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_cells_h5ad" ]; then
        skip "$filtered_cells_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_cells_h5ad && \
	./scanpy-filter-cells.py -i $input_h5ad \
				 -o $filtered_cells_h5ad \
				 -p $FC_parameters \
				 -l $FC_min_genes \
				 -j $FC_max_genes

    [ "$status" -eq 0 ]
    [ -f  "$filtered_cells_h5ad" ]
}

@test "Filter cells from a raw object (loom)" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_cells_loom" ]; then
        skip "$filtered_cells_loom exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_cells_loom && \
	./scanpy-filter-cells.py -i $input_loom \
				 -o $filtered_cells_loom \
				 -p $FC_parameters \
				 -l $FC_min_genes \
				 -j $FC_max_genes

    [ "$status" -eq 0 ]
    [ -f  "$filtered_cells_loom" ]
}

# Run scanpy-filter-genes.py

@test "Filter genes from a cell-filtered object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_genes_h5ad" ]; then
        skip "$filtered_genes_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_genes_h5ad && \
	./scanpy-filter-genes.py -i $filtered_cells_h5ad \
				 -o $filtered_genes_h5ad \
				 -p $FT_parameters \
				 -l $FT_min_cells \
				 -j $FT_max_cells

    [ "$status" -eq 0 ]
    [ -f  "$filtered_genes_h5ad" ]
}

@test "Filter genes from a cell-filtered object (loom)" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_genes_loom" ]; then
        skip "$filtered_genes_loom exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_genes_loom && \
	./scanpy-filter-genes.py -i $filtered_cells_loom \
				 -o $filtered_genes_loom \
				 -p $FT_parameters \
				 -l $FT_min_cells \
				 -j $FT_max_cells

    [ "$status" -eq 0 ]
    [ -f  "$filtered_genes_loom" ]
}

# Run scanpy-normalise-data.py

@test "Normalise expression values per cell" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$normalised_h5ad" ]; then
        skip "$normalised_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $normalised_h5ad && \
	./scanpy-normalise-data.py -i $filtered_genes_h5ad \
				   -o $normalised_h5ad \
				   -s $ND_scale_factor

    [ "$status" -eq 0 ]
    [ -f  "$normalised_h5ad" ]
}

@test "Normalise expression values per cell (loom)" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$normalised_loom" ]; then
        skip "$normalised_loom exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $normalised_loom && \
	./scanpy-normalise-data.py -i $filtered_genes_loom \
				   -o $normalised_loom \
				   -s $ND_scale_factor

    [ "$status" -eq 0 ]
    [ -f  "$normalised_loom" ]
}

# Run scanpy-find-variable-genes.R

@test "Find variable genes" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$variable_genes_h5ad" ]; then
        skip "$variable_genes_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $variable_genes_h5ad $variable_image_file && \
	./scanpy-find-variable-genes.py -i $normalised_h5ad \
					-o $variable_genes_h5ad \
					--flavor $FVG_flavor \
					-b $FVG_nbins \
					-p $FVG_parameters \
					-l $FVG_low_mean,$FVG_low_disp \
					-j $FVG_high_mean,$FVG_high_disp \
					-P $variable_image_file

    [ "$status" -eq 0 ]
    [ -f  "$variable_genes_h5ad" ] && [ -f "$variable_image_file" ]
}

@test "Find variable genes (loom)" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$variable_genes_loom" ]; then
        skip "$variable_genes_loom exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $variable_genes_loom $variable_image_file && \
	./scanpy-find-variable-genes.py -i $normalised_loom \
					-o $variable_genes_loom \
					--flavor $FVG_flavor \
					-b $FVG_nbins \
					-p $FVG_parameters \
					-l $FVG_low_mean,$FVG_low_disp \
					-j $FVG_high_mean,$FVG_high_disp \
					-P $variable_image_file

    [ "$status" -eq 0 ]
    [ -f  "$variable_genes_loom" ] && [ -f "$variable_image_file" ]
}


# Scale expression values

@test "Scale expression values" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$scaled_h5ad" ]; then
        skip "$scaled_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $scaled_h5ad && \
	./scanpy-scale-data.py -i $variable_genes_h5ad \
			       -x $SD_scale_max \
			       -o $scaled_h5ad \
			       -V $SD_vars_to_regress \
			       $SD_zero_center

    [ "$status" -eq 0 ]
    [ -f  "$scaled_h5ad" ]
}

@test "Scale expression values (loom)" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$scaled_loom" ]; then
        skip "$scaled_loom exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $scaled_loom && \
	./scanpy-scale-data.py -i $variable_genes_loom \
			       -x $SD_scale_max \
			       -o $scaled_loom \
			       -V $SD_vars_to_regress \
			       $SD_zero_center

    [ "$status" -eq 0 ]
    [ -f  "$scaled_loom" ]
}

# Run PCA

@test "Run principal component analysis" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$pca_h5ad" ]; then
        skip "$pca_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $pca_h5ad $pca_image_file && \
	./scanpy-run-pca.py -i $scaled_h5ad \
			    -o $pca_h5ad \
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
    [ -f  "$pca_h5ad" ] && [ -f "$pca_image_file" ] && \
	[ -f "$pca_embeddings_file" ] && [ -f "$pca_loadings_file" ] && \
	[ -f "$pca_stdev_file" ] && [ -f "$pca_var_ratio_file" ]
}

@test "Run principal component analysis (loom)" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$pca_loom" ]; then
        skip "$pca_loom exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $pca_loom $pca_image_file2 && \
	./scanpy-run-pca.py -i $scaled_loom \
			    -o $pca_loom \
			    --output-embeddings-file $pca_embeddings_file2 \
			    --output-loadings-file $pca_loadings_file2 \
			    --output-stdev-file $pca_stdev_file2 \
			    --output-var-ratio-file $pca_var_ratio_file2 \
			    -n $PCA_npcs \
			    --svd-solver $PCA_svd_solver \
			    -s $PCA_random_seed \
			    -P $pca_image_file2 \
			    --color $PCA_color \
			    --projection $PCA_projection \
			    $PCA_frameon

    [ "$status" -eq 0 ]
    [ -f  "$pca_loom" ] && [ -f "$pca_image_file2" ] && \
	[ -f "$pca_embeddings_file2" ] && [ -f "$pca_loadings_file2" ] && \
	[ -f "$pca_stdev_file2" ] && [ -f "$pca_var_ratio_file2" ]
}

# Run compute graph

@test "Run compute neighbor graph" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$graph_h5ad" ]; then
        skip "$scaled_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $graph_h5ad $graph_image_file && \
	./scanpy-neighbours.py -i $pca_h5ad \
			       -o $graph_h5ad \
			       -N $CG_nneighbor \
			       -n $CG_npcs \
			       -s $CG_random_seed \
			       --method $CG_method \
			       $CG_knn

    [ "$status" -eq 0 ]
    [ -f  "$graph_h5ad" ]
}

@test "Run compute neighbor graph (loom)" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$graph_loom" ]; then
        skip "$scaled_loom exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $graph_loom $graph_image_file && \
	./scanpy-neighbours.py -i $pca_h5ad \
			       -o $graph_loom \
			       -N $CG_nneighbor \
			       -n $CG_npcs \
			       -s $CG_random_seed \
			       --method $CG_method \
			       $CG_knn

    [ "$status" -eq 0 ]
    [ -f  "$graph_loom" ]
}

# Run find cluster

@test "Run find cluster" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$cluster_h5ad" ]; then
        skip "$cluster_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $cluster_h5ad $cluster_text_file && \
	./scanpy-find-cluster.py -i $graph_h5ad \
				 -o $cluster_h5ad \
				 --output-text-file $cluster_text_file \
				 --flavor $FC_flavor \
				 --resolution $FC_resolution \
				 --key-added $FC_key_added \
				 -s $FC_random_seed \
				 $FC_use_weight

    [ "$status" -eq 0 ]
    [ -f  "$cluster_h5ad" ] && [ -f "$cluster_text_file" ]
}

@test "Run find cluster (loom)" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$cluster_loom" ]; then
        skip "$cluster_loom exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $cluster_loom $cluster_text_file2 && \
	./scanpy-find-cluster.py -i $graph_loom \
				 -o $cluster_loom \
				 --output-text-file $cluster_text_file2 \
				 --flavor $FC_flavor \
				 --resolution $FC_resolution \
				 --key-added $FC_key_added \
				 -s $FC_random_seed \
				 $FC_use_weight

    [ "$status" -eq 0 ]
    [ -f  "$cluster_loom" ] && [ -f "$cluster_text_file2" ] && [[ -z $(diff -q "$cluster_text_file" "$cluster_text_file2") ]]
}

# Run UMAP

@test "Run UMAP analysis" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$umap_h5ad" ]; then
        skip "$umap_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $umap_h5ad $umap_image_file $umap_embeddings_file && \
	./scanpy-run-umap.py -i $cluster_h5ad -o $umap_h5ad \
			     --output-embeddings-file $umap_embeddings_file \
			     -s $UMAP_random_seed \
			     -n $UMAP_ncomp \
			     --min-dist $UMAP_min_dist \
			     --spread $UMAP_spread \
			     --alpha $UMAP_alpha \
			     --gamma $UMAP_gamma \
			     --init-pos $UMAP_initpos \
			     -P $umap_image_file \
			     --color $UMAP_color \
			     --projection $UMAP_projection \
			     $UMAP_frameon

    [ "$status" -eq 0 ]
    [ -f  "$umap_h5ad" ] && [ -f "$umap_image_file" ] && [ -f "$umap_embeddings_file" ]
}

@test "Run UMAP analysis (loom)" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$umap_loom" ]; then
        skip "$umap_loom exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $umap_loom $umap_image_file2 $umap_embeddings_file && \
	./scanpy-run-umap.py -i $cluster_loom -o $umap_loom \
			     --output-embeddings-file $umap_embeddings_file \
			     -s $UMAP_random_seed \
			     -n $UMAP_ncomp \
			     --min-dist $UMAP_min_dist \
			     --spread $UMAP_spread \
			     --alpha $UMAP_alpha \
			     --gamma $UMAP_gamma \
			     --init-pos $UMAP_initpos \
			     -P $umap_image_file2 \
			     --color $UMAP_color \
			     --projection $UMAP_projection \
			     $UMAP_frameon

    [ "$status" -eq 0 ]
    [ -f  "$umap_loom" ] && [ -f "$umap_image_file2" ] && [ -f "$umap_embeddings_file" ] && [[ -z $(diff -q "$umap_image_file" "umap_image_file2") ]]
}

# Run TSNE

@test "Run TSNE analysis" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$tsne_h5ad" ]; then
        skip "$tsne_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $tsne_h5ad $tsne_image_file $tsne_embeddings_file && \
	./scanpy-run-tsne.py -i $cluster_h5ad -o $tsne_h5ad \
			     --output-embeddings-file $tsne_embeddings_file \
			     -s $TSNE_random_seed \
			     --perplexity $TSNE_perplexity \
			     --early-exaggeration $TSNE_early_exaggeration \
			     --learning-rate $TSNE_learning_rate \
			     -P $tsne_image_file \
			     --color $TSNE_color \
			     --projection $TSNE_projection \
			     $TSNE_frameon

    [ "$status" -eq 0 ]
    [ -f  "$tsne_h5ad" ] && [ -f "$tsne_image_file" ] && [ -f "$tsne_embeddings_file" ]
}

@test "Run TSNE analysis (loom)" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$tsne_loom" ]; then
        skip "$tsne_loom exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $tsne_loom $tsne_image_file2 $tsne_embeddings_file && \
	./scanpy-run-tsne.py -i $cluster_loom -o $tsne_loom \
			     --output-embeddings-file $tsne_embeddings_file \
			     -s $TSNE_random_seed \
			     --perplexity $TSNE_perplexity \
			     --early-exaggeration $TSNE_early_exaggeration \
			     --learning-rate $TSNE_learning_rate \
			     -P $tsne_image_file2 \
			     --color $TSNE_color \
			     --projection $TSNE_projection \
			     $TSNE_frameon

    [ "$status" -eq 0 ]
    [ -f  "$tsne_loom" ] && [ -f "$tsne_image_file2" ] && [ -f "$tsne_embeddings_file" ] && [[ -z $(diff -q "$tsne_image_file" "tsne_image_file2") ]]
}

# Run find markers

@test "Run find markers" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$marker_h5ad" ]; then
        skip "$marker_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $marker_h5ad $marker_image_file $marker_text_file && \
	./scanpy-find-markers.py -i $cluster_h5ad -o $marker_h5ad \
				 --output-text-file $marker_text_file \
				 --groupby $FM_groupby \
				 --groups $FM_groups \
				 --reference $FM_reference \
				 --n-genes $FM_n_genes \
				 --method $FM_method \
				 -P $marker_image_file \
				 --show-n-genes $FM_show_n_genes \
				 --debug \
				 --key $FM_key

    [ "$status" -eq 0 ]
    [ -f  "$marker_h5ad" ] && [ -f "$marker_image_file" ] && [ -f "$marker_text_file" ]
}

@test "Run find markers (loom)" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$marker2_h5ad" ]; then
        skip "$marker2_h5ad exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $marker_loom $marker_image_file2 $marker_text_file2 && \
	./scanpy-find-markers.py -i $cluster_loom -o $marker2_h5ad \
				 --output-text-file $marker_text_file2 \
				 --groupby $FM_groupby \
				 --groups $FM_groups \
				 --reference $FM_reference \
				 --n-genes $FM_n_genes \
				 --method $FM_method \
				 -P $marker_image_file2 \
				 --show-n-genes $FM_show_n_genes \
				 --debug \
				 --key $FM_key

    [ "$status" -eq 0 ]
    [ -f  "$marker_h5ad" ] && [ -f "$marker_image_file2" ] && [ -f "$marker_text_file2" ] && [[ -z $(diff -q "$marker_text_file" "$marker_text_file2") ]]
}

# Local Variables:
# mode: sh
# End:
