#!/usr/bin/env bats

# Extract the test data
setup() {
    scanpy="scanpy-cli"
    test_dir="post_install_tests"
    data_dir="${test_dir}/data"
    output_dir="${test_dir}/outputs"
    test_data_url='https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
    test_data_archive="${test_dir}/$(basename $test_data_url)"
    raw_matrix="${data_dir}/matrix.mtx"
    read_opt="-x $data_dir --show-obj stdout"
    read_obj="${output_dir}/read.h5ad"
    filter_opt="-p n_genes 200 2500 -p cell:n_counts 0 50000 -p n_cells 3 inf -p pct_counts_mito 0 0.2 --show-obj stdout"
    filter_obj="${output_dir}/filter.h5ad"
    norm_opt="-r yes -t 10000 --show-obj stdout"
    norm_obj="${output_dir}/norm.h5ad"
    hvg_opt="-m 0.0125 3 -d 0.5 inf -s --show-obj stdout"
    hvg_obj="${output_dir}/hvg.h5ad"
    regress_opt="-k n_counts --show-obj stdout"
    regress_obj="${output_dir}/regress.h5ad"
    scale_opt="-m 10 --show-obj stdout"
    scale_obj="${output_dir}/scale.h5ad"
    pca_opt="-n 50 -V arpack --show-obj stdout"
    pca_obj="${output_dir}/pca.h5ad"
    neighbor_opt="-k 5,10,20 -n 25 -m umap --show-obj stdout"
    neighbor_obj="${output_dir}/neighbor.h5ad"
}

@test "Downlaod and extract .mtx matrix" {
    if [ -f "$raw_matrix" ]; then
        skip "$raw_matrix exists"
    fi

    run wget $test_data_url -P $test_dir && tar -xvzf $test_data_archive --strip-components 2 -C $data_dir

    [ "$status" -eq 0 ]
    [ -f "$raw_matrix" ]
}

# Read 10x dataset

@test "Scanpy object creation from 10x" {
    if [ "$resume" = 'true' ] && [ -f "$read_obj" ]; then
        skip "$read_obj exists and resume is set to 'true'"
    fi

    run rm -f $read_obj && $scanpy read $read_opt $read_obj

    [ "$status" -eq 0 ]
    [ -f  "$read_obj" ]
}

# Filter

@test "Filter cells and genes from a raw object" {
    if [ "$resume" = 'true' ] && [ -f "$filter_obj" ]; then
        skip "$filter_obj exists and resume is set to 'true'"
    fi

    run rm -f $filter_obj && $scanpy filter $filter_opt $read_obj $filter_obj

    [ "$status" -eq 0 ]
    [ -f  "$filter_obj" ]
}

# Normalise

@test "Normalise expression values per cell" {
    if [ "$resume" = 'true' ] && [ -f "$norm_obj" ]; then
        skip "$norm_obj exists and resume is set to 'true'"
    fi

    run rm -f $norm_obj && $scanpy norm $norm_opt $filter_obj $norm_obj

    [ "$status" -eq 0 ]
    [ -f  "$norm_obj" ]
}

# Find variable genes

@test "Find variable genes" {
    if [ "$resume" = 'true' ] && [ -f "$hvg_obj" ]; then
        skip "$hvg_obj exists and resume is set to 'true'"
    fi

    run rm -f $hvg_obj $hvg_obj && $scanpy hvg $hvg_opt $norm_obj $hvg_obj

    [ "$status" -eq 0 ]
    [ -f  "$hvg_obj" ]
}

# Regress out variables

@test "Regress out unwanted variable" {
    if [ "$resume" = 'true' ] && [ -f "$regress_obj" ]; then
        skip "$regress_obj exists and resume is set to 'true'"
    fi

    run rm -f $regress_obj && $scanpy regress $regress_opt $hvg_obj $regress_obj

    [ "$status" -eq 0 ]
    [ -f  "$regress_obj" ]
}

# Scale expression values

@test "Scale expression values" {
    if [ "$resume" = 'true' ] && [ -f "$scale_obj" ]; then
        skip "$scale_obj exists and resume is set to 'true'"
    fi

    run rm -f $scale_obj && $scanpy scale $scale_opt $hvg_obj $scale_obj

    [ "$status" -eq 0 ]
    [ -f  "$scale_obj" ]
}

# Run PCA

@test "Run principal component analysis" {
    if [ "$resume" = 'true' ] && [ -f "$pca_obj" ]; then
        skip "$pca_obj exists and resume is set to 'true'"
    fi

    run rm -f $pca_obj && $scanpy pca $pca_opt $scale_obj $pca_obj

    [ "$status" -eq 0 ]
    [ -f  "$pca_obj" ]
}

# Compute graph

@test "Run compute neighbor graph" {
    if [ "$resume" = 'true' ] && [ -f "$neighbor_obj" ]; then
        skip "$scaled_object exists and resume is set to 'true'"
    fi

    run rm -f $neighbor_obj && $scanpy neighbor $neighbor_opt $pca_obj $neighbor_obj

    [ "$status" -eq 0 ]
    [ -f  "$neighbor_obj" ]
}

# # Find clusters

# @test "Run find cluster" {
#     if [ "$resume" = 'true' ] && [ -f "$cluster_object" ]; then
#         skip "$cluster_object exists and resume is set to 'true'"
#     fi

#     run rm -f $cluster_object $cluster_text_file && \
#         scanpy-find-cluster -i $graph_object \
#             -o $cluster_object \
#             --output-text-file $cluster_text_file \
#             --flavor $FC_flavor \
#             --resolution $FC_resolution \
#             --key-added $FC_key_added \
#             -s $FC_random_seed \
#             $FC_use_weight

#     [ "$status" -eq 0 ]
#     [ -f  "$cluster_object" ] && [ -f "$cluster_text_file" ]
# }

# # Run UMAP

# @test "Run UMAP analysis" {
#     if [ "$resume" = 'true' ] && [ -f "$umap_object" ]; then
#         skip "$umap_object exists and resume is set to 'true'"
#     fi

#     run rm -f $umap_object $umap_image_file $umap_embeddings_file && \
#         scanpy-run-umap -i $cluster_object -o $umap_object \
#             --output-embeddings-file $umap_embeddings_file \
#             -s $UMAP_random_seed \
#             -n $UMAP_ncomp \
#             --min-dist $UMAP_min_dist \
#             --spread $UMAP_spread \
#             --alpha $UMAP_alpha \
#             --gamma $UMAP_gamma \
#             --init-pos $UMAP_initpos \
#             -P $umap_image_file \
#             --color $UMAP_color \
#             --projection $UMAP_projection \
#             $UMAP_frameon

#     [ "$status" -eq 0 ]
#     [ -f  "$umap_object" ] && [ -f "$umap_image_file" ] && [ -f "$umap_embeddings_file" ]
# }

# # Run TSNE

# @test "Run TSNE analysis" {
#     if [ "$resume" = 'true' ] && [ -f "$tsne_object" ]; then
#         skip "$tsne_object exists and resume is set to 'true'"
#     fi

#     run rm -f $tsne_object $tsne_image_file $tsne_embeddings_file && \
#         scanpy-run-tsne -i $cluster_object -o $tsne_object \
#             --output-embeddings-file $tsne_embeddings_file \
#             -s $TSNE_random_seed \
#             --perplexity $TSNE_perplexity \
#             --early-exaggeration $TSNE_early_exaggeration \
#             --learning-rate $TSNE_learning_rate \
#             -P $tsne_image_file \
#             --color $TSNE_color \
#             --projection $TSNE_projection \
#             $TSNE_frameon

#     [ "$status" -eq 0 ]
#     [ -f  "$tsne_object" ] && [ -f "$tsne_image_file" ] && [ -f "$tsne_embeddings_file" ]
# }

# # Find markers

# @test "Run find markers" {
#     if [ "$resume" = 'true' ] && [ -f "$marker_object" ]; then
#         skip "$marker_object exists and resume is set to 'true'"
#     fi

#     run rm -f $marker_object $marker_image_file $marker_text_file && \
#         scanpy-find-markers -i $cluster_object -o $marker_object \
#             --output-text-file $marker_text_file \
#             --groupby $FM_groupby \
#             --groups $FM_groups \
#             --reference $FM_reference \
#             --n-genes $FM_n_genes \
#             --method $FM_method \
#             -P $marker_image_file \
#             --show-n-genes $FM_show_n_genes \
#             --debug \
#             --key $FM_key

#     [ "$status" -eq 0 ]
#     [ -f  "$marker_object" ] && [ -f "$marker_image_file" ] && [ -f "$marker_text_file" ]
# }

# Local Variables:
# mode: sh
# End:
