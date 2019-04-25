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
    norm_mtx="${output_dir}/norm"
    norm_opt="-r yes -t 10000 -X ${norm_mtx} --show-obj stdout"
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
    tsne_opt="-n 25 --use-rep X_pca --learning-rate 200"
    tsne_obj="${output_dir}/tsne.h5ad"
    umap_opt="--use-graph neighbors_k10 --min-dist 0.75 --alpha 1 --gamma 1"
    umap_obj="${output_dir}/umap.h5ad"
    louvain_opt="-r 1 --use-graph neighbors_k10 --key-added louvain_k10_r1_0"
    louvain_obj="${output_dir}/louvain.h5ad"
    leiden_opt="-r 0.7 --use-graph neighbors_k10 --key-added leiden_k10_r0_7"
    leiden_obj="${output_dir}/leiden.h5ad"
    diffexp_tsv="${output_dir}/diffexp.tsv"
    diffexp_opt="-g leiden_k10_r0_7 --reference rest --save ${diffexp_tsv}"
    diffexp_obj="${output_dir}/diffexp.h5ad"
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
    [ -f  "$norm_obj" ] && [ -f "${norm_mtx}_matrix.mtx" ]
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

# Run TSNE

@test "Run TSNE analysis" {
    if [ "$resume" = 'true' ] && [ -f "$tsne_obj" ]; then
        skip "$tsne_obj exists and resume is set to 'true'"
    fi

    run rm -f $tsne_obj && $scanpy embed tsne $tsne_opt $pca_obj $tsne_obj

    [ "$status" -eq 0 ]
    [ -f  "$tsne_obj" ]
}

# Run UMAP

@test "Run UMAP analysis" {
    if [ "$resume" = 'true' ] && [ -f "$umap_obj" ]; then
        skip "$umap_obj exists and resume is set to 'true'"
    fi

    run rm -f $umap_obj && $scanpy embed umap $umap_opt $neighbor_obj $umap_obj

    [ "$status" -eq 0 ]
    [ -f  "$umap_obj" ]
}

# Find clusters Louvain

@test "Run find cluster (louvain)" {
    if [ "$resume" = 'true' ] && [ -f "$louvain_obj" ]; then
        skip "$louvain_obj exists and resume is set to 'true'"
    fi

    run rm -f $louvain_obj && $scanpy cluster louvain $louvain_opt $umap_obj $louvain_obj

    [ "$status" -eq 0 ]
    [ -f  "$louvain_obj" ]
}

# Find clusters Leiden

@test "Run find cluster (leiden)" {
    if [ "$resume" = 'true' ] && [ -f "$leiden_obj" ]; then
        skip "$leiden_obj exists and resume is set to 'true'"
    fi

    run rm -f $leiden_obj && $scanpy cluster leiden $leiden_opt $umap_obj $leiden_obj

    [ "$status" -eq 0 ]
    [ -f  "$leiden_obj" ]
}

# Find markers

@test "Run find markers" {
    if [ "$resume" = 'true' ] && [ -f "$diffexp_obj" ]; then
        skip "$diffexp_obj exists and resume is set to 'true'"
    fi

    run rm -f $diffexp_obj $diffexp_tsv && $scanpy diffexp $diffexp_opt $leiden_obj $diffexp_obj

    [ "$status" -eq 0 ]
    [ -f  "$diffexp_obj" ] && [ -f "$diffexp_tsv" ]
}

# Local Variables:
# mode: sh
# End:
