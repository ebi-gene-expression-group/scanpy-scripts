#!/usr/bin/env bats

# Extract the test data
setup() {
    scanpy="scanpy-cli"
    test_dir="post_install_tests"
    data_dir="${test_dir}/data"
    output_dir="${test_dir}/outputs"
    raw_matrix="${data_dir}/matrix.mtx"
    singlet_obs="${data_dir}/singlet_obs.txt"
    read_opt="-x $data_dir --show-obj stdout"
    read_obj="${output_dir}/read.h5ad"
    filter_opt="--save-raw -p n_genes 200 2500 -p c:n_counts 0 50000 -p n_cells 3 inf -p pct_counts_mito 0 0.2 -c mito '!True' --show-obj stdout"
    filter_obj="${output_dir}/filter.h5ad"
    scrublet_tsv="${output_dir}/scrublet.tsv"
    scrublet_png="${output_dir}/scrublet.png"
    scrublet_obj="${output_dir}/scrublet.h5ad"
    scrublet_simulate_obj="${output_dir}/scrublet_simulate.h5ad"
    scrublet_opt="--input-obj-sim ${scrublet_simulate_obj} --filter --export-table ${scrublet_tsv}"
    norm_mtx="${output_dir}/norm"
    norm_opt="--save-layer filtered -t 10000 -l all -n after -X ${norm_mtx} --show-obj stdout"
    norm_obj="${output_dir}/norm.h5ad"
    hvg_opt="-m 0.0125 3 -d 0.5 inf -s --show-obj stdout"
    hvg_obj="${output_dir}/hvg.h5ad"
    regress_opt="-k n_counts --show-obj stdout"
    regress_obj="${output_dir}/regress.h5ad"
    scale_opt="--save-layer normalised -m 10 --show-obj stdout"
    scale_obj="${output_dir}/scale.h5ad"
    pca_embed="${output_dir}/pca.tsv"
    pca_opt="--n-comps 50 -V auto --show-obj stdout -E ${pca_embed}"
    pca_obj="${output_dir}/pca.h5ad"
    neighbor_opt="-k 5,10,20 -n 25 -m umap -t euclidean --show-obj stdout"
    neighbor_obj="${output_dir}/neighbor.h5ad"
    tsne_embed="${output_dir}/tsne.tsv"
    tsne_opt="-n 25 --use-rep X_pca --learning-rate 200 -E ${tsne_embed}"
    tsne_obj="${output_dir}/tsne.h5ad"
    umap_embed="${output_dir}/umap.tsv"
    umap_opt="--neighbors-key k10 --min-dist 0.75 --alpha 1 --gamma 1 -E ${umap_embed}"
    umap_obj="${output_dir}/umap.h5ad"
    fdg_embed="${output_dir}/fdg.tsv"
    fdg_opt="--neighbors-key k10 --layout fr -E ${fdg_embed} --init-pos paga"
    fdg_obj="${output_dir}/fdg.h5ad"
    louvain_tsv="${output_dir}/louvain.tsv"
    louvain_opt="-r 0.5,1 --neighbors-key k10 --key-added k10 --export-cluster ${louvain_tsv}"
    louvain_obj="${output_dir}/louvain.h5ad"
    leiden_tsv="${output_dir}/leiden.tsv"
    leiden_opt="-r 0.3,0.7 --neighbors-key k10 --key-added k10 -F loom --loom-write-obsm-varm --export-cluster ${leiden_tsv}"
    leiden_obj="${output_dir}/leiden.loom"
    test_clustering='louvain_k10_r0_5'
    diffexp_tsv="${output_dir}/diffexp.tsv"
    diffexp_opt="-g ${test_clustering} --reference rest --filter-params min_in_group_fraction:0.25,min_fold_change:1.5 --save ${diffexp_tsv}"
    diffexp_obj="${output_dir}/diffexp.h5ad"
    test_singlet_clustering='groupby_with_singlet'
    diffexp_singlet_tsv="${output_dir}/diffexp_singlet.tsv"
    diffexp_singlet_opt="-g ${test_singlet_clustering} --reference rest --filter-params min_in_group_fraction:0.25,min_fold_change:1.5 --save ${diffexp_singlet_tsv}"
    diffexp_singlet_obj="${output_dir}/diffexp_singlet.h5ad"
    paga_opt="--neighbors-key k10 --key-added ${test_clustering} --groups ${test_clustering} --model v1.2"
    paga_obj="${output_dir}/paga.h5ad"
    diffmap_embed="${output_dir}/diffmap.tsv"
    diffmap_opt="--neighbors-key k10 --n-comps 10 -E ${diffmap_embed}"
    diffmap_obj="${output_dir}/diffmap.h5ad"
    dpt_opt="--neighbors-key k10 --key-added k10 --n-dcs 10 --disallow-kendall-tau-shift --root ${test_clustering} 0"
    dpt_obj="${output_dir}/dpt.h5ad"
    plt_embed_opt="--projection 2d --color ${test_clustering} --title test"
    plt_embed_pdf="${output_dir}/umap_${test_clustering}.pdf"
    plt_paga_pdf="${output_dir}/paga_k10_r0_7.pdf"
    plt_paga_obj="${output_dir}/paga_k10_r0_7.h5ad"
    plt_paga_opt="--use-key paga_${test_clustering} --node-size-scale 2 --edge-width-scale 0.5 --basis diffmap --color dpt_pseudotime_k10 --frameoff --output-obj $plt_paga_obj"
    test_markers='LDHB,CD3D,CD3E'
    diffexp_plot_opt="--var-names $test_markers --use-raw --dendrogram --groupby ${test_clustering}"
    plt_stacked_violin_opt="${diffexp_plot_opt} --no-jitter --swap-axes"
    plt_stacked_violin_pdf="${output_dir}/sviolin_${test_clustering}_LDHB_CD3D_CD3E.pdf"
    plt_dotplot_pdf="${output_dir}/dot_${test_clustering}_LDHB_CD3D_CD3E.pdf"
    plt_matrixplot_pdf="${output_dir}/matrix_${test_clustering}_LDHB_CD3D_CD3E.pdf"
    plt_heatmap_pdf="${output_dir}/heatmap_${test_clustering}_LDHB_CD3D_CD3E.pdf"
    plt_rank_genes_groups_opt="--rgg --groups 3,4"
    plt_rank_genes_groups_singlet_opt="--rgg"
    plt_rank_genes_groups_stacked_violin_pdf="${output_dir}/rggsviolin_${test_clustering}.pdf"
    plt_rank_genes_groups_matrix_pdf="${output_dir}/rggmatrix_${test_clustering}.pdf"
    plt_rank_genes_groups_dot_pdf="${output_dir}/rggdot_${test_clustering}.pdf"
    plt_rank_genes_groups_dot_singlet_pdf="${output_dir}/rggdot_${test_singlet_clustering}.pdf"
    plt_rank_genes_groups_heatmap_pdf="${output_dir}/rggheatmap_${test_clustering}.pdf"
    harmony_integrate_obj="${output_dir}/harmony_integrate.h5ad"
    harmony_integrate_opt="--batch-key ${test_clustering}"
    harmony_plt_embed_opt="--projection 2d --color ${test_clustering} --title 'PCA embeddings after harmony' --basis 'X_pca_harmony'"
    noharmony_plt_embed_opt="--projection 2d --color ${test_clustering} --title 'PCA embeddings before harmony' --basis 'X_pca'"
    harmony_integrated_pca_pdf="${output_dir}/harmony_pca_${test_clustering}.pdf"
    noharmony_integrated_pca_pdf="${output_dir}/pca_${test_clustering}.pdf"
    bbknn_obj="${output_dir}/bbknn.h5ad"
    bbknn_opt="--batch-key ${test_clustering} --key-added bbknn"
    mnn_obj="${output_dir}/mnn.h5ad"
    mnn_opt="--save-layer uncorrected --batch-key ${test_clustering}"
    combat_obj="${output_dir}/combat.h5ad"
    combat_opt="--batch-key ${test_clustering}"
    

    if [ ! -d "$data_dir" ]; then
        mkdir -p $data_dir
    fi

    if [ ! -d "$output_dir" ]; then
        mkdir -p $output_dir
    fi
}

@test "Extract test data from Scanpy" {
    if [ "$resume" = 'true' ] && [ -f "$raw_matrix" ]; then
        skip "$raw_matrix exists"
    fi

    run rm -rf ${data_dir}/* && eval "echo -e \"import scanpy as sc\nfrom scanpy_scripts.cmd_utils import write_mtx\nimport os\nos.makedirs('$data_dir', exist_ok=True)\nwrite_mtx(sc.datasets.pbmc3k(), '$data_dir/')\" | python"

    [ "$status" -eq 0 ]
    [ -f "$raw_matrix" ]
}

@test "Make .obs with a singlet cell group" {

    if [ "$resume" = 'true' ] && [ -f "$singlet_obs" ]; then
        skip "$singlet_obs exists"
    fi

    run rm -rf $singlet_obs && eval "echo -e \"index\tgroupby_with_singlet\" > $singlet_obs && head -n 1 $data_dir/barcodes.tsv | awk -v cluster='cluster1' '{print \$1\"\t\"cluster}' >> $singlet_obs && sed -n '2,100p;101q' $data_dir/barcodes.tsv | awk -v cluster='cluster3' '{print \$1\"\t\"cluster}' >> $singlet_obs && tail -n +101 $data_dir/barcodes.tsv | awk -v cluster='cluster2' '{print \$1\"\t\"cluster}' >> $singlet_obs"

    [ "$status" -eq 0 ]
    [ -f "$singlet_obs" ]
}

# Read 10x dataset

@test "Scanpy object creation from 10x" {
    if [ "$resume" = 'true' ] && [ -f "$read_obj" ]; then
        skip "$read_obj exists and resume is set to 'true'"
    fi

    run rm -f $read_obj && eval "$scanpy read --extra-obs $singlet_obs $read_opt $read_obj"

    [ "$status" -eq 0 ]
    [ -f  "$read_obj" ]
}

# Filter

@test "Filter cells and genes from a raw object" {
    if [ "$resume" = 'true' ] && [ -f "$filter_obj" ]; then
        skip "$filter_obj exists and resume is set to 'true'"
    fi

    run rm -f $filter_obj && eval "$scanpy filter $filter_opt $read_obj $filter_obj"

    [ "$status" -eq 0 ]
    [ -f  "$filter_obj" ]
}

# Normalise

@test "Normalise expression values per cell" {
    if [ "$resume" = 'true' ] && [ -f "$norm_obj" ]; then
        skip "$norm_obj exists and resume is set to 'true'"
    fi

    run rm -f $norm_obj && eval "$scanpy norm $norm_opt $filter_obj $norm_obj"

    [ "$status" -eq 0 ]
    [ -f  "$norm_obj" ] && [ -f "${norm_mtx}_matrix.mtx" ]
}

# Find variable genes

@test "Find variable genes" {
    if [ "$resume" = 'true' ] && [ -f "$hvg_obj" ]; then
        skip "$hvg_obj exists and resume is set to 'true'"
    fi

    run rm -f $hvg_obj $hvg_obj && eval "$scanpy hvg $hvg_opt $norm_obj $hvg_obj"

    [ "$status" -eq 0 ]
    [ -f  "$hvg_obj" ]
}

# Do separate doublet simulation step (normally we'd just let the main scrublet
# process do this).

@test "Run Scrublet doublet simulation" {
    if [ "$resume" = 'true' ] && [ -f "$scrublet_simulate_obj" ]; then
        skip "$scrublet_simulate_obj exists and resume is set to 'true'"
    fi

    run rm -f $srublet_simulate_obj && eval "$scanpy multiplet scrublet_simulate_doublets $hvg_obj $scrublet_simulate_obj"

    [ "$status" -eq 0 ]
    [ -f  "$scrublet_simulate_obj" ]
}

# Detect multiplets with Scrublet

@test "Run Scrublet for multiplet detection" {
    if [ "$resume" = 'true' ] && [ -f "$scrublet_obj" ]; then
        skip "$scrublet_obj exists and resume is set to 'true'"
    fi

    run rm -f $scrublet_obj && eval "$scanpy multiplet scrublet $scrublet_opt $hvg_obj $scrublet_obj"

    [ "$status" -eq 0 ]
    [ -f  "$scrublet_obj" ] && [ -f "$scrublet_tsv" ]
}

# Run the doublet plot from Scrublet

@test "Run Scrublet score distribution plot" {
    if [ "$resume" = 'true' ] && [ -f "$scrublet_png" ]; then
        skip "$scrublet_png exists and resume is set to 'true'"
    fi

    run rm -f $scrublet_png && eval "$scanpy plot scrublet $scrublet_obj $scrublet_png"

    [ "$status" -eq 0 ]
    [ -f  "$scrublet_png" ]
}

# Regress out variables

@test "Regress out unwanted variable" {
    if [ "$resume" = 'true' ] && [ -f "$regress_obj" ]; then
        skip "$regress_obj exists and resume is set to 'true'"
    fi

    run rm -f $regress_obj && eval "$scanpy regress $regress_opt $hvg_obj $regress_obj"

    [ "$status" -eq 0 ]
    [ -f  "$regress_obj" ]
}

# Scale expression values

@test "Scale expression values" {
    if [ "$resume" = 'true' ] && [ -f "$scale_obj" ]; then
        skip "$scale_obj exists and resume is set to 'true'"
    fi

    run rm -f $scale_obj && eval "$scanpy scale $scale_opt $hvg_obj $scale_obj"

    [ "$status" -eq 0 ]
    [ -f  "$scale_obj" ]
}

# Run PCA

@test "Run principal component analysis" {
    if [ "$resume" = 'true' ] && [ -f "$pca_obj" ]; then
        skip "$pca_obj exists and resume is set to 'true'"
    fi

    run rm -f $pca_obj && eval "$scanpy pca $pca_opt $scale_obj $pca_obj"

    [ "$status" -eq 0 ]
    [ -f  "$pca_obj" ]
}

# Compute graph

@test "Run compute neighbor graph" {
    if [ "$resume" = 'true' ] && [ -f "$neighbor_obj" ]; then
        skip "$scaled_object exists and resume is set to 'true'"
    fi

    run rm -f $neighbor_obj && eval "$scanpy neighbor $neighbor_opt $pca_obj $neighbor_obj"

    [ "$status" -eq 0 ]
    [ -f  "$neighbor_obj" ]
}

# Run TSNE

@test "Run TSNE analysis" {
    if [ "$resume" = 'true' ] && [ -f "$tsne_obj" ]; then
        skip "$tsne_obj exists and resume is set to 'true'"
    fi

    run rm -f $tsne_obj && eval "$scanpy embed tsne $tsne_opt $pca_obj $tsne_obj"

    [ "$status" -eq 0 ]
    [ -f  "$tsne_obj" ] && [ -f "$tsne_embed" ]
}

# Run UMAP

@test "Run UMAP analysis" {
    if [ "$resume" = 'true' ] && [ -f "$umap_obj" ]; then
        skip "$umap_obj exists and resume is set to 'true'"
    fi

    run rm -f $umap_obj && eval "$scanpy embed umap $umap_opt $neighbor_obj $umap_obj"

    [ "$status" -eq 0 ]
    [ -f  "$umap_obj" ] && [ -f "$umap_embed" ]
}

# Find clusters Louvain

@test "Run find cluster (louvain)" {
    if [ "$resume" = 'true' ] && [ -f "$louvain_obj" ]; then
        skip "$louvain_obj exists and resume is set to 'true'"
    fi

    run rm -f $louvain_obj && eval "$scanpy cluster louvain $louvain_opt $umap_obj $louvain_obj"

    [ "$status" -eq 0 ]
    [ -f  "$louvain_obj" ] && [ -f "$louvain_tsv" ]
}

# Find clusters Leiden

@test "Run find cluster (leiden)" {
    if [ "$resume" = 'true' ] && [ -f "$leiden_obj" ]; then
        skip "$leiden_obj exists and resume is set to 'true'"
    fi

    run rm -f $leiden_obj && eval "$scanpy cluster leiden $leiden_opt $umap_obj $leiden_obj"

    [ "$status" -eq 0 ]
    [ -f  "$leiden_obj" ] && [ -f "$leiden_tsv" ]
}

# Find markers

@test "Run find markers" {
    if [ "$resume" = 'true' ] && [ -f "$diffexp_obj" ]; then
        skip "$diffexp_obj exists and resume is set to 'true'"
    fi

    run rm -f $diffexp_obj $diffexp_tsv && eval "$scanpy diffexp $diffexp_opt $louvain_obj $diffexp_obj"

    [ "$status" -eq 0 ]
    [ -f  "$diffexp_obj" ] && [ -f "$diffexp_tsv" ]
}

# Find markers, with singlet group

@test "Run find markers, with singlet group ignored" {
    if [ "$resume" = 'true' ] && [ -f "$diffexp_singlet_obj" ]; then
        skip "$diffexp_singlet_obj exists and resume is set to 'true'"
    fi

    run rm -f $diffexp_singlet_obj $diffexp_singlet_tsv && eval "$scanpy diffexp $diffexp_singlet_opt $louvain_obj $diffexp_singlet_obj"

    [ "$status" -eq 0 ]
    [ -f  "$diffexp_singlet_obj" ] && [ -f "$diffexp_singlet_tsv" ]
}

# Run PAGA

@test "Run PAGA" {
    if [ "$resume" = 'true' ] && [ -f "$paga_obj" ]; then
        skip "$paga_obj exists and resume is set to 'true'"
    fi

    run rm -f $paga_obj && eval "$scanpy paga $paga_opt $louvain_obj $paga_obj"

    [ "$status" -eq 0 ]
    [ -f  "$paga_obj" ]
}

# Run Diffmap

@test "Run Diffmap" {
    if [ "$resume" = 'true' ] && [ -f "$diffmap_obj" ]; then
        skip "$diffmap_obj exists and resume is set to 'true'"
    fi

    run rm -f $diffmap_obj && eval "$scanpy embed diffmap $diffmap_opt $paga_obj $diffmap_obj"

    [ "$status" -eq 0 ]
    [ -f  "$diffmap_obj" ] && [ -f "$diffmap_embed" ]
}

# Run DPT

@test "Run DPT" {
    if [ "$resume" = 'true' ] && [ -f "$dpt_obj" ]; then
        skip "$dpt_obj exists and resume is set to 'true'"
    fi

    run rm -f $dpt_obj && eval "$scanpy dpt $dpt_opt $diffmap_obj $dpt_obj"

    [ "$status" -eq 0 ]
    [ -f  "$dpt_obj" ]
}

# Run Plot embedding

@test "Run Plot embedding" {
    if [ "$resume" = 'true' ] && [ -f "$plt_embed_pdf" ]; then
        skip "$plt_embed_pdf exists and resume is set to 'true'"
    fi

    run rm -f $plt_embed_pdf && eval "$scanpy plot embed $plt_embed_opt $louvain_obj $plt_embed_pdf"

    [ "$status" -eq 0 ]
    [ -f  "$plt_embed_pdf" ]
}

# Run Plot paga

@test "Run Plot trajectory" {
    if [ "$resume" = 'true' ] && [ -f "$plt_paga_pdf" ]; then
        skip "$plt_paga_pdf exists and resume is set to 'true'"
    fi

    run rm -f $plt_paga_pdf && eval "$scanpy plot paga $plt_paga_opt $dpt_obj $plt_paga_pdf"

    [ "$status" -eq 0 ]
    [ -f  "$plt_paga_pdf" ] && [ -f  "$plt_paga_obj" ]
}

# Run FDG, with initial coordinates from paga plotting

@test "Run FDG analysis" {
    if [ "$resume" = 'true' ] && [ -f "$fdg_obj" ]; then
        skip "$fdg_obj exists and resume is set to 'true'"
    fi

    run rm -f $fdg_obj && eval "$scanpy embed fdg $fdg_opt $plt_paga_obj $fdg_obj"

    [ "$status" -eq 0 ]
    [ -f  "$fdg_obj" ] && [ -f "$fdg_embed" ]
}


# Plot a stacked violin plot for markers

@test "Run Plot stacked violins" {
    if [ "$resume" = 'true' ] && [ -f "$plt_stacked_violin_pdf" ]; then
        skip "$plt_stacked_violin_pdf exists and resume is set to 'true'"
    fi

    run rm -f $plt_stacked_violin_pdf && eval "$scanpy plot sviol $plt_stacked_violin_opt $diffexp_obj $plt_stacked_violin_pdf"

    [ "$status" -eq 0 ]
    [ -f  "$plt_stacked_violin_pdf" ]
}

# Plot ranking of genes using a stacked violin plot for markers

@test "Run Plot ranking of genes using stacked_violin plot" {
    if [ "$resume" = 'true' ] && [ -f "$plt_rank_genes_groups_stacked_violin_pdf" ]; then
        skip "$plt_rank_genes_groups_stacked_violin_pdf exists and resume is set to 'true'"
    fi

    run rm -f $plt_rank_genes_groups_stacked_violin_pdf && eval "$scanpy plot sviol $plt_rank_genes_groups_opt $diffexp_obj $plt_rank_genes_groups_stacked_violin_pdf"

    [ "$status" -eq 0 ]
    [ -f  "$plt_rank_genes_groups_stacked_violin_pdf" ]
}

# Plot a dot plot for markers

@test "Run Plot dotplot" {
    if [ "$resume" = 'true' ] && [ -f "$plt_dotplot_pdf" ]; then
        skip "$plt_dotplot_pdf exists and resume is set to 'true'"
    fi

    run rm -f $plt_dotplot_pdf && eval "$scanpy plot dot $diffexp_plot_opt $diffexp_obj $plt_dotplot_pdf"
    
    [ "$status" -eq 0 ]
    [ -f  "$plt_dotplot_pdf" ]
}

# Plot ranking of genes using a dot plot for markers

@test "Run Plot ranking of genes using a dot plot" {
    if [ "$resume" = 'true' ] && [ -f "$plt_rank_genes_groups_dot_pdf" ]; then
        skip "$plt_rank_genes_groups_dot_pdf exists and resume is set to 'true'"
    fi

    run rm -f $plt_rank_genes_groups_dot_pdf && eval "$scanpy plot dot $plt_rank_genes_groups_opt $diffexp_obj $plt_rank_genes_groups_dot_pdf"

    [ "$status" -eq 0 ]
    [ -f  "$plt_rank_genes_groups_dot_pdf" ]
}

# Plot ranking of genes using a dot plot for markers, high resolution clustering

@test "Run Plot ranking of genes using a dot plot, high resolution clustering" {
    if [ "$resume" = 'true' ] && [ -f "$plt_rank_genes_groups_dot_singlet_pdf" ]; then
        skip "$plt_rank_genes_groups_dot_singlet_pdf exists and resume is set to 'true'"
    fi

    run rm -f $plt_rank_genes_groups_dot_singlet_pdf && eval "$scanpy plot dot $plt_rank_genes_groups_singlet_opt $diffexp_singlet_obj $plt_rank_genes_groups_dot_singlet_pdf"

    [ "$status" -eq 0 ]
    [ -f  "$plt_rank_genes_groups_dot_singlet_pdf" ]
}

# Plot a matrix plot for markers

@test "Run Plot matrix" {
    if [ "$resume" = 'true' ] && [ -f "$plt_matrixplot_pdf" ]; then
        skip "$plt_matrixplot_pdf exists and resume is set to 'true'"
    fi

    run rm -f $plt_matrixplot_pdf && eval "$scanpy plot matrix $diffexp_plot_opt $diffexp_obj $plt_matrixplot_pdf"
    
    [ "$status" -eq 0 ]
    [ -f  "$plt_matrixplot_pdf" ]
}

# Plot ranking of genes using a matrix plot for markers

@test "Run Plot ranking of genes using a matrix plot" {
    if [ "$resume" = 'true' ] && [ -f "$plt_rank_genes_groups_matrix_pdf" ]; then
        skip "$plt_rank_genes_groups_matrix_pdf exists and resume is set to 'true'"
    fi

    run rm -f $plt_rank_genes_groups_matrix_pdf && eval "$scanpy plot matrix $plt_rank_genes_groups_opt $diffexp_obj $plt_rank_genes_groups_matrix_pdf"

    [ "$status" -eq 0 ]
    [ -f  "$plt_rank_genes_groups_matrix_pdf" ]
}

# Plot a matrix plot for markers

@test "Run Heatmap" {
    if [ "$resume" = 'true' ] && [ -f "$plt_heatmap_pdf" ]; then
        skip "$plt_matrixplot_pdf exists and resume is set to 'true'"
    fi

    run rm -f $plt_heatmap_pdf && eval "$scanpy plot heat $diffexp_plot_opt $diffexp_obj $plt_heatmap_pdf"
    
    [ "$status" -eq 0 ]
    [ -f  "$plt_heatmap_pdf" ]
}

# Plot ranking of genes using a matrix plot for markers

@test "Run Plot ranking of genes using a heatmap" {
    if [ "$resume" = 'true' ] && [ -f "$plt_rank_genes_groups_heatmap_pdf" ]; then
        skip "$plt_rank_genes_groups_heatmap_pdf exists and resume is set to 'true'"
    fi

    run rm -f $plt_rank_genes_groups_heatmap_pdf && eval "$scanpy plot heat $plt_rank_genes_groups_opt $diffexp_obj $plt_rank_genes_groups_heatmap_pdf"

    [ "$status" -eq 0 ]
    [ -f  "$plt_rank_genes_groups_matrix_pdf" ]
}

# Do harmony batch correction, using clustering as batch (just for test purposes)

@test "Run Harmony batch integration using clustering as batch" {
    if [ "$resume" = 'true' ] && [ -f "$harmony_integrate_obj" ]; then
        skip "$harmony_integrate_obj exists and resume is set to 'true'"
    fi

    run rm -f $harmony_integrate_obj && eval "$scanpy integrate harmony $harmony_integrate_opt $louvain_obj $harmony_integrate_obj"

    [ "$status" -eq 0 ]
    [ -f  "$plt_rank_genes_groups_matrix_pdf" ]

}

# Run Plot PCA embedding before harmony

@test "Run Plot PCA embedding before Harmony" {
    if [ "$resume" = 'true' ] && [ -f "$noharmony_integrated_pca_pdf" ]; then
        skip "$noharmony_integrated_pca_pdf exists and resume is set to 'true'"
    fi

    run rm -f $noharmony_integrated_pca_pdf && eval "$scanpy plot embed $noharmony_plt_embed_opt $louvain_obj $noharmony_integrated_pca_pdf"

    [ "$status" -eq 0 ]
    [ -f  "$noharmony_integrated_pca_pdf" ]
}

# Run Plot PCA embedding after harmony

@test "Run Plot PCA embedding after Harmony" {
    if [ "$resume" = 'true' ] && [ -f "$harmony_integrated_pca_pdf" ]; then
        skip "$harmony_integrated_pca_pdf exists and resume is set to 'true'"
    fi

    run rm -f $harmony_integrated_pca_pdf && eval "$scanpy plot embed $harmony_plt_embed_opt $harmony_integrate_obj $harmony_integrated_pca_pdf"

    [ "$status" -eq 0 ]
    [ -f  "$harmony_integrated_pca_pdf" ]
}

# Do bbknn batch correction, using clustering as batch (just for test purposes)

@test "Run BBKNN batch integration using clustering as batch" {
    if [ "$resume" = 'true' ] && [ -f "$bbknn_obj" ]; then
        skip "$bbknn_obj exists and resume is set to 'true'"
    fi

    run rm -f $bbknn_obj && eval "$scanpy integrate bbknn $bbknn_opt $louvain_obj $bbknn_obj"

    [ "$status" -eq 0 ]
    [ -f  "$plt_rank_genes_groups_matrix_pdf" ]
}

# Do MNN batch correction, using clustering as batch (just for test purposes)

@test "Run MNN batch integration using clustering as batch" {
    if [ "$resume" = 'true' ] && [ -f "$mnn_obj" ]; then
        skip "$mnn_obj exists and resume is set to 'true'"
    fi

    run rm -f $mnn_obj && eval "$scanpy integrate mnn $mnn_opt $louvain_obj $mnn_obj"

    [ "$status" -eq 0 ]
    [ -f  "$mnn_obj" ]
}

# Do ComBat batch correction, using clustering as batch (just for test purposes)

@test "Run Combat batch integration using clustering as batch" {
    if [ "$resume" = 'true' ] && [ -f "$combat_obj" ]; then
        skip "$combat_obj exists and resume is set to 'true'"
    fi

    run rm -f $combat_obj && eval "$scanpy integrate combat $combat_opt $louvain_obj $combat_obj"

    [ "$status" -eq 0 ]
    [ -f  "$combat_obj" ]
}

# Local Variables:
# mode: sh
# End:
