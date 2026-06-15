.libPaths(c(normalizePath("Rlibs", winslash = "/", mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(digest)
})

lineage <- tribble(
  ~panel, ~direct_input, ~producer, ~root_sources,
  "A,D,E", "results/fig1_inputs/generated/GSE254877_14PCD_ssGSEA_zscores.csv", "scripts/fig1/01_expression_inputs/10_build_fig1_expression_inputs.R", "GSE254877 RNA counts + pcd_gene_sets.csv",
  "B-top", "results/intNMF_reproducible_two_block/CPI_summary.csv", "scripts/fig1/02_intnmf/20_train_reproducible_two_block_intnmf.R", "raw RNA/lncRNA/small-RNA + locked feature specification",
  "B-bottom,C", "results/intNMF_reproducible_two_block/final_consensus_matrix.csv", "scripts/fig1/02_intnmf/20_train_reproducible_two_block_intnmf.R", "raw RNA/lncRNA/small-RNA + locked feature specification",
  "D,E,F,G,H", "results/intNMF_reproducible_two_block/two_block_final_K4_labels.csv", "scripts/fig1/04_annotate_reproducible_K4.R", "reproducible K4 fit + locked pyroptosis annotation rule",
  "F,G", "results/fig1_final/all84_pyro_module_Euclidean_K2_main/all84_pyro_module_scores_by_sample.csv", "scripts/fig1/03_pyroptosis/32_all84_pyro_module_main_figure.R", "GSE254877 RNA counts + final K4 labels + locked module genes",
  "F,G", "results/fig1_final/all84_pyro_module_Euclidean_K2_main/all84_K2_validation_scores_by_sample.csv", "scripts/fig1/03_pyroptosis/32_all84_pyro_module_main_figure.R", "GSE254877 RNA counts + final K4 labels + held-out genes",
  "G,H", "results/fig1_final/all84_pyro_module_Euclidean_K2_main/all84_sample_metadata_used.csv", "scripts/fig1/03_pyroptosis/32_all84_pyro_module_main_figure.R", "GEO phenotype + expression-inferred sex + final K4 labels",
  "H", "results/fig1_inputs/generated/GSE254877_pyroptosis_marker_gene_zmatrix.csv", "scripts/fig1/01_expression_inputs/10_build_fig1_expression_inputs.R", "GSE254877 RNA counts + fixed marker list"
)

missing_inputs <- lineage$direct_input[!file.exists(lineage$direct_input)]
producer_paths <- unique(unlist(strsplit(lineage$producer, "; ", fixed = TRUE)))
missing_producers <- producer_paths[!file.exists(producer_paths)]

root_files <- c(
  "data/raw/RNA_seq/GSE254877_raw_counts_expression.csv",
  "data/raw/lncRNA_seq/GSE254877_lncRNA_raw_counts_expression.csv",
  "data/raw/small_RNA_seq/GSE254878_smallRNAs_raw_counts_expression.csv",
  "data/raw/RNA_seq/GSE254877_pheno.csv",
  "scripts/fig1/01_expression_inputs/pcd_gene_sets.csv",
  "scripts/fig1/01_expression_inputs/locked_intnmf_feature_spec.csv",
  "results/intNMF_reproducible_two_block/final_fit.rds",
  "results/intNMF_reproducible_two_block/CPI_all_runs.csv"
)
missing_roots <- root_files[!file.exists(root_files)]

manifest_paths <- c(
  "results/fig1_inputs/generated/expression_input_lineage_manifest.csv",
  "results/fig1_inputs/generated/intNMF_input_lineage_manifest.csv",
  "results/intNMF_reproducible_two_block/training_lineage_manifest.csv",
  "results/intNMF_reproducible_two_block/annotation_lineage_manifest.csv"
)
missing_manifests <- manifest_paths[!file.exists(manifest_paths)]

verify_manifest <- function(path) {
  manifest <- read_csv(path, show_col_types = FALSE)
  manifest %>%
    rowwise() %>%
    mutate(
      exists = file.exists(path),
      current_sha256 = if (exists) digest(path, algo = "sha256", file = TRUE) else NA_character_,
      hash_matches = exists && current_sha256 == sha256
    ) %>%
    ungroup()
}

hash_audit <- bind_rows(lapply(manifest_paths[file.exists(manifest_paths)], verify_manifest))

if (length(missing_inputs) || length(missing_producers) || length(missing_roots) ||
    length(missing_manifests) || any(!hash_audit$hash_matches)) {
  stop(
    "Fig. 1 lineage audit failed.\n",
    "Missing direct inputs: ", paste(missing_inputs, collapse = ", "), "\n",
    "Missing producers: ", paste(missing_producers, collapse = ", "), "\n",
    "Missing roots: ", paste(missing_roots, collapse = ", "), "\n",
    "Missing manifests: ", paste(missing_manifests, collapse = ", "), "\n",
    "Hash failures: ", paste(hash_audit$path[!hash_audit$hash_matches], collapse = ", ")
  )
}

dir.create("results/fig1_inputs", recursive = TRUE, showWarnings = FALSE)
write_csv(lineage, "results/fig1_inputs/Fig1_panel_data_lineage.csv")
write_csv(hash_audit, "results/fig1_inputs/Fig1_source_hash_audit.csv")
cat("Fig. 1 lineage audit passed for panels A-H.\n")
