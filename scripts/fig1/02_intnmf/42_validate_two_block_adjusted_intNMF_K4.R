.libPaths(c(normalizePath("Rlibs", winslash="/", mustWork=TRUE), .libPaths()))

suppressPackageStartupMessages({
  library(IntNMF)
  library(NMF)
  library(cluster)
  library(mclust)
  library(readr)
  library(dplyr)
  library(tibble)
  library(pheatmap)
})

set.seed(42)

input_dir <- "results/fig1_inputs/generated"
outdir <- "results/intNMF_adjusted_two_block_K4_final"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(input_dir, "two_block_intNMF_inputs.rds")
qc_path <- file.path(input_dir, "two_block_sample_QC.csv")
meta_path <- "results/fig1_inputs/generated/intNMF_sample_metadata.csv"

stopifnot(
  file.exists(input_path),
  file.exists(qc_path),
  file.exists(meta_path)
)

dat <- readRDS(input_path)

qc <- read_csv(qc_path, show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

meta <- read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(
    sample_id = as.character(sample_id),
    subgroup = factor(subgroup)
  )

sample_order <- meta$sample_id

if (length(sample_order) != 84 || anyDuplicated(sample_order) > 0) {
  stop("Expected exactly 84 unique samples.")
}

for (nm in names(dat)) {
  if (!identical(rownames(dat[[nm]]), sample_order)) {
    stop(nm, ": sample order mismatch.")
  }
}

cat("===== Two-block input =====\n")
for (nm in names(dat)) {
  cat(
    nm, ":",
    nrow(dat[[nm]]), "samples x",
    ncol(dat[[nm]]), "features\n"
  )
}

seeds <- c(
  11, 23, 37, 41, 53,
  67, 79, 83, 97, 101,
  113, 127, 139, 149, 163,
  173, 181, 193, 211, 227
)

labels_matrix <- matrix(
  NA_integer_,
  nrow = length(sample_order),
  ncol = length(seeds),
  dimnames = list(
    sample_order,
    paste0("seed_", seeds)
  )
)

run_metrics <- vector("list", length(seeds))

for (i in seq_along(seeds)) {
  s <- seeds[i]

  cat(
    "K4 run",
    i,
    "of",
    length(seeds),
    "- seed",
    s,
    "\n"
  )

  set.seed(s)

  fit <- nmf.mnnals(
    dat = dat,
    k = 4,
    maxiter = 500,
    st.count = 50,
    n.ini = 20,
    ini.nndsvd = TRUE,
    seed = TRUE,
    wt = c(1, 1)
  )

  clusters <- as.integer(fit$clusters)
  labels_matrix[, i] <- clusters

  consensus <- as.matrix(fit$consensus)
  rownames(consensus) <- sample_order
  colnames(consensus) <- sample_order

  sil <- silhouette(
    clusters,
    as.dist(1 - consensus)
  )

  tab <- table(clusters)

  run_metrics[[i]] <- tibble(
    seed = s,
    mean_silhouette = mean(sil[, "sil_width"]),
    cophenetic_correlation = tryCatch(
      NMF::cophcor(consensus),
      error = function(e) NA_real_
    ),
    dispersion = tryCatch(
      NMF::dispersion(consensus),
      error = function(e) NA_real_
    ),
    minimum_cluster_size = min(tab),
    maximum_cluster_size = max(tab),
    cluster_sizes = paste(
      paste0("C", names(tab)),
      as.integer(tab),
      sep = "=",
      collapse = ";"
    )
  )

  rm(fit, consensus)
  invisible(gc())
}

run_metrics_df <- bind_rows(run_metrics)

write_csv(
  run_metrics_df,
  file.path(outdir, "two_block_K4_seed_run_metrics.csv")
)

write_csv(
  as.data.frame(labels_matrix) %>%
    rownames_to_column("sample_id"),
  file.path(outdir, "two_block_K4_labels_across_seeds.csv")
)

# Pairwise ARI
n_runs <- ncol(labels_matrix)

ari_matrix <- matrix(
  1,
  nrow = n_runs,
  ncol = n_runs,
  dimnames = list(
    colnames(labels_matrix),
    colnames(labels_matrix)
  )
)

for (i in seq_len(n_runs)) {
  for (j in seq_len(n_runs)) {
    if (j <= i) next

    ari_value <- adjustedRandIndex(
      labels_matrix[, i],
      labels_matrix[, j]
    )

    ari_matrix[i, j] <- ari_value
    ari_matrix[j, i] <- ari_value
  }
}

pairwise_ari <- ari_matrix[upper.tri(ari_matrix)]

mean_ari_by_seed <- (
  rowSums(ari_matrix) - 1
) / (n_runs - 1)

medoid_index <- which.max(mean_ari_by_seed)
medoid_seed <- seeds[medoid_index]

stability <- tibble(
  number_of_seed_runs = n_runs,
  pairwise_ARI_mean = mean(pairwise_ari),
  pairwise_ARI_median = median(pairwise_ari),
  pairwise_ARI_q25 = as.numeric(quantile(pairwise_ari, 0.25)),
  pairwise_ARI_q75 = as.numeric(quantile(pairwise_ari, 0.75)),
  pairwise_ARI_min = min(pairwise_ari),
  pairwise_ARI_max = max(pairwise_ari),
  medoid_seed = medoid_seed
)

write_csv(
  stability,
  file.path(outdir, "two_block_K4_stability_summary.csv")
)

write_csv(
  as.data.frame(ari_matrix) %>%
    rownames_to_column("run"),
  file.path(outdir, "two_block_K4_pairwise_ARI_matrix.csv")
)

cat("\n===== K4 seed stability =====\n")
print(stability)

# Refit using medoid seed
cat("\n===== Final K4 fit; medoid seed =", medoid_seed, "=====\n")

set.seed(medoid_seed)

final_fit <- nmf.mnnals(
  dat = dat,
  k = 4,
  maxiter = 500,
  st.count = 50,
  n.ini = 20,
  ini.nndsvd = TRUE,
  seed = TRUE,
  wt = c(1, 1)
)

final_clusters <- as.integer(final_fit$clusters)

final_labels <- tibble(
  sample_id = sample_order,
  adjusted_two_block_K4 = paste0("C", final_clusters)
) %>%
  left_join(meta, by = "sample_id") %>%
  arrange(
    adjusted_two_block_K4,
    subgroup,
    sample_id
  )

write_csv(
  final_labels,
  file.path(outdir, "two_block_final_K4_labels.csv")
)

final_consensus <- as.matrix(final_fit$consensus)
rownames(final_consensus) <- sample_order
colnames(final_consensus) <- sample_order

write_csv(
  as.data.frame(final_consensus) %>%
    rownames_to_column("sample_id"),
  file.path(outdir, "two_block_final_K4_consensus_matrix.csv")
)

saveRDS(
  final_fit,
  file.path(outdir, "two_block_final_K4_fit.rds")
)

final_sil <- silhouette(
  final_clusters,
  as.dist(1 - final_consensus)
)

cluster_table <- table(final_labels$adjusted_two_block_K4)

final_metrics <- tibble(
  K = 4,
  medoid_seed = medoid_seed,
  mean_silhouette = mean(final_sil[, "sil_width"]),
  cophenetic_correlation = tryCatch(
    NMF::cophcor(final_consensus),
    error = function(e) NA_real_
  ),
  dispersion = tryCatch(
    NMF::dispersion(final_consensus),
    error = function(e) NA_real_
  ),
  minimum_cluster_size = min(cluster_table),
  maximum_cluster_size = max(cluster_table),
  cluster_sizes = paste(
    names(cluster_table),
    as.integer(cluster_table),
    sep = "=",
    collapse = ";"
  )
)

write_csv(
  final_metrics,
  file.path(outdir, "two_block_final_K4_metrics.csv")
)

# QC audit
qc_final <- qc %>%
  left_join(
    final_labels %>%
      select(sample_id, adjusted_two_block_K4),
    by = "sample_id"
  )

qc_metrics <- c(
  "log10_RNA_library_size",
  "RNA_detected_features",
  "RNA_zero_fraction",
  "log10_smallRNA_library_size",
  "smallRNA_detected_features",
  "smallRNA_zero_fraction"
)

qc_tests <- lapply(
  qc_metrics,
  function(metric) {
    form <- as.formula(
      paste(metric, "~ adjusted_two_block_K4")
    )

    tibble(
      metric = metric,
      pvalue = kruskal.test(
        form,
        data = qc_final
      )$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(pvalue, method = "BH")
  ) %>%
  arrange(FDR)

write_csv(
  qc_tests,
  file.path(outdir, "two_block_final_K4_QC_tests.csv")
)

# Cancer-type composition
cancer_table <- table(
  final_labels$adjusted_two_block_K4,
  final_labels$subgroup
)

cancer_test <- fisher.test(cancer_table)

write_csv(
  tibble(
    test = "K4 versus cancer type",
    pvalue = cancer_test$p.value
  ),
  file.path(outdir, "two_block_final_K4_cancer_type_test.csv")
)

# Consensus heatmap
sample_plot_order <- final_labels$sample_id

annotation <- final_labels %>%
  select(
    sample_id,
    adjusted_two_block_K4,
    subgroup
  ) %>%
  as.data.frame()

rownames(annotation) <- annotation$sample_id
annotation$sample_id <- NULL

pdf(
  file.path(outdir, "two_block_final_K4_consensus_heatmap.pdf"),
  width = 7,
  height = 7
)

pheatmap(
  final_consensus[
    sample_plot_order,
    sample_plot_order
  ],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = annotation[
    sample_plot_order,
    ,
    drop = FALSE
  ],
  annotation_col = annotation[
    sample_plot_order,
    ,
    drop = FALSE
  ],
  border_color = NA,
  main = "Two-block technical-adjusted intNMF, K = 4"
)

dev.off()

cat("\n===== FINAL K4 RESULTS =====\n")

cat("\nCluster sizes:\n")
print(cluster_table)

cat("\nMetrics:\n")
print(final_metrics)

cat("\nQC tests:\n")
print(qc_tests)

cat("\nK4 by cancer type:\n")
print(cancer_table)

cat("\nCancer-type Fisher P:", cancer_test$p.value, "\n")

cat("\nOutput directory:", outdir, "\n")
print(list.files(outdir, full.names = TRUE))
