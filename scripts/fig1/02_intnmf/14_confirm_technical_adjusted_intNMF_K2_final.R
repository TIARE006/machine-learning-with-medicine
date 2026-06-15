.libPaths(c(normalizePath("Rlibs", winslash="/", mustWork=TRUE), .libPaths()))

local_lib <- normalizePath("Rlibs", winslash="/", mustWork=TRUE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

need_pkgs <- c(
  "IntNMF", "NMF", "cluster", "mclust",
  "readr", "dplyr", "tidyr", "tibble",
  "stringr", "ggplot2", "pheatmap"
)

for (p in need_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      lib = local_lib,
      dependencies = c("Depends", "Imports", "LinkingTo")
    )
  }
}

suppressPackageStartupMessages({
  library(IntNMF)
  library(NMF)
  library(cluster)
  library(mclust)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
})

set.seed(42)

pilot_dir <- "results/intNMF_strict84_technical_adjusted_pilot"
outdir <- "results/intNMF_strict84_technical_adjusted_final"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(
  pilot_dir,
  "technical_adjusted_intNMF_input_blocks.rds"
)

meta_path <- "results/fig1_inputs/generated/intNMF_sample_metadata.csv"

qc_path <- file.path(
  pilot_dir,
  "technical_adjustment_QC_table.csv"
)

stopifnot(file.exists(input_path))
stopifnot(file.exists(meta_path))
stopifnot(file.exists(qc_path))

dat <- readRDS(input_path)

meta <- read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(
    sample_id = as.character(sample_id),
    subgroup = factor(subgroup)
  )

qc <- read_csv(qc_path, show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

sample_order <- meta$sample_id

if (
  length(sample_order) != 84 ||
  anyDuplicated(sample_order) > 0
) {
  stop("Expected exactly 84 unique samples.")
}

for (nm in names(dat)) {
  if (!identical(rownames(dat[[nm]]), sample_order)) {
    stop(nm, ": sample order mismatch.")
  }

  if (
    any(!is.finite(dat[[nm]])) ||
    any(dat[[nm]] < 0)
  ) {
    stop(nm, ": invalid values in adjusted input block.")
  }
}

cat("===== Technical-adjusted intNMF inputs =====\n")

for (nm in names(dat)) {
  cat(
    nm, ":",
    nrow(dat[[nm]]), "samples x",
    ncol(dat[[nm]]), "features\n"
  )
}

# ============================================================
# Part 1: Full CPI confirmation
# ============================================================

cat("\n===== Full CPI confirmation =====\n")
cat("K range: 2 to 6\n")
cat("n.runs: 30\n")
cat("n.fold: 5\n")

k_range <- 2:6

set.seed(20260607)

pdf(
  file.path(outdir, "technical_adjusted_final_CPI_native_plot.pdf"),
  width = 6.5,
  height = 5
)

cpi <- nmf.opt.k(
  dat = dat,
  n.runs = 30,
  n.fold = 5,
  k.range = k_range,
  result = TRUE,
  make.plot = TRUE,
  progress = TRUE,
  st.count = 30,
  maxiter = 300,
  wt = rep(1, length(dat))
)

dev.off()

cpi <- as.matrix(cpi)

if (
  nrow(cpi) != length(k_range) &&
  ncol(cpi) == length(k_range)
) {
  cpi <- t(cpi)
}

if (nrow(cpi) != length(k_range)) {
  stop(
    "Unexpected CPI dimensions: ",
    paste(dim(cpi), collapse = " x ")
  )
}

rownames(cpi) <- paste0("K", k_range)

write_csv(
  as.data.frame(cpi) %>%
    rownames_to_column("K"),
  file.path(outdir, "technical_adjusted_final_CPI_all_runs.csv")
)

cpi_summary <- tibble(
  K = k_range,
  CPI_mean = rowMeans(cpi, na.rm = TRUE),
  CPI_median = apply(cpi, 1, median, na.rm = TRUE),
  CPI_sd = apply(cpi, 1, sd, na.rm = TRUE),
  CPI_q25 = apply(cpi, 1, quantile, probs = 0.25, na.rm = TRUE),
  CPI_q75 = apply(cpi, 1, quantile, probs = 0.75, na.rm = TRUE),
  CPI_min = apply(cpi, 1, min, na.rm = TRUE),
  CPI_max = apply(cpi, 1, max, na.rm = TRUE)
)

best_k <- cpi_summary$K[
  which.max(cpi_summary$CPI_mean)
]

cpi_summary <- cpi_summary %>%
  mutate(final_best_K = K == best_k)

write_csv(
  cpi_summary,
  file.path(outdir, "technical_adjusted_final_CPI_summary.csv")
)

cat("\n===== Final CPI summary =====\n")
print(cpi_summary)
cat("\nBest K by mean CPI:", best_k, "\n")

cpi_long <- as.data.frame(cpi) %>%
  rownames_to_column("K") %>%
  mutate(K = as.integer(str_remove(K, "^K"))) %>%
  pivot_longer(
    cols = -K,
    names_to = "run",
    values_to = "CPI"
  )

p_cpi <- ggplot(
  cpi_long,
  aes(x = factor(K), y = CPI)
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.4,
    alpha = 0.55
  ) +
  geom_line(
    data = cpi_summary,
    aes(
      x = factor(K),
      y = CPI_mean,
      group = 1
    ),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_point(
    data = cpi_summary,
    aes(
      x = factor(K),
      y = CPI_mean
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Technical-adjusted intNMF: CPI-based K selection",
    subtitle = paste0("30 runs; best K = ", best_k),
    x = "Number of clusters",
    y = "Cluster Prediction Index"
  )

ggsave(
  file.path(outdir, "technical_adjusted_final_CPI_by_K.pdf"),
  p_cpi,
  width = 6.5,
  height = 5
)

ggsave(
  file.path(outdir, "technical_adjusted_final_CPI_by_K.png"),
  p_cpi,
  width = 6.5,
  height = 5,
  dpi = 300
)

# ============================================================
# Part 2: Repeated K=2 fits across random seeds
# ============================================================

cat("\n===== Repeated K=2 stability analysis =====\n")

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

run_metrics <- vector(
  "list",
  length(seeds)
)

for (i in seq_along(seeds)) {
  s <- seeds[i]

  cat(
    "K2 stability run",
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
    k = 2,
    maxiter = 500,
    st.count = 50,
    n.ini = 20,
    ini.nndsvd = TRUE,
    seed = TRUE,
    wt = rep(1, length(dat))
  )

  cluster_i <- as.integer(fit$clusters)

  if (length(cluster_i) != length(sample_order)) {
    stop("Unexpected label count for seed ", s)
  }

  labels_matrix[, i] <- cluster_i

  consensus_i <- as.matrix(fit$consensus)
  rownames(consensus_i) <- sample_order
  colnames(consensus_i) <- sample_order

  sil_i <- silhouette(
    cluster_i,
    as.dist(1 - consensus_i)
  )

  tab_i <- table(cluster_i)

  run_metrics[[i]] <- tibble(
    seed = s,
    mean_silhouette = mean(sil_i[, "sil_width"]),
    cophenetic_correlation = tryCatch(
      NMF::cophcor(consensus_i),
      error = function(e) NA_real_
    ),
    dispersion = tryCatch(
      NMF::dispersion(consensus_i),
      error = function(e) NA_real_
    ),
    minimum_cluster_size = min(tab_i),
    maximum_cluster_size = max(tab_i),
    cluster_sizes = paste(
      paste0("C", names(tab_i)),
      as.integer(tab_i),
      sep = "=",
      collapse = ";"
    )
  )

  rm(fit, consensus_i)
  invisible(gc())
}

run_metrics_df <- bind_rows(run_metrics)

write_csv(
  run_metrics_df,
  file.path(outdir, "technical_adjusted_K2_seed_run_metrics.csv")
)

write_csv(
  as.data.frame(labels_matrix) %>%
    rownames_to_column("sample_id"),
  file.path(outdir, "technical_adjusted_K2_labels_across_seeds.csv")
)

# Pairwise adjusted Rand index
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

mean_ari_by_run <- rowMeans(
  ari_matrix[
    ,
    ,
    drop = FALSE
  ]
)

diag_adjustment <- 1 / n_runs

mean_ari_excluding_self <- (
  rowSums(ari_matrix) - 1
) / (n_runs - 1)

stability_summary <- tibble(
  seed = seeds,
  mean_ARI_to_other_runs = as.numeric(mean_ari_excluding_self)
) %>%
  left_join(
    run_metrics_df,
    by = "seed"
  ) %>%
  arrange(desc(mean_ARI_to_other_runs))

medoid_seed <- stability_summary$seed[1]

pairwise_ari <- ari_matrix[
  upper.tri(ari_matrix)
]

overall_stability <- tibble(
  number_of_seed_runs = n_runs,
  pairwise_ARI_mean = mean(pairwise_ari),
  pairwise_ARI_median = median(pairwise_ari),
  pairwise_ARI_q25 = quantile(pairwise_ari, 0.25),
  pairwise_ARI_q75 = quantile(pairwise_ari, 0.75),
  pairwise_ARI_min = min(pairwise_ari),
  pairwise_ARI_max = max(pairwise_ari),
  medoid_seed = medoid_seed
)

write_csv(
  as.data.frame(ari_matrix) %>%
    rownames_to_column("run"),
  file.path(outdir, "technical_adjusted_K2_pairwise_ARI_matrix.csv")
)

write_csv(
  stability_summary,
  file.path(outdir, "technical_adjusted_K2_seed_stability.csv")
)

write_csv(
  overall_stability,
  file.path(outdir, "technical_adjusted_K2_overall_stability.csv")
)

cat("\n===== K2 seed stability =====\n")
print(overall_stability)
cat("\nSelected medoid seed:", medoid_seed, "\n")

# ============================================================
# Part 3: Final K=2 fit using the medoid seed
# ============================================================

cat("\n===== Final medoid-seed K=2 fit =====\n")

set.seed(medoid_seed)

final_fit <- nmf.mnnals(
  dat = dat,
  k = 2,
  maxiter = 500,
  st.count = 50,
  n.ini = 20,
  ini.nndsvd = TRUE,
  seed = TRUE,
  wt = rep(1, length(dat))
)

final_clusters <- as.integer(
  final_fit$clusters
)

final_labels <- tibble(
  sample_id = sample_order,
  adjusted_intNMF_K2 = paste0("C", final_clusters)
) %>%
  left_join(
    meta,
    by = "sample_id"
  ) %>%
  arrange(
    adjusted_intNMF_K2,
    subgroup,
    sample_id
  )

write_csv(
  final_labels,
  file.path(outdir, "technical_adjusted_final_K2_labels.csv")
)

final_consensus <- as.matrix(
  final_fit$consensus
)

rownames(final_consensus) <- sample_order
colnames(final_consensus) <- sample_order

write_csv(
  as.data.frame(final_consensus) %>%
    rownames_to_column("sample_id"),
  file.path(outdir, "technical_adjusted_final_K2_consensus_matrix.csv")
)

final_silhouette <- silhouette(
  final_clusters,
  as.dist(1 - final_consensus)
)

final_tab <- table(
  final_labels$adjusted_intNMF_K2
)

final_metrics <- tibble(
  K = 2,
  medoid_seed = medoid_seed,
  mean_silhouette = mean(
    final_silhouette[, "sil_width"]
  ),
  cophenetic_correlation = tryCatch(
    NMF::cophcor(final_consensus),
    error = function(e) NA_real_
  ),
  dispersion = tryCatch(
    NMF::dispersion(final_consensus),
    error = function(e) NA_real_
  ),
  minimum_cluster_size = min(final_tab),
  maximum_cluster_size = max(final_tab),
  cluster_sizes = paste(
    names(final_tab),
    as.integer(final_tab),
    sep = "=",
    collapse = ";"
  )
)

write_csv(
  final_metrics,
  file.path(outdir, "technical_adjusted_final_K2_metrics.csv")
)

saveRDS(
  final_fit,
  file.path(outdir, "technical_adjusted_final_K2_fit.rds")
)

# ============================================================
# Part 4: Final QC association audit
# ============================================================

qc_final <- qc %>%
  left_join(
    final_labels %>%
      select(
        sample_id,
        adjusted_intNMF_K2
      ),
    by = "sample_id"
  )

qc_metrics <- colnames(qc_final)[
  grepl(
    "log10_|detected_features|zero_fraction",
    colnames(qc_final)
  )
]

qc_tests <- lapply(
  qc_metrics,
  function(metric) {
    x1 <- qc_final[[metric]][
      qc_final$adjusted_intNMF_K2 == "C1"
    ]

    x2 <- qc_final[[metric]][
      qc_final$adjusted_intNMF_K2 == "C2"
    ]

    wt <- wilcox.test(
      x2,
      x1,
      exact = FALSE
    )

    tibble(
      metric = metric,
      C1_median = median(x1, na.rm = TRUE),
      C2_median = median(x2, na.rm = TRUE),
      difference_C2_minus_C1 =
        median(x2, na.rm = TRUE) -
        median(x1, na.rm = TRUE),
      pvalue = wt$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(
      pvalue,
      method = "BH"
    )
  ) %>%
  arrange(FDR)

cancer_test <- fisher.test(
  table(
    qc_final$adjusted_intNMF_K2,
    qc_final$subgroup
  )
)

write_csv(
  qc_tests,
  file.path(outdir, "technical_adjusted_final_K2_QC_tests.csv")
)

write_csv(
  tibble(
    test = "K2 versus cancer type",
    pvalue = cancer_test$p.value
  ),
  file.path(outdir, "technical_adjusted_final_K2_cancer_type_test.csv")
)

cat("\n===== Final K2 cluster table =====\n")
print(final_tab)

cat("\n===== Final K2 metrics =====\n")
print(final_metrics)

cat("\n===== Final QC tests =====\n")
print(qc_tests)

cat("\n===== Final K2 by cancer type =====\n")
print(
  table(
    final_labels$adjusted_intNMF_K2,
    final_labels$subgroup
  )
)

# ============================================================
# Part 5: Consensus heatmap
# ============================================================

sample_plot_order <- final_labels %>%
  arrange(
    adjusted_intNMF_K2,
    subgroup,
    sample_id
  ) %>%
  pull(sample_id)

annotation <- final_labels %>%
  select(
    sample_id,
    adjusted_intNMF_K2,
    subgroup
  ) %>%
  as.data.frame()

rownames(annotation) <- annotation$sample_id
annotation$sample_id <- NULL

pdf(
  file.path(outdir, "technical_adjusted_final_K2_consensus_heatmap.pdf"),
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
  main = "Technical-adjusted intNMF consensus, K=2"
)

dev.off()

# ============================================================
# Part 6: Publication decision
# ============================================================

passes_best_k <- best_k == 2

passes_seed_stability <-
  overall_stability$pairwise_ARI_median >= 0.80

passes_cluster_size <-
  final_metrics$minimum_cluster_size >= 10

passes_qc <-
  all(
    qc_tests$FDR >= 0.05 |
      is.na(qc_tests$FDR)
  )

publication_status <- if (
  passes_best_k &&
  passes_seed_stability &&
  passes_cluster_size &&
  passes_qc
) {
  "PASS: technical-adjusted K=2 is supported for the main analysis."
} else {
  "FAIL OR INCONCLUSIVE: do not present K=2 as the final biological subtype solution."
}

decision <- tibble(
  criterion = c(
    "K2 has highest mean CPI",
    "Median pairwise seed ARI >= 0.80",
    "Minimum cluster size >= 10",
    "No measured QC association at FDR < 0.05"
  ),
  passed = c(
    passes_best_k,
    passes_seed_stability,
    passes_cluster_size,
    passes_qc
  )
)

write_csv(
  decision,
  file.path(outdir, "technical_adjusted_final_decision_criteria.csv")
)

writeLines(
  c(
    publication_status,
    paste0("Best K by 30-run CPI: ", best_k),
    paste0(
      "Median pairwise K2 ARI: ",
      overall_stability$pairwise_ARI_median
    ),
    paste0(
      "Final K2 cluster sizes: ",
      final_metrics$cluster_sizes
    ),
    paste0(
      "Smallest QC FDR: ",
      min(qc_tests$FDR, na.rm = TRUE)
    )
  ),
  file.path(outdir, "technical_adjusted_final_decision.txt")
)

cat("\n============================================================\n")
cat(publication_status, "\n")
cat("Best K:", best_k, "\n")
cat(
  "Median pairwise K2 ARI:",
  overall_stability$pairwise_ARI_median,
  "\n"
)
cat("Output directory:", outdir, "\n")
cat("============================================================\n")
