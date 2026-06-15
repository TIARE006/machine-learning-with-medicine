.libPaths(c(normalizePath("Rlibs", winslash="/", mustWork=TRUE), .libPaths()))

local_lib <- normalizePath("Rlibs", winslash="/", mustWork=TRUE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_need <- c("IntNMF", "NMF", "cluster", "readr", "dplyr", "tidyr", "tibble",
               "stringr", "matrixStats", "ggplot2", "pheatmap", "snow")
for (p in cran_need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, lib = local_lib, dependencies = TRUE)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = local_lib)
}

bio_need <- c("edgeR", "limma")
for (p in bio_need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, lib = local_lib, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(IntNMF)
  library(NMF)
  library(cluster)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(edgeR)
  library(limma)
  library(matrixStats)
  library(ggplot2)
  library(pheatmap)
})

set.seed(42)

outdir <- "results/intNMF_strict84_technical_adjusted_pilot"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta_path <- "results/fig1_inputs/generated/intNMF_sample_metadata.csv"
mrna_path <- "data/raw/RNA_seq/GSE254877_raw_counts_expression.csv"
lncrna_path <- "data/raw/lncRNA_seq/GSE254877_lncRNA_raw_counts_expression.csv"
smallrna_path <- "data/raw/small_RNA_seq/GSE254878_smallRNAs_raw_counts_expression.csv"

stopifnot(file.exists(meta_path))
stopifnot(file.exists(mrna_path))
stopifnot(file.exists(lncrna_path))
stopifnot(file.exists(smallrna_path))

meta <- read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(
    sample_id = as.character(sample_id),
    subgroup = relevel(factor(subgroup), ref = "Colorectal")
  )

sample_order <- meta$sample_id

if (!identical(sample_order, paste0("AB_", 1:84))) {
  stop("Sample order is not AB_1 to AB_84.")
}

read_count_block <- function(path, label, allowed_types = NULL, uppercase_features = FALSE) {
  raw <- read_csv(
    path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    name_repair = "unique"
  )
  
  first_row <- as.character(unlist(raw[1, , drop = FALSE], use.names = FALSE))
  sample_cols <- which(str_detect(first_row, "^AB_[0-9]+$"))
  
  if (length(sample_cols) != 84) {
    stop(label, ": expected 84 AB sample columns, detected ", length(sample_cols))
  }
  if (!all(diff(sample_cols) == 1)) {
    stop(label, ": sample columns are not contiguous.")
  }
  
  dat <- raw[-1, , drop = FALSE]
  feature_col <- colnames(dat)[1]
  type_col <- if (min(sample_cols) >= 3) colnames(dat)[2] else NULL
  
  features <- as.character(dat[[feature_col]])
  if (uppercase_features) features <- toupper(features)
  
  keep <- !is.na(features) & features != ""
  
  if (!is.null(allowed_types)) {
    if (is.null(type_col)) stop(label, ": no type column detected.")
    feature_type <- as.character(dat[[type_col]])
    keep <- keep & feature_type %in% allowed_types
  }
  
  count_char <- as.matrix(dat[, sample_cols, drop = FALSE])
  count_mat <- apply(count_char, 2, function(x) suppressWarnings(as.numeric(x)))
  if (is.null(dim(count_mat))) count_mat <- matrix(count_mat, ncol = 1)
  count_mat[is.na(count_mat)] <- 0
  
  rownames(count_mat) <- features
  colnames(count_mat) <- first_row[sample_cols]
  
  count_mat <- count_mat[keep, , drop = FALSE]
  
  missing_samples <- setdiff(sample_order, colnames(count_mat))
  if (length(missing_samples) > 0) {
    stop(label, " missing samples: ", paste(missing_samples, collapse = ", "))
  }
  
  count_mat <- count_mat[, sample_order, drop = FALSE]
  count_mat <- rowsum(count_mat, group = rownames(count_mat), reorder = FALSE)
  count_mat <- as.matrix(count_mat)
  storage.mode(count_mat) <- "numeric"
  
  cat(label, ":", nrow(count_mat), "features x", ncol(count_mat), "samples\n")
  count_mat
}

make_qc <- function(count_mat, prefix) {
  tibble(
    sample_id = colnames(count_mat),
    !!paste0(prefix, "_library_size") := colSums(count_mat),
    !!paste0("log10_", prefix, "_library_size") := log10(colSums(count_mat) + 1),
    !!paste0(prefix, "_detected_features") := colSums(count_mat > 0),
    !!paste0(prefix, "_zero_fraction") := colMeans(count_mat == 0)
  )
}

cat("===== Read count blocks =====\n")

mrna_counts <- read_count_block(
  mrna_path,
  label = "mRNA",
  allowed_types = "protein_coding",
  uppercase_features = TRUE
)

lncrna_counts <- read_count_block(
  lncrna_path,
  label = "lncRNA",
  allowed_types = NULL,
  uppercase_features = TRUE
)

smallrna_counts <- read_count_block(
  smallrna_path,
  label = "smallRNA",
  allowed_types = NULL,
  uppercase_features = FALSE
)

qc <- make_qc(mrna_counts, "mRNA") %>%
  left_join(make_qc(lncrna_counts, "lncRNA"), by = "sample_id") %>%
  left_join(make_qc(smallrna_counts, "smallRNA"), by = "sample_id") %>%
  left_join(meta, by = "sample_id")

write_csv(qc, file.path(outdir, "technical_adjustment_QC_table.csv"))

prepare_adjusted_block <- function(counts, label, prefix, min_count, min_samples, top_n) {
  counts <- counts[, sample_order, drop = FALSE]
  
  keep <- rowSums(counts >= min_count) >= min_samples
  counts <- counts[keep, , drop = FALSE]
  
  dge <- DGEList(counts = round(counts))
  dge <- normLibSizes(dge)
  logcpm <- cpm(dge, log = TRUE, prior.count = 1)
  
  covar <- qc %>%
    arrange(match(sample_id, sample_order)) %>%
    transmute(
      detected = as.numeric(scale(.data[[paste0(prefix, "_detected_features")]])),
      loglib = as.numeric(scale(.data[[paste0("log10_", prefix, "_library_size")]]))
    ) %>%
    as.matrix()
  
  design_preserve <- model.matrix(~ subgroup, data = qc %>% arrange(match(sample_id, sample_order)))
  
  adjusted <- removeBatchEffect(
    logcpm,
    covariates = covar,
    design = design_preserve
  )
  
  vars <- rowVars(adjusted)
  valid <- is.finite(vars) & vars > 0
  adjusted <- adjusted[valid, , drop = FALSE]
  vars <- vars[valid]
  
  ord <- order(vars, decreasing = TRUE)
  n_use <- min(top_n, length(ord))
  selected <- adjusted[ord[seq_len(n_use)], , drop = FALSE]
  
  x <- t(selected)
  
  min_x <- min(x, na.rm = TRUE)
  if (min_x < 0) x <- x + abs(min_x)
  x <- pmax(x, .Machine$double.eps)
  x <- x / max(x, na.rm = TRUE)
  
  rownames(x) <- sample_order
  
  if (any(!is.finite(x)) || any(x < 0) || !identical(rownames(x), sample_order)) {
    stop(label, ": invalid adjusted intNMF input.")
  }
  
  cat(label, "adjusted block:", dim(x)[1], "samples x", dim(x)[2], "features\n")
  
  list(
    matrix = x,
    raw_features = nrow(counts),
    selected_features = colnames(x)
  )
}

cat("\n===== Build technical-adjusted intNMF blocks =====\n")

mrna <- prepare_adjusted_block(
  mrna_counts,
  label = "mRNA",
  prefix = "mRNA",
  min_count = 10,
  min_samples = 5,
  top_n = 3000
)

lncrna <- prepare_adjusted_block(
  lncrna_counts,
  label = "lncRNA",
  prefix = "lncRNA",
  min_count = 5,
  min_samples = 5,
  top_n = 1000
)

smallrna <- prepare_adjusted_block(
  smallrna_counts,
  label = "smallRNA",
  prefix = "smallRNA",
  min_count = 5,
  min_samples = 5,
  top_n = 500
)

dat <- list(
  mRNA = mrna$matrix,
  lncRNA = lncrna$matrix,
  smallRNA = smallrna$matrix
)

saveRDS(dat, file.path(outdir, "technical_adjusted_intNMF_input_blocks.rds"))

feature_summary <- tibble(
  modality = c("mRNA", "lncRNA", "smallRNA"),
  features_used = c(ncol(dat$mRNA), ncol(dat$lncRNA), ncol(dat$smallRNA)),
  min_value = c(min(dat$mRNA), min(dat$lncRNA), min(dat$smallRNA)),
  max_value = c(max(dat$mRNA), max(dat$lncRNA), max(dat$smallRNA))
)

write_csv(feature_summary, file.path(outdir, "technical_adjusted_feature_summary.csv"))
print(feature_summary)

cat("\n===== Technical-adjusted intNMF pilot K selection =====\n")

k_range <- 2:6

pdf(file.path(outdir, "technical_adjusted_intNMF_CPI_native_plot.pdf"), width = 6.5, height = 5)

cpi <- nmf.opt.k(
  dat = dat,
  n.runs = 5,
  n.fold = 5,
  k.range = k_range,
  result = TRUE,
  make.plot = TRUE,
  progress = TRUE,
  st.count = 20,
  maxiter = 200,
  wt = rep(1, length(dat))
)

dev.off()

cpi <- as.matrix(cpi)

if (nrow(cpi) != length(k_range) && ncol(cpi) == length(k_range)) {
  cpi <- t(cpi)
}
if (nrow(cpi) != length(k_range)) {
  stop("Unexpected CPI dimensions: ", paste(dim(cpi), collapse = " x "))
}

rownames(cpi) <- paste0("K", k_range)

write_csv(
  as.data.frame(cpi) %>% rownames_to_column("K"),
  file.path(outdir, "technical_adjusted_intNMF_CPI_all_runs.csv")
)

cpi_summary <- tibble(
  K = k_range,
  CPI_mean = rowMeans(cpi, na.rm = TRUE),
  CPI_median = apply(cpi, 1, median, na.rm = TRUE),
  CPI_sd = apply(cpi, 1, sd, na.rm = TRUE),
  CPI_min = apply(cpi, 1, min, na.rm = TRUE),
  CPI_max = apply(cpi, 1, max, na.rm = TRUE)
)

best_k <- cpi_summary$K[which.max(cpi_summary$CPI_mean)]

cpi_summary <- cpi_summary %>%
  mutate(pilot_best_K = K == best_k)

write_csv(
  cpi_summary,
  file.path(outdir, "technical_adjusted_intNMF_CPI_summary.csv")
)

print(cpi_summary)
cat("Technical-adjusted pilot best K:", best_k, "\n")

fit_and_audit <- function(k, prefix) {
  cat("\n===== Fit technical-adjusted intNMF K=", k, "=====\n")
  
  set.seed(42)
  
  fit <- nmf.mnnals(
    dat = dat,
    k = k,
    maxiter = 400,
    st.count = 40,
    n.ini = 30,
    ini.nndsvd = TRUE,
    seed = TRUE,
    wt = rep(1, length(dat))
  )
  
  clusters <- as.integer(fit$clusters)
  
  labels_new <- tibble(
    sample_id = sample_order,
    adjusted_intNMF_cluster = paste0("C", clusters)
  ) %>%
    left_join(meta, by = "sample_id")
  
  write_csv(labels_new, file.path(outdir, paste0(prefix, "_labels.csv")))
  
  consensus <- as.matrix(fit$consensus)
  rownames(consensus) <- sample_order
  colnames(consensus) <- sample_order
  
  write_csv(
    as.data.frame(consensus) %>% rownames_to_column("sample_id"),
    file.path(outdir, paste0(prefix, "_consensus_matrix.csv"))
  )
  
  dist_cons <- as.dist(1 - consensus)
  sil <- silhouette(clusters, dist_cons)
  
  metrics <- tibble(
    K = k,
    mean_silhouette = mean(sil[, "sil_width"]),
    cophenetic_correlation = tryCatch(NMF::cophcor(consensus), error = function(e) NA_real_),
    dispersion = tryCatch(NMF::dispersion(consensus), error = function(e) NA_real_),
    min_cluster_size = min(table(labels_new$adjusted_intNMF_cluster)),
    max_cluster_size = max(table(labels_new$adjusted_intNMF_cluster)),
    cluster_sizes = paste(
      names(table(labels_new$adjusted_intNMF_cluster)),
      as.integer(table(labels_new$adjusted_intNMF_cluster)),
      sep = "=",
      collapse = ";"
    )
  )
  
  write_csv(metrics, file.path(outdir, paste0(prefix, "_metrics.csv")))
  
  qc_join <- qc %>%
    select(sample_id, starts_with("log10_"), ends_with("_detected_features"), ends_with("_zero_fraction"), subgroup) %>%
    left_join(labels_new %>% select(sample_id, adjusted_intNMF_cluster), by = "sample_id")
  
  qc_metrics <- colnames(qc_join)[
    grepl("log10_|detected_features|zero_fraction", colnames(qc_join))
  ]
  
  qc_tests <- lapply(qc_metrics, function(metric) {
    form <- as.formula(paste(metric, "~ adjusted_intNMF_cluster"))
    p <- tryCatch(kruskal.test(form, data = qc_join)$p.value, error = function(e) NA_real_)
    tibble(metric = metric, pvalue = p)
  }) %>%
    bind_rows() %>%
    mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
    arrange(FDR)
  
  write_csv(qc_tests, file.path(outdir, paste0(prefix, "_QC_association_tests.csv")))
  
  cat("\nCluster sizes:\n")
  print(table(labels_new$adjusted_intNMF_cluster))
  cat("\nMetrics:\n")
  print(metrics)
  cat("\nTop QC associations:\n")
  print(head(qc_tests, 10))
  
  sample_plot_order <- labels_new %>%
    arrange(adjusted_intNMF_cluster, subgroup, sample_id) %>%
    pull(sample_id)
  
  ann <- labels_new %>%
    select(sample_id, adjusted_intNMF_cluster, subgroup) %>%
    as.data.frame()
  rownames(ann) <- ann$sample_id
  ann$sample_id <- NULL
  
  pdf(file.path(outdir, paste0(prefix, "_consensus_heatmap.pdf")), width = 7, height = 7)
  pheatmap(
    consensus[sample_plot_order, sample_plot_order],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_row = ann[sample_plot_order, , drop = FALSE],
    annotation_col = ann[sample_plot_order, , drop = FALSE],
    border_color = NA,
    main = paste0("Technical-adjusted intNMF consensus, K=", k)
  )
  dev.off()
  
  saveRDS(fit, file.path(outdir, paste0(prefix, "_fit.rds")))
  
  invisible(list(labels = labels_new, metrics = metrics, qc_tests = qc_tests))
}

fit_best <- fit_and_audit(best_k, paste0("technical_adjusted_bestK", best_k))

if (best_k != 2) {
  fit_k2 <- fit_and_audit(2, "technical_adjusted_K2_for_comparison")
}

cat("\n===== DONE =====\n")
cat("Output directory:", outdir, "\n")
print(list.files(outdir, full.names = TRUE))
