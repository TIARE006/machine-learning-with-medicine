.libPaths(c(normalizePath("Rlibs", winslash = "/", mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages({
  library(IntNMF)
  library(NMF)
  library(cluster)
  library(mclust)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(digest)
})

input_path <- "results/fig1_inputs/generated/two_block_intNMF_inputs.rds"
meta_path <- "results/fig1_inputs/generated/intNMF_sample_metadata.csv"
feature_spec_path <- "scripts/fig1/01_expression_inputs/locked_intnmf_feature_spec.csv"
outdir <- "results/intNMF_reproducible_two_block"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(input_path)
meta <- read_csv(meta_path, show_col_types = FALSE)

if (!all(vapply(dat, function(x) identical(rownames(x), meta$sample_id), logical(1)))) {
  stop("intNMF blocks do not follow the locked metadata sample order.")
}

package_versions <- tibble(
  package = c("R", "IntNMF", "NMF", "edgeR", "limma", "matrixStats"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    vapply(c("IntNMF", "NMF", "edgeR", "limma", "matrixStats"),
           function(x) as.character(packageVersion(x)), character(1))
  )
)
write_csv(package_versions, file.path(outdir, "package_versions.csv"))

set.seed(20260615)
k_range <- 2:6

pdf(file.path(outdir, "CPI_native_plot.pdf"), width = 6.5, height = 5)
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
  wt = c(1, 1)
)
dev.off()

cpi <- as.matrix(cpi)
if (nrow(cpi) != length(k_range) && ncol(cpi) == length(k_range)) cpi <- t(cpi)
if (nrow(cpi) != length(k_range)) stop("Unexpected CPI dimensions.")
rownames(cpi) <- paste0("K", k_range)

cpi_runs <- as.data.frame(cpi) %>% rownames_to_column("K")
write_csv(cpi_runs, file.path(outdir, "CPI_all_runs.csv"))

cpi_summary <- tibble(
  K = k_range,
  CPI_mean = rowMeans(cpi),
  CPI_median = apply(cpi, 1, median),
  CPI_sd = apply(cpi, 1, sd),
  CPI_q25 = apply(cpi, 1, quantile, probs = 0.25),
  CPI_q75 = apply(cpi, 1, quantile, probs = 0.75),
  CPI_min = apply(cpi, 1, min),
  CPI_max = apply(cpi, 1, max)
) %>%
  mutate(selected_K = K == K[which.max(CPI_mean)])
write_csv(cpi_summary, file.path(outdir, "CPI_summary.csv"))
print(cpi_summary)

best_k <- cpi_summary$K[which.max(cpi_summary$CPI_mean)]
seeds <- c(11, 23, 37, 41, 53, 67, 79, 83, 97, 101,
           113, 127, 139, 149, 163, 173, 181, 193, 211, 227)

labels_matrix <- matrix(
  NA_integer_, nrow = nrow(meta), ncol = length(seeds),
  dimnames = list(meta$sample_id, paste0("seed_", seeds))
)
fits <- vector("list", length(seeds))
run_metrics <- vector("list", length(seeds))

for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  fit <- nmf.mnnals(
    dat = dat, k = best_k, maxiter = 500, st.count = 50,
    n.ini = 20, ini.nndsvd = TRUE, seed = TRUE, wt = c(1, 1)
  )
  labels_matrix[, i] <- as.integer(fit$clusters)
  consensus <- as.matrix(fit$consensus)
  sil <- silhouette(as.integer(fit$clusters), as.dist(1 - consensus))
  run_metrics[[i]] <- tibble(
    seed = seeds[i],
    K = best_k,
    mean_silhouette = mean(sil[, "sil_width"]),
    cophenetic_correlation = NMF::cophcor(consensus),
    dispersion = NMF::dispersion(consensus),
    cluster_sizes = paste(as.integer(table(fit$clusters)), collapse = ";")
  )
  fits[[i]] <- fit
  cat("Completed stable-fit seed", seeds[i], "\n")
}

ari <- outer(
  seq_along(seeds), seq_along(seeds),
  Vectorize(function(i, j) adjustedRandIndex(labels_matrix[, i], labels_matrix[, j]))
)
dimnames(ari) <- list(paste0("seed_", seeds), paste0("seed_", seeds))
mean_ari <- (rowSums(ari) - 1) / (length(seeds) - 1)
medoid_index <- which.max(mean_ari)
final_fit <- fits[[medoid_index]]

labels <- tibble(
  sample_id = meta$sample_id,
  reproducible_two_block_cluster = paste0("C", as.integer(final_fit$clusters))
) %>% left_join(meta, by = "sample_id")

consensus <- as.matrix(final_fit$consensus)
rownames(consensus) <- meta$sample_id
colnames(consensus) <- meta$sample_id

write_csv(bind_rows(run_metrics), file.path(outdir, "seed_run_metrics.csv"))
write_csv(as.data.frame(labels_matrix) %>% rownames_to_column("sample_id"), file.path(outdir, "labels_across_seeds.csv"))
write_csv(as.data.frame(ari) %>% rownames_to_column("run"), file.path(outdir, "pairwise_ARI.csv"))
write_csv(labels, file.path(outdir, "final_labels.csv"))
write_csv(as.data.frame(consensus) %>% rownames_to_column("sample_id"), file.path(outdir, "final_consensus_matrix.csv"))
saveRDS(final_fit, file.path(outdir, "final_fit.rds"))

summary <- tibble(
  selected_K = best_k,
  medoid_seed = seeds[medoid_index],
  pairwise_ARI_mean = mean(ari[upper.tri(ari)]),
  pairwise_ARI_min = min(ari[upper.tri(ari)]),
  input_sha256 = digest(input_path, algo = "sha256", file = TRUE),
  feature_spec_sha256 = digest(feature_spec_path, algo = "sha256", file = TRUE)
)
write_csv(summary, file.path(outdir, "reproducibility_summary.csv"))
print(summary)

manifest_files <- c(
  input_path,
  meta_path,
  feature_spec_path,
  file.path(outdir, "CPI_all_runs.csv"),
  file.path(outdir, "CPI_summary.csv"),
  file.path(outdir, "final_fit.rds"),
  file.path(outdir, "final_labels.csv"),
  file.path(outdir, "final_consensus_matrix.csv"),
  file.path(outdir, "reproducibility_summary.csv")
)
write_csv(
  tibble(
    role = c("model_input", "sample_metadata", "locked_feature_spec", "CPI_runs", "CPI_summary", "final_fit", "raw_cluster_labels", "final_consensus", "stability_summary"),
    path = manifest_files,
    sha256 = vapply(manifest_files, digest, algo = "sha256", file = TRUE, FUN.VALUE = character(1))
  ),
  file.path(outdir, "training_lineage_manifest.csv")
)
