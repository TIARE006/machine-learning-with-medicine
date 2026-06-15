.libPaths(c(normalizePath("Rlibs", winslash = "/", mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages({
  library(IntNMF)
  library(mclust)
  library(readr)
  library(tibble)
})

set.seed(11)

input_path <- "results/fig1_inputs/generated/two_block_intNMF_inputs.rds"
archived_fit_path <- "results/intNMF_reproducible_two_block/final_fit.rds"
outdir <- "results/intNMF_reproducible_two_block/reproduction_validation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(input_path)
archived <- readRDS(archived_fit_path)

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

ari <- adjustedRandIndex(as.integer(archived$clusters), as.integer(fit$clusters))
consensus_max_abs_diff <- max(abs(as.matrix(archived$consensus) - as.matrix(fit$consensus)))

result <- tibble(
  seed = 11,
  adjusted_rand_index = ari,
  consensus_max_abs_diff = consensus_max_abs_diff,
  reference_cluster_sizes = paste(as.integer(table(archived$clusters)), collapse = ";"),
  reproduced_cluster_sizes = paste(as.integer(table(fit$clusters)), collapse = ";")
)

write_csv(result, file.path(outdir, "raw_to_fit_reproduction_metrics.csv"))
saveRDS(fit, file.path(outdir, "reproduced_K4_fit_seed11.rds"))
print(result)

if (ari < 1 || consensus_max_abs_diff > 1e-12) quit(status = 2)
