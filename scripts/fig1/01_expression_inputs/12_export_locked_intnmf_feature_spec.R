.libPaths(c(normalizePath("Rlibs", winslash = "/", mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages({
  library(tibble)
  library(readr)
  library(dplyr)
})

fit_path <- "results/intNMF_adjusted_two_block_K4_final/two_block_final_K4_fit.rds"
out_path <- "scripts/fig1/01_expression_inputs/locked_intnmf_feature_spec.csv"

fit <- readRDS(fit_path)
h1 <- colnames(fit$H[[1]])
h2 <- colnames(fit$H[[2]])

spec <- bind_rows(
  tibble(
    block = "RNA_lncRNA",
    modality = sub("__.*$", "", h1),
    feature = sub("^[^_]+__", "", h1),
    model_name = h1,
    position = seq_along(h1)
  ),
  tibble(
    block = "smallRNA",
    modality = "smallRNA",
    feature = h2,
    model_name = h2,
    position = seq_along(h2)
  )
)

write_csv(spec, out_path)
cat("Wrote locked intNMF feature specification:", out_path, "\n")
