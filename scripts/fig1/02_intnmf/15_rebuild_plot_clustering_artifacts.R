.libPaths(c(normalizePath("Rlibs", winslash = "/", mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(digest)
})

fit_path <- "results/intNMF_adjusted_two_block_K4_final/two_block_final_K4_fit.rds"
cpi_runs_path <- "results/intNMF_strict84_technical_adjusted_final/technical_adjusted_final_CPI_all_runs.csv"
pheno_path <- "data/raw/RNA_seq/GSE254877_pheno.csv"
k4_outdir <- "results/intNMF_adjusted_two_block_K4_final"
cpi_outdir <- "results/intNMF_strict84_technical_adjusted_final"
manifest_path <- "results/fig1_inputs/clustering_plot_artifact_lineage_manifest.csv"

stopifnot(file.exists(fit_path), file.exists(cpi_runs_path), file.exists(pheno_path))

fit <- readRDS(fit_path)
sample_ids <- names(fit$clusters)
if (length(sample_ids) != 84 || anyDuplicated(sample_ids) > 0) {
  stop("Archived K4 fit does not contain 84 unique named clusters.")
}

pheno <- read_csv(pheno_path, show_col_types = FALSE, name_repair = "unique") %>%
  transmute(
    sample_id = sub(",.*$", "", as.character(title)),
    geo_accession = as.character(geo_accession),
    title = as.character(title),
    subgroup = trimws(sub("^[^,]+,", "", as.character(title)))
  ) %>%
  filter(sample_id %in% sample_ids)

labels <- tibble(
  sample_id = sample_ids,
  adjusted_two_block_K4 = paste0("C", as.integer(fit$clusters))
) %>%
  left_join(pheno, by = "sample_id") %>%
  arrange(adjusted_two_block_K4, sample_id)

if (nrow(labels) != 84 || any(is.na(labels$subgroup))) {
  stop("Failed to attach GEO phenotype to archived K4 fit samples.")
}

consensus <- as.matrix(fit$consensus)
rownames(consensus) <- sample_ids
colnames(consensus) <- sample_ids

label_path <- file.path(k4_outdir, "two_block_final_K4_labels.csv")
consensus_path <- file.path(k4_outdir, "two_block_final_K4_consensus_matrix.csv")
write_csv(labels, label_path)
write_csv(as.data.frame(consensus) %>% rownames_to_column("sample_id"), consensus_path)

cpi_runs <- read_csv(cpi_runs_path, show_col_types = FALSE) %>%
  mutate(K = as.integer(sub("^K", "", K)))

cpi_summary <- cpi_runs %>%
  pivot_longer(-K, names_to = "run", values_to = "CPI") %>%
  group_by(K) %>%
  summarise(
    CPI_mean = mean(CPI),
    CPI_median = median(CPI),
    CPI_sd = sd(CPI),
    CPI_q25 = quantile(CPI, 0.25),
    CPI_q75 = quantile(CPI, 0.75),
    CPI_min = min(CPI),
    CPI_max = max(CPI),
    .groups = "drop"
  ) %>%
  mutate(final_best_K = K == K[which.max(CPI_mean)])

cpi_summary_path <- file.path(cpi_outdir, "technical_adjusted_final_CPI_summary.csv")
write_csv(cpi_summary, cpi_summary_path)

files <- c(fit_path, cpi_runs_path, pheno_path, label_path, consensus_path, cpi_summary_path)
write_csv(
  tibble(
    role = c("archived_K4_fit", "archived_CPI_runs", "GEO_phenotype", "regenerated_K4_labels", "regenerated_K4_consensus", "regenerated_CPI_summary"),
    path = files,
    sha256 = vapply(files, digest, algo = "sha256", file = TRUE, FUN.VALUE = character(1))
  ),
  manifest_path
)

cat("Rebuilt exact plotting artifacts from archived intNMF model/run objects.\n")
