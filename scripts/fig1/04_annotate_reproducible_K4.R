.libPaths(c(normalizePath("Rlibs", winslash = "/", mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(digest)
})

label_path <- "results/intNMF_reproducible_two_block/final_labels.csv"
score_path <- "results/fig1_final/all84_pyro_module_Euclidean_K2_main/all84_pyro_module_scores_by_sample.csv"
validation_path <- "results/fig1_final/all84_pyro_module_Euclidean_K2_main/all84_K2_validation_scores_by_sample.csv"
outdir <- "results/intNMF_reproducible_two_block"

labels <- read_csv(label_path, show_col_types = FALSE)
scores <- read_csv(score_path, show_col_types = FALSE)
validation <- read_csv(validation_path, show_col_types = FALSE)

audit <- labels %>%
  left_join(scores %>% select(sample_id, integrated_pyro_module_score), by = "sample_id") %>%
  left_join(validation %>% select(sample_id, Heldout_pyroptosis), by = "sample_id") %>%
  group_by(raw_cluster = reproducible_two_block_cluster) %>%
  summarise(
    n = n(),
    mean_integrated_pyro = mean(integrated_pyro_module_score),
    mean_heldout_pyro = mean(Heldout_pyroptosis),
    .groups = "drop"
  ) %>%
  arrange(mean_integrated_pyro)

if (nrow(audit) != 4) stop("Expected four reproducible intNMF clusters.")

# The two clusters with the highest locked pyroptosis-module means are PA.
# Renumbering is deterministic and keeps the manuscript convention PA=C2+C4.
mapping <- audit %>%
  mutate(
    PA_PQ = c("PQ", "PQ", "PA", "PA"),
    adjusted_two_block_K4 = c("C1", "C3", "C4", "C2")
  )

annotated <- labels %>%
  left_join(
    mapping %>% select(raw_cluster, adjusted_two_block_K4, PA_PQ),
    by = c("reproducible_two_block_cluster" = "raw_cluster")
  ) %>%
  select(sample_id, adjusted_two_block_K4, geo_accession, title, subgroup, PA_PQ, reproducible_two_block_cluster)

write_csv(mapping, file.path(outdir, "K4_PA_PQ_annotation_audit.csv"))
write_csv(annotated, file.path(outdir, "two_block_final_K4_labels.csv"))

manifest_files <- c(
  label_path,
  score_path,
  validation_path,
  file.path(outdir, "K4_PA_PQ_annotation_audit.csv"),
  file.path(outdir, "two_block_final_K4_labels.csv")
)
write_csv(
  tibble(
    role = c("raw_K4_labels", "pyroptosis_module_scores", "heldout_validation_scores", "annotation_audit", "annotated_K4_labels"),
    path = manifest_files,
    sha256 = vapply(manifest_files, digest, algo = "sha256", file = TRUE, FUN.VALUE = character(1))
  ),
  file.path(outdir, "annotation_lineage_manifest.csv")
)
print(mapping)
