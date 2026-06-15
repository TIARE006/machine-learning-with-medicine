options(stringsAsFactors = FALSE, timeout = 1800)

.libPaths(c(
  normalizePath("Rlibs", winslash = "/", mustWork = TRUE),
  .libPaths()
))

local_lib <- normalizePath("Rlibs", winslash = "/", mustWork = TRUE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_need <- c(
  "readr", "dplyr", "tidyr", "tibble", "ggplot2",
  "pheatmap", "magick", "circlize", "RColorBrewer"
)

for (p in cran_need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      lib = local_lib,
      dependencies = c("Depends", "Imports", "LinkingTo")
    )
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = local_lib)
}

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install(
    "ComplexHeatmap",
    lib = local_lib,
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(magick)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

# ============================================================
# 0. Paths
# ============================================================

outdir <- "results/fig1_final/intNMF_K4_PA_PQ_final_main_figure"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

label_path <- "results/intNMF_reproducible_two_block/two_block_final_K4_labels.csv"
consensus_path <- "results/intNMF_reproducible_two_block/final_consensus_matrix.csv"

pcd_z_path <- "results/fig1_inputs/generated/GSE254877_14PCD_ssGSEA_zscores.csv"

module_path <- paste0(
  "results/fig1_final/all84_pyro_module_Euclidean_K2_main/",
  "all84_pyro_module_scores_by_sample.csv"
)

validation_path <- paste0(
  "results/fig1_final/all84_pyro_module_Euclidean_K2_main/",
  "all84_K2_validation_scores_by_sample.csv"
)

metadata_path <- paste0(
  "results/fig1_final/all84_pyro_module_Euclidean_K2_main/",
  "all84_sample_metadata_used.csv"
)

marker_z_path <- "results/fig1_inputs/generated/GSE254877_pyroptosis_marker_gene_zmatrix.csv"

required_files <- c(
  label_path,
  consensus_path,
  pcd_z_path,
  module_path,
  validation_path,
  metadata_path,
  marker_z_path
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ============================================================
# 1. Labels: stable intNMF K4, PA=C2+C4, PQ=C1+C3
# ============================================================

labels_raw <- read_csv(label_path, show_col_types = FALSE)

if (!all(c("sample_id", "adjusted_two_block_K4") %in% names(labels_raw))) {
  stop("Label file must contain sample_id and adjusted_two_block_K4.")
}

labels <- labels_raw %>%
  transmute(
    sample_id = as.character(sample_id),
    intNMF_K4 = as.character(adjusted_two_block_K4)
  ) %>%
  mutate(
    intNMF_K4 = factor(intNMF_K4, levels = c("C1", "C2", "C3", "C4")),
    PA_PQ = case_when(
      intNMF_K4 %in% c("C2", "C4") ~ "PA",
      intNMF_K4 %in% c("C1", "C3") ~ "PQ",
      TRUE ~ NA_character_
    ),
    PA_PQ_plot = factor(PA_PQ, levels = c("PA", "PQ")),
    PA_PQ_lm = factor(PA_PQ, levels = c("PQ", "PA"))
  )

if (nrow(labels) != 84 || anyDuplicated(labels$sample_id) > 0) {
  stop("Expected 84 unique intNMF K4 labels.")
}

cat("\n===== intNMF K4 / PA-PQ sizes =====\n")
print(table(labels$intNMF_K4))
print(table(labels$PA_PQ_plot))

state_sizes <- table(labels$PA_PQ_plot)
pa_n <- unname(state_sizes["PA"])
pq_n <- unname(state_sizes["PQ"])

# ============================================================
# 2. Read scores and metadata
# ============================================================

pcd_programs <- c(
  "Pyroptosis",
  "Ferroptosis",
  "Necroptosis",
  "Alkaliptosis",
  "Anoikis",
  "NETosis",
  "Lysosome_dependent_cell_death",
  "Autophagy",
  "Parthanatos",
  "Entosis",
  "Apoptosis",
  "Disulfidptosis",
  "Cuproptosis",
  "Oxeiptosis"
)

pcd_z <- read_csv(pcd_z_path, show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

pcd_cols <- intersect(pcd_programs, names(pcd_z))

if (length(pcd_cols) < 10) {
  stop("Too few PCD program columns found in pcd_z file.")
}

pcd_z <- pcd_z %>%
  dplyr::select(sample_id, all_of(pcd_cols)) %>%
  mutate(across(all_of(pcd_cols), ~ suppressWarnings(as.numeric(.x))))

module_scores <- read_csv(module_path, show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

validation_scores <- read_csv(validation_path, show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

metadata <- read_csv(metadata_path, show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

analysis_df <- labels %>%
  left_join(
    module_scores %>%
      dplyr::select(
        sample_id,
        Inflammasome_sensor_adapter,
        Inflammatory_caspase_axis,
        NOD_gasdermin_related_axis,
        integrated_pyro_module_score
      ),
    by = "sample_id"
  ) %>%
  left_join(
    validation_scores %>%
      dplyr::select(
        sample_id,
        Heldout_pyroptosis,
        Inflammation,
        Myeloid_infiltration,
        Muscle_atrophy
      ),
    by = "sample_id"
  ) %>%
  left_join(
    metadata %>%
      dplyr::select(sample_id, sex, subgroup),
    by = "sample_id"
  ) %>%
  left_join(
    pcd_z,
    by = "sample_id"
  ) %>%
  mutate(
    sex = factor(as.character(sex)),
    subgroup = factor(as.character(subgroup), levels = c("Colorectal", "Pancreas"))
  )

if (nrow(analysis_df) != 84 || anyDuplicated(analysis_df$sample_id) > 0) {
  stop("Merged analysis table does not contain 84 unique samples.")
}

write_csv(
  analysis_df,
  file.path(outdir, "intNMF_K4_PA_PQ_analysis_table.csv")
)

# ============================================================
# Utility functions
# ============================================================

zscore_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

clip_mat <- function(mat, limit = 2) {
  mat[mat > limit] <- limit
  mat[mat < -limit] <- -limit
  mat
}

star_from_fdr <- function(fdr) {
  ifelse(
    fdr < 0.001, "****",
    ifelse(
      fdr < 0.01, "***",
      ifelse(
        fdr < 0.05, "**",
        ifelse(fdr < 0.10, "*", "ns")
      )
    )
  )
}

# Color settings
k4_cols <- c(
  C1 = "#F8766D",
  C2 = "#00BA38",
  C3 = "#619CFF",
  C4 = "#C77CFF"
)

pa_cols <- c(
  PA = "#E64B35",
  PQ = "#4DBBD5"
)

subgroup_cols <- c(
  Colorectal = "#7CAE00",
  Pancreas = "#00BFC4"
)

heat_col <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "#F7F7F7", "#B2182B")
)

# ============================================================
# Panel A: 14-PCD Spearman correlation matrix
# ============================================================

pcd_mat_for_cor <- analysis_df %>%
  dplyr::select(all_of(pcd_cols)) %>%
  as.matrix()

cor_mat <- cor(
  pcd_mat_for_cor,
  use = "pairwise.complete.obs",
  method = "spearman"
)

A_png <- file.path(outdir, "Fig1A_14PCD_spearman_correlation.png")
A_pdf <- file.path(outdir, "Fig1A_14PCD_spearman_correlation.pdf")

png(A_png, width = 2200, height = 1900, res = 300)

pheatmap(
  cor_mat,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  display_numbers = matrix(sprintf("%.2f", cor_mat), nrow = nrow(cor_mat)),
  number_color = "black",
  fontsize_number = 5.5,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = "white",
  main = "Spearman correlation among 14 PCD programs"
)

dev.off()

pdf(A_pdf, width = 7.2, height = 6.2, useDingbats = FALSE)

pheatmap(
  cor_mat,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  display_numbers = matrix(sprintf("%.2f", cor_mat), nrow = nrow(cor_mat)),
  number_color = "black",
  fontsize_number = 5.5,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = "white",
  main = "Spearman correlation among 14 PCD programs"
)

dev.off()

# ============================================================
# Panel B top: K robustness / K selection
# ============================================================

metric_files <- list.files(
  "results",
  pattern = "two_block_final_K[0-9]+_metrics\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

metric_tbl <- list()

for (f in metric_files) {
  dat <- tryCatch(read_csv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) next

  if (!all(c("K", "cophenetic_correlation") %in% names(dat))) next

  metric_tbl[[length(metric_tbl) + 1]] <- dat %>%
    transmute(
      K = as.integer(K),
      cophenetic_correlation = as.numeric(cophenetic_correlation),
      mean_silhouette = if ("mean_silhouette" %in% names(dat)) {
        as.numeric(mean_silhouette)
      } else {
        NA_real_
      },
      minimum_cluster_size = if ("minimum_cluster_size" %in% names(dat)) {
        as.numeric(minimum_cluster_size)
      } else {
        NA_real_
      }
    )
}

metric_tbl <- bind_rows(metric_tbl) %>%
  distinct(K, .keep_all = TRUE) %>%
  arrange(K)

Btop_png <- file.path(outdir, "Fig1B_top_K_robustness_curve.png")

if (nrow(metric_tbl) >= 2) {

  pBtop <- ggplot(
    metric_tbl,
    aes(x = K, y = cophenetic_correlation)
  ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.6) +
    geom_vline(xintercept = 4, linetype = "dashed", linewidth = 0.5) +
    theme_bw(base_size = 10) +
    scale_x_continuous(breaks = sort(unique(metric_tbl$K))) +
    labs(
      title = "Robustness of two-block adjusted intNMF supports K=4",
      x = "Number of subtypes (K)",
      y = "Cophenetic correlation"
    )

} else {

  k4_metric <- read_csv(
    "results/intNMF_adjusted_two_block_K4_final/two_block_final_K4_metrics.csv",
    show_col_types = FALSE
  )

  pBtop <- ggplot(
    k4_metric,
    aes(x = K, y = cophenetic_correlation)
  ) +
    geom_point(size = 3) +
    geom_vline(xintercept = 4, linetype = "dashed", linewidth = 0.5) +
    annotate(
      "text",
      x = 4,
      y = k4_metric$cophenetic_correlation,
      label = paste0("K=4, coph.=", sprintf("%.3f", k4_metric$cophenetic_correlation)),
      vjust = -1,
      size = 3.5
    ) +
    theme_bw(base_size = 10) +
    scale_x_continuous(breaks = 4, limits = c(2, 6)) +
    ylim(0, 1) +
    labs(
      title = "Stable two-block adjusted intNMF K=4",
      x = "Number of subtypes (K)",
      y = "Cophenetic correlation"
    )
}

ggsave(Btop_png, pBtop, width = 6.2, height = 2.2, dpi = 300)

# ============================================================
# Panel B bottom: true K4 consensus matrix
# ============================================================

cons_raw <- read.csv(
  consensus_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_id_col <- names(cons_raw)[1]

consensus <- as.matrix(cons_raw[, -1, drop = FALSE])
storage.mode(consensus) <- "double"
rownames(consensus) <- as.character(cons_raw[[sample_id_col]])

if (is.null(colnames(consensus)) || !setequal(rownames(consensus), colnames(consensus))) {
  stop("Consensus matrix row and column IDs do not match.")
}

order_k4 <- analysis_df %>%
  arrange(intNMF_K4, subgroup, sample_id) %>%
  pull(sample_id)

ann_B <- analysis_df %>%
  dplyr::select(sample_id, intNMF_K4, PA_PQ_plot, subgroup) %>%
  mutate(
    intNMF_K4 = factor(intNMF_K4, levels = c("C1", "C2", "C3", "C4")),
    PA_PQ_plot = factor(PA_PQ_plot, levels = c("PA", "PQ"))
  ) %>%
  as.data.frame()

rownames(ann_B) <- ann_B$sample_id
ann_B$sample_id <- NULL
names(ann_B) <- c("K4_cluster", "PA_PQ", "subgroup")

ann_cols_B <- list(
  K4_cluster = k4_cols,
  PA_PQ = pa_cols,
  subgroup = subgroup_cols
)

Bbottom_png <- file.path(outdir, "Fig1B_bottom_true_K4_consensus.png")
Bbottom_pdf <- file.path(outdir, "Fig1B_bottom_true_K4_consensus.pdf")

png(Bbottom_png, width = 2200, height = 2100, res = 300)

pheatmap(
  consensus[order_k4, order_k4, drop = FALSE],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = ann_B[order_k4, , drop = FALSE],
  annotation_col = ann_B[order_k4, , drop = FALSE],
  annotation_colors = ann_cols_B,
  border_color = NA,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  breaks = seq(0, 1, length.out = 101),
  main = "True K4 consensus matrix (84 samples)"
)

dev.off()

pdf(Bbottom_pdf, width = 7.2, height = 7.0, useDingbats = FALSE)

pheatmap(
  consensus[order_k4, order_k4, drop = FALSE],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = ann_B[order_k4, , drop = FALSE],
  annotation_col = ann_B[order_k4, , drop = FALSE],
  annotation_colors = ann_cols_B,
  border_color = NA,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  breaks = seq(0, 1, length.out = 101),
  main = "True K4 consensus matrix (84 samples)"
)

dev.off()

# ============================================================
# Panel C: MDS from consensus distance
# ============================================================

cmd <- cmdscale(
  as.dist(1 - consensus[order_k4, order_k4]),
  k = 2,
  eig = TRUE
)

eig <- cmd$eig
positive_eig <- eig[eig > 0]
var_pct <- positive_eig[1:2] / sum(positive_eig) * 100

mds_df <- tibble(
  sample_id = order_k4,
  MDS1 = cmd$points[, 1],
  MDS2 = cmd$points[, 2]
) %>%
  left_join(
    analysis_df %>%
      dplyr::select(sample_id, intNMF_K4, subgroup, PA_PQ_plot),
    by = "sample_id"
  )

C_png <- file.path(outdir, "Fig1C_intNMF_K4_MDS.png")
C_pdf <- file.path(outdir, "Fig1C_intNMF_K4_MDS.pdf")

pC <- ggplot(
  mds_df,
  aes(
    x = MDS1,
    y = MDS2,
    color = intNMF_K4,
    shape = subgroup
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey75") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey75") +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = k4_cols, name = "K4 subtype") +
  theme_bw(base_size = 12) +
  labs(
    title = "Low-dimensional separation of K4 subtypes (MDS)",
    x = paste0("MDS1 (", sprintf("%.1f", var_pct[1]), "% variance)"),
    y = paste0("MDS2 (", sprintf("%.1f", var_pct[2]), "% variance)"),
    shape = "Subgroup"
  ) +
  annotate(
    "label",
    x = max(mds_df$MDS1, na.rm = TRUE),
    y = min(mds_df$MDS2, na.rm = TRUE),
    hjust = 1,
    vjust = 0,
    label = "Candidate PA clusters:\nC2 and C4",
    size = 3.2
  )

ggsave(C_png, pC, width = 7.2, height = 5.8, dpi = 300)
ggsave(C_pdf, pC, width = 7.2, height = 5.8)

# ============================================================
# Panel D: 14 PCD scores across K4
# ============================================================

pcd_long <- analysis_df %>%
  dplyr::select(sample_id, intNMF_K4, all_of(pcd_cols)) %>%
  pivot_longer(
    cols = all_of(pcd_cols),
    names_to = "PCD_program",
    values_to = "score"
  ) %>%
  mutate(
    PCD_program = factor(PCD_program, levels = pcd_programs[pcd_programs %in% pcd_cols])
  )

D_png <- file.path(outdir, "Fig1D_14PCD_across_intNMF_K4.png")
D_pdf <- file.path(outdir, "Fig1D_14PCD_across_intNMF_K4.pdf")

pD <- ggplot(
  pcd_long,
  aes(x = intNMF_K4, y = score, fill = intNMF_K4)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.78) +
  geom_jitter(width = 0.12, size = 0.65, alpha = 0.45) +
  facet_wrap(~ PCD_program, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = k4_cols, name = "K4 subtype") +
  theme_bw(base_size = 9) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    title = "Programmed cell-death (PCD) scores across K4 subtypes",
    y = "Standardized ssGSEA score"
  )

ggsave(D_png, pD, width = 10.5, height = 6.8, dpi = 300)
ggsave(D_pdf, pD, width = 10.5, height = 6.8)

# ============================================================
# Panel E: PA vs PQ 14-PCD effect forest
# ============================================================

pcd_stats <- list()

for (program in pcd_cols) {

  dat <- analysis_df %>%
    transmute(
      score = .data[[program]],
      PA_PQ_lm,
      sex,
      subgroup
    ) %>%
    filter(!is.na(score), !is.na(PA_PQ_lm), !is.na(sex), !is.na(subgroup))

  fit <- lm(score ~ sex + subgroup + PA_PQ_lm, data = dat)
  sm <- summary(fit)$coefficients

  coef_name <- grep("PA_PQ_lmPA", rownames(sm), value = TRUE)

  pcd_stats[[length(pcd_stats) + 1]] <- tibble(
    PCD_program = program,
    beta_PA_vs_PQ = sm[coef_name, "Estimate"],
    SE = sm[coef_name, "Std. Error"],
    P_value = sm[coef_name, "Pr(>|t|)"]
  )
}

pcd_stats <- bind_rows(pcd_stats) %>%
  mutate(
    FDR = p.adjust(P_value, method = "BH"),
    CI_low = beta_PA_vs_PQ - 1.96 * SE,
    CI_high = beta_PA_vs_PQ + 1.96 * SE
  ) %>%
  arrange(desc(beta_PA_vs_PQ))

write_csv(
  pcd_stats,
  file.path(outdir, "PA_C2C4_vs_PQ_C1C3_14PCD_effect_statistics.csv")
)

pcd_stats_plot <- pcd_stats %>%
  arrange(beta_PA_vs_PQ) %>%
  mutate(
    PCD_program = factor(PCD_program, levels = PCD_program),
    FDR_label = ifelse(
      FDR < 0.001,
      "FDR<0.001",
      paste0("FDR=", signif(FDR, 2))
    )
  )

E_png <- file.path(outdir, "Fig1E_PA_vs_PQ_14PCD_effect_forest.png")
E_pdf <- file.path(outdir, "Fig1E_PA_vs_PQ_14PCD_effect_forest.pdf")

pE <- ggplot(
  pcd_stats_plot,
  aes(x = beta_PA_vs_PQ, y = PCD_program)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_segment(
    aes(x = CI_low, xend = CI_high, y = PCD_program, yend = PCD_program),
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  geom_text(aes(label = FDR_label), hjust = -0.05, size = 3.0) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.35))) +
  theme_bw(base_size = 11) +
  labs(
    title = "Adjusted effects: PA (C2+C4) vs PQ (C1+C3)",
    x = "Adjusted difference in standardized PCD score\n(PA - PQ)",
    y = NULL
  )

ggsave(E_png, pE, width = 7.2, height = 6.2, dpi = 300)
ggsave(E_pdf, pE, width = 7.2, height = 6.2)

# ============================================================
# Panel F: PA vs PQ confirmatory pyroptosis signatures
# ============================================================

sig_features <- c(
  integrated_pyro_module_score = "Integrated pyro module score",
  Heldout_pyroptosis = "Held-out pyroptosis signature",
  Inflammasome_sensor_adapter = "Inflammasome / sensor-adapter axis",
  Inflammatory_caspase_axis = "Inflammatory caspase axis"
)

sig_long <- analysis_df %>%
  dplyr::select(sample_id, PA_PQ_plot, PA_PQ_lm, sex, subgroup, all_of(names(sig_features))) %>%
  pivot_longer(
    cols = all_of(names(sig_features)),
    names_to = "feature",
    values_to = "score"
  ) %>%
  mutate(
    feature_label = sig_features[feature],
    feature_label = factor(feature_label, levels = sig_features)
  )

sig_stats <- list()

for (feature_name in names(sig_features)) {

  dat <- analysis_df %>%
    transmute(
      score = .data[[feature_name]],
      PA_PQ_lm,
      sex,
      subgroup
    ) %>%
    filter(!is.na(score), !is.na(PA_PQ_lm), !is.na(sex), !is.na(subgroup))

  fit <- lm(score ~ sex + subgroup + PA_PQ_lm, data = dat)
  sm <- summary(fit)$coefficients
  coef_name <- grep("PA_PQ_lmPA", rownames(sm), value = TRUE)

  sig_stats[[length(sig_stats) + 1]] <- tibble(
    feature = feature_name,
    feature_label = sig_features[[feature_name]],
    beta = sm[coef_name, "Estimate"],
    SE = sm[coef_name, "Std. Error"],
    P_value = sm[coef_name, "Pr(>|t|)"]
  )
}

sig_stats <- bind_rows(sig_stats) %>%
  mutate(
    FDR = p.adjust(P_value, method = "BH"),
    star = star_from_fdr(FDR)
  )

write_csv(
  sig_stats,
  file.path(outdir, "PA_C2C4_vs_PQ_C1C3_confirmatory_signature_statistics.csv")
)

F_png <- file.path(outdir, "Fig1F_PA_vs_PQ_confirmatory_pyroptosis_signatures.png")
F_pdf <- file.path(outdir, "Fig1F_PA_vs_PQ_confirmatory_pyroptosis_signatures.pdf")

pF <- ggplot(
  sig_long,
  aes(x = PA_PQ_plot, y = score, fill = PA_PQ_plot)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.78) +
  geom_jitter(width = 0.11, size = 0.75, alpha = 0.55) +
  facet_wrap(~ feature_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(
    values = pa_cols,
    labels = c(
      paste0("PA (C2+C4)\n(n=", pa_n, ")"),
      paste0("PQ (C1+C3)\n(n=", pq_n, ")")
    ),
    name = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank()
  ) +
  labs(
    title = "Confirmatory pyroptosis signatures: PA vs PQ",
    y = "Standardized ssGSEA score"
  )

ggsave(F_png, pF, width = 7.2, height = 5.8, dpi = 300)
ggsave(F_pdf, pF, width = 7.2, height = 5.8)

# ============================================================
# Panel G: Integrated heatmap across K4
# ============================================================

g_features <- c(
  "integrated_pyro_module_score",
  "Heldout_pyroptosis",
  "Inflammasome_sensor_adapter",
  "Inflammatory_caspase_axis",
  "NOD_gasdermin_related_axis",
  pcd_cols
)

g_mat <- analysis_df %>%
  arrange(intNMF_K4, subgroup, desc(integrated_pyro_module_score), sample_id) %>%
  dplyr::select(sample_id, all_of(g_features)) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

g_mat <- t(g_mat)
g_mat <- clip_mat(zscore_rows(g_mat), 2)

g_row_names <- rownames(g_mat)
g_row_names[g_row_names == "integrated_pyro_module_score"] <- "Integrated pyro module"
g_row_names[g_row_names == "Heldout_pyroptosis"] <- "Held-out pyroptosis"
g_row_names[g_row_names == "Inflammasome_sensor_adapter"] <- "Inflammasome / sensor-adapter"
g_row_names[g_row_names == "Inflammatory_caspase_axis"] <- "Inflammatory caspase axis"
g_row_names[g_row_names == "NOD_gasdermin_related_axis"] <- "NOD / gasdermin axis"
rownames(g_mat) <- g_row_names

g_col_order <- colnames(g_mat)

g_ann <- analysis_df %>%
  dplyr::select(sample_id, intNMF_K4, PA_PQ_plot, subgroup) %>%
  as.data.frame()

rownames(g_ann) <- g_ann$sample_id
g_ann <- g_ann[g_col_order, , drop = FALSE]

g_col_split <- factor(
  as.character(g_ann$intNMF_K4),
  levels = c("C1", "C2", "C3", "C4")
)

g_top <- HeatmapAnnotation(
  K4_cluster = g_ann$intNMF_K4,
  PA_PQ = g_ann$PA_PQ_plot,
  subgroup = g_ann$subgroup,
  col = list(
    K4_cluster = k4_cols,
    PA_PQ = pa_cols,
    subgroup = subgroup_cols
  ),
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  simple_anno_size = unit(3.5, "mm")
)

g_row_split <- factor(
  c(
    rep("Confirmatory pyroptosis signatures", 5),
    rep("PCD programs (14)", length(pcd_cols))
  ),
  levels = c("Confirmatory pyroptosis signatures", "PCD programs (14)")
)

G_png <- file.path(outdir, "Fig1G_integrated_K4_PA_PQ_heatmap.png")
G_pdf <- file.path(outdir, "Fig1G_integrated_K4_PA_PQ_heatmap.pdf")

make_G <- function() {
  Heatmap(
    g_mat,
    name = "z-score",
    col = heat_col,
    top_annotation = g_top,
    column_split = g_col_split,
    row_split = g_row_split,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 7.2),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_title_gp = gpar(fontsize = 8.5, fontface = "bold"),
    column_gap = unit(2.5, "mm"),
    row_gap = unit(3, "mm"),
    use_raster = FALSE,
    heatmap_legend_param = list(
      title = "z-score",
      at = c(-2, 0, 2)
    )
  )
}

pdf(G_pdf, width = 11, height = 7.0, useDingbats = FALSE)
draw(make_G(), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(G_png, width = 3600, height = 2300, res = 300)
draw(make_G(), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# ============================================================
# Panel H: Core pyroptosis-gene heatmap PA vs PQ
# ============================================================

marker_z <- read_csv(marker_z_path, show_col_types = FALSE) %>%
  column_to_rownames("gene") %>%
  as.matrix()

priority_genes <- c(
  "GSDMD",
  "NOD1",
  "PYCARD",
  "GSDME",
  "NLRP1",
  "GSDMB",
  "IL18",
  "CASP4",
  "CASP1",
  "NLRP3",
  "GSDMC",
  "CASP5",
  "NLRC4",
  "NOD2",
  "IL1B",
  "AIM2"
)

gene_order <- priority_genes[priority_genes %in% rownames(marker_z)]

h_col_order <- analysis_df %>%
  arrange(PA_PQ_plot, subgroup, intNMF_K4, sample_id) %>%
  pull(sample_id)

missing_marker_samples <- setdiff(h_col_order, colnames(marker_z))

if (length(missing_marker_samples) > 0) {
  stop("Marker z-matrix is missing samples: ", paste(missing_marker_samples, collapse = ", "))
}

h_mat <- marker_z[gene_order, h_col_order, drop = FALSE]
h_mat <- clip_mat(h_mat, 2)

h_ann <- analysis_df %>%
  dplyr::select(sample_id, PA_PQ_plot, subgroup) %>%
  as.data.frame()

rownames(h_ann) <- h_ann$sample_id
h_ann <- h_ann[h_col_order, , drop = FALSE]

h_top <- HeatmapAnnotation(
  PA_PQ = h_ann$PA_PQ_plot,
  subgroup = h_ann$subgroup,
  col = list(
    PA_PQ = pa_cols,
    subgroup = subgroup_cols
  ),
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  simple_anno_size = unit(3.5, "mm")
)

h_col_split <- factor(
  as.character(h_ann$PA_PQ_plot),
  levels = c("PA", "PQ")
)

H_png <- file.path(outdir, "Fig1H_core_pyroptosis_genes_PA_vs_PQ.png")
H_pdf <- file.path(outdir, "Fig1H_core_pyroptosis_genes_PA_vs_PQ.pdf")

make_H <- function() {
  Heatmap(
    h_mat,
    name = "z-score",
    col = heat_col,
    top_annotation = h_top,
    column_split = h_col_split,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8.5),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_gap = unit(3, "mm"),
    use_raster = FALSE,
    heatmap_legend_param = list(
      title = "z-score",
      at = c(-2, 0, 2)
    )
  )
}

pdf(H_pdf, width = 8.2, height = 6.0, useDingbats = FALSE)
draw(make_H(), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(H_png, width = 2700, height = 2000, res = 300)
draw(make_H(), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# ============================================================
# Assemble full figure
# ============================================================

read_clean <- function(path) {
  img <- image_read(path)
  img <- image_trim(img)
  img <- image_background(img, color = "white", flatten = TRUE)
  image_border(img, color = "white", geometry = "20x20")
}

fit_box <- function(img, width, height) {
  img <- image_resize(img, paste0(width, "x", height))
  image_extent(
    img,
    geometry = paste0(width, "x", height),
    gravity = "center",
    color = "white"
  )
}

add_panel_label <- function(img, panel_label) {
  image_annotate(
    img,
    text = panel_label,
    gravity = "northwest",
    location = "+18+10",
    size = 78,
    weight = 700,
    color = "black"
  )
}

single_w <- 1600
single_h <- 1120
bottom_h <- 1450

panel_A <- read_clean(A_png) %>%
  fit_box(single_w, single_h) %>%
  add_panel_label("A")

panel_B_top <- read_clean(Btop_png) %>%
  fit_box(single_w, 330)

panel_B_bottom <- read_clean(Bbottom_png) %>%
  fit_box(single_w, 790)

panel_B <- image_append(
  c(panel_B_top, panel_B_bottom),
  stack = TRUE
) %>%
  fit_box(single_w, single_h) %>%
  add_panel_label("B")

panel_C <- read_clean(C_png) %>%
  fit_box(single_w, single_h) %>%
  add_panel_label("C")

panel_D <- read_clean(D_png) %>%
  fit_box(single_w, single_h) %>%
  add_panel_label("D")

panel_E <- read_clean(E_png) %>%
  fit_box(single_w, single_h) %>%
  add_panel_label("E")

panel_F <- read_clean(F_png) %>%
  fit_box(single_w, single_h) %>%
  add_panel_label("F")

panel_G <- read_clean(G_png) %>%
  fit_box(single_w * 2, bottom_h) %>%
  add_panel_label("G")

panel_H <- read_clean(H_png) %>%
  fit_box(single_w, bottom_h) %>%
  add_panel_label("H")

row_1 <- image_append(c(panel_A, panel_B, panel_C), stack = FALSE)
row_2 <- image_append(c(panel_D, panel_E, panel_F), stack = FALSE)
row_3 <- image_append(c(panel_G, panel_H), stack = FALSE)

row_1 <- image_border(row_1, "white", "20x20")
row_2 <- image_border(row_2, "white", "20x20")
row_3 <- image_border(row_3, "white", "20x20")

figure_width <- single_w * 3 + 40

title_bar <- image_blank(
  width = figure_width,
  height = 250,
  color = "white"
)

title_bar <- image_annotate(
  title_bar,
  text = paste0(
    "Fig. 1 | Programmed cell-death landscape and PA/PQ state identification ",
    "in the 84-sample discovery cohort"
  ),
  gravity = "northwest",
  location = "+55+35",
  size = 58,
  weight = 700,
  color = "black"
)

title_bar <- image_annotate(
  title_bar,
  text = paste0(
    "Stable two-block adjusted intNMF identifies four molecular subtypes; ",
    "C2+C4 define the pyroptosis-activated state (PA) and C1+C3 define the pyroptosis-quiet state (PQ)."
  ),
  gravity = "northwest",
  location = "+58+125",
  size = 34,
  color = "#333333"
)

final_figure <- image_append(
  c(title_bar, row_1, row_2, row_3),
  stack = TRUE
)

final_figure <- image_border(
  final_figure,
  color = "white",
  geometry = "35x35"
)

final_png <- file.path(
  outdir,
  "Fig1_intNMF_K4_PA_PQ_final_main_figure.png"
)

final_pdf <- file.path(
  outdir,
  "Fig1_intNMF_K4_PA_PQ_final_main_figure.pdf"
)

image_write(final_figure, path = final_png, format = "png")
image_write(final_figure, path = final_pdf, format = "pdf")

summary_lines <- c(
  "Final Fig.1 generated from stable two-block adjusted intNMF K4.",
  "",
  "K4 subtype sizes:",
  paste(capture.output(print(table(labels$intNMF_K4))), collapse = " "),
  "",
  "Merged PA/PQ definition:",
  paste0("PA = C2 + C4, n = ", pa_n, "."),
  paste0("PQ = C1 + C3, n = ", pq_n, "."),
  "",
  "Panels:",
  "A: 14-PCD Spearman correlation matrix.",
  "B: K robustness / K selection plus true K4 consensus matrix.",
  "C: MDS of intNMF K4 discovery subtypes.",
  "D: 14 PCD scores across C1-C4.",
  "E: PA vs PQ 14-PCD adjusted effect forest.",
  "F: PA vs PQ confirmatory pyroptosis signatures.",
  "G: Integrated K4/PA-PQ/PCD/pyroptosis heatmap.",
  "H: Core pyroptosis-gene heatmap across PA vs PQ.",
  "",
  paste0("Final PNG: ", final_png),
  paste0("Final PDF: ", final_pdf)
)

writeLines(
  summary_lines,
  file.path(outdir, "Fig1_intNMF_K4_PA_PQ_final_main_figure_summary.txt")
)

cat("\n============================================================\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\nOutput directory:", outdir)
cat("\n============================================================\n")
