#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage:\n")
  cat("  Rscript R/figure2_correlation_piezo_tead_gsdmd.R <RUN_DIR> [SEED_TAG]\n\n")
  cat("Reads:\n")
  cat("  <RUN_DIR>/labels/cluster_results_RNA_<SEED_TAG>.csv\n")
  cat("  <RUN_DIR>/degs_deseq2/vst_matrix.csv\n\n")
  cat("Writes:\n")
  cat("  <RUN_DIR>/plots/Figure2A_PIEZO1_vs_TEAD1_all.png\n")
  cat("  <RUN_DIR>/plots/Figure2B_TEAD1_vs_GSDMD_all.png\n")
  cat("  <RUN_DIR>/plots/Figure2C_clusterwise_corr_heatmap.png\n")
  cat("  <RUN_DIR>/plots/Figure2D_clusterwise_scatter_facets.png\n\n")
  cat("Example:\n")
  cat("  Rscript R/figure2_correlation_piezo_tead_gsdmd.R results/clustering/RNA_only/run_20251225-180901__0c54699caa seed42\n")
  quit(status = 1)
}

run_dir  <- args[1]
seed_tag <- ifelse(length(args) >= 2, args[2], "seed42")

label_file <- file.path(run_dir, "labels", paste0("cluster_results_RNA_", seed_tag, ".csv"))
vst_file   <- file.path(run_dir, "degs_deseq2", "vst_matrix.csv")
out_dir    <- file.path(run_dir, "plots")

if (!file.exists(label_file)) stop("Missing label_file: ", label_file)
if (!file.exists(vst_file))   stop("Missing vst_file: ", vst_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# Load labels (your confirmed schema: Sample_ID, Cluster)
# Also drop spurious row: "Unnamed: 1"
# -----------------------------------------------------------------------------
labels_df <- read_csv(label_file, show_col_types = FALSE) %>%
  rename(sample_id = Sample_ID, cluster = Cluster) %>%
  mutate(
    sample_id = trimws(sample_id),
    cluster   = factor(cluster)
  ) %>%
  filter(sample_id != "Unnamed: 1")

# -----------------------------------------------------------------------------
# Load vst matrix (gene x sample) and convert to sample x gene
# -----------------------------------------------------------------------------
vst_mat <- read_csv(vst_file, show_col_types = FALSE)

if (!("gene" %in% colnames(vst_mat))) {
  # Be tolerant if first column isn't named gene
  colnames(vst_mat)[1] <- "gene"
}

expr_df <- vst_mat %>%
  pivot_longer(cols = -gene, names_to = "sample_id", values_to = "expression") %>%
  pivot_wider(names_from = gene, values_from = expression)

expr_df$sample_id <- trimws(expr_df$sample_id)

# -----------------------------------------------------------------------------
# Merge
# -----------------------------------------------------------------------------
df <- inner_join(labels_df, expr_df, by = "sample_id")

cat("Merged samples:", nrow(df), "\n")

# -----------------------------------------------------------------------------
# Genes (case-sensitive)
# -----------------------------------------------------------------------------
genes_needed <- c("PIEZO1", "TEAD1", "GSDMD")
missing_genes <- setdiff(genes_needed, colnames(df))
if (length(missing_genes) > 0) {
  stop("Missing genes in vst matrix: ", paste(missing_genes, collapse = ", "),
       "\nCheck gene symbols in vst_matrix.csv (case-sensitive).")
}

# =============================================================================
# Helper: correlation annotation text
# =============================================================================
corr_text <- function(x, y) {
  # Spearman
  sp <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  # Pearson
  pe <- suppressWarnings(cor.test(x, y, method = "pearson"))
  paste0(
    "Spearman ρ = ", signif(unname(sp$estimate), 3), ", p = ", format(sp$p.value, digits = 3, scientific = TRUE), "\n",
    "Pearson  r = ", signif(unname(pe$estimate), 3), ", p = ", format(pe$p.value, digits = 3, scientific = TRUE)
  )
}

# =============================================================================
# Figure 2A: PIEZO1 vs TEAD1 (all samples)
# =============================================================================
ann_A <- corr_text(df$PIEZO1, df$TEAD1)

p2a <- ggplot(df, aes(x = PIEZO1, y = TEAD1)) +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate("text",
           x = min(df$PIEZO1, na.rm = TRUE),
           y = max(df$TEAD1, na.rm = TRUE),
           label = ann_A, hjust = 0, vjust = 1, size = 4) +
  labs(
    title = "PIEZO1 vs TEAD1 (all samples)",
    x = "PIEZO1 expression (DESeq2 vst)",
    y = "TEAD1 expression (DESeq2 vst)"
  ) +
  theme_bw(base_size = 14)

ggsave(file.path(out_dir, "Figure2A_PIEZO1_vs_TEAD1_all.png"),
       p2a, width = 6.5, height = 5.2, dpi = 600)

# =============================================================================
# Figure 2B: TEAD1 vs GSDMD (all samples)
# =============================================================================
ann_B <- corr_text(df$TEAD1, df$GSDMD)

p2b <- ggplot(df, aes(x = TEAD1, y = GSDMD)) +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate("text",
           x = min(df$TEAD1, na.rm = TRUE),
           y = max(df$GSDMD, na.rm = TRUE),
           label = ann_B, hjust = 0, vjust = 1, size = 4) +
  labs(
    title = "TEAD1 vs GSDMD (all samples)",
    x = "TEAD1 expression (DESeq2 vst)",
    y = "GSDMD expression (DESeq2 vst)"
  ) +
  theme_bw(base_size = 14)

ggsave(file.path(out_dir, "Figure2B_TEAD1_vs_GSDMD_all.png"),
       p2b, width = 6.5, height = 5.2, dpi = 600)

# =============================================================================
# Figure 2C: Cluster-wise correlations (heatmap of Spearman rho)
# =============================================================================
corr_cluster <- function(d, xcol, ycol) {
  if (nrow(d) < 4) return(data.frame(rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(d[[xcol]], d[[ycol]], method = "spearman", exact = FALSE))
  data.frame(rho = unname(ct$estimate), p = ct$p.value)
}

pairs <- tribble(
  ~pair,                ~x,       ~y,
  "PIEZO1–TEAD1",        "PIEZO1",  "TEAD1",
  "TEAD1–GSDMD",         "TEAD1",   "GSDMD"
)

cluster_stats <- df %>%
  group_by(cluster) %>%
  group_modify(~{
    out <- lapply(seq_len(nrow(pairs)), function(i) {
      st <- corr_cluster(.x, pairs$x[i], pairs$y[i])
      data.frame(pair = pairs$pair[i], rho = st$rho, p = st$p)
    })
    bind_rows(out)
  }) %>%
  ungroup() %>%
  mutate(
    padj = p.adjust(p, method = "fdr"),
    sig  = case_when(
      is.na(padj) ~ "",
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

p2c <- ggplot(cluster_stats, aes(x = pair, y = cluster, fill = rho)) +
  geom_tile() +
  geom_text(aes(label = paste0(signif(rho, 2), sig)), size = 4) +
  labs(
    title = "Cluster-specific Spearman correlations",
    x = NULL,
    y = "Cluster",
    fill = "Spearman ρ"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(file.path(out_dir, "Figure2C_clusterwise_corr_heatmap.png"),
       p2c, width = 7.0, height = 3.8, dpi = 600)

# =============================================================================
# Figure 2D: Faceted scatter plots by cluster (optional)
# =============================================================================
df_long <- df %>%
  select(sample_id, cluster, PIEZO1, TEAD1, GSDMD) %>%
  pivot_longer(cols = c(PIEZO1, TEAD1, GSDMD), names_to = "gene", values_to = "expr")

# Two facet scatter plots combined into one image
p2d1 <- ggplot(df, aes(x = PIEZO1, y = TEAD1)) +
  geom_point(alpha = 0.75, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ cluster, scales = "free") +
  labs(
    title = "PIEZO1 vs TEAD1 (by cluster)",
    x = "PIEZO1 (vst)",
    y = "TEAD1 (vst)"
  ) +
  theme_bw(base_size = 12)

p2d2 <- ggplot(df, aes(x = TEAD1, y = GSDMD)) +
  geom_point(alpha = 0.75, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ cluster, scales = "free") +
  labs(
    title = "TEAD1 vs GSDMD (by cluster)",
    x = "TEAD1 (vst)",
    y = "GSDMD (vst)"
  ) +
  theme_bw(base_size = 12)

p2d <- cowplot::plot_grid(p2d1, p2d2, ncol = 1, labels = c("D", "E"))

ggsave(file.path(out_dir, "Figure2D_clusterwise_scatter_facets.png"),
       p2d, width = 7.5, height = 9.5, dpi = 600)

cat("Saved Figure 2 outputs to:", out_dir, "\n")
