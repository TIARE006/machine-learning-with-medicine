#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage:\n  Rscript R/Figure2_corr_plus_GDF15.R <RUN_DIR> <SEED_TAG>\n\n")
  cat("Example:\n  Rscript R/Figure2_corr_plus_GDF15.R results/clustering/RNA_only/run_20251225-180901__0c54699caa seed42\n")
  quit(status = 1)
}

run_dir  <- args[1]
seed_tag <- args[2]

labels_csv <- file.path(run_dir, "labels", paste0("cluster_results_RNA_", seed_tag, ".csv"))
vst_csv    <- file.path(run_dir, "degs_deseq2", "vst_matrix.csv")
out_dir    <- file.path(run_dir, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Load labels (Sample_ID, Cluster)
# ----------------------------
labels_df <- read.csv(labels_csv, check.names = FALSE) %>%
  dplyr::rename(sample_id = Sample_ID, cluster = Cluster) %>%
  dplyr::mutate(cluster = factor(cluster))

# ----------------------------
# Load vst matrix (gene x samples) and convert to sample x gene
# Assumption: first column named 'gene' (or similar); remaining columns are sample IDs
# ----------------------------
vst_raw <- read.csv(vst_csv, check.names = FALSE)
gene_col <- colnames(vst_raw)[1]
colnames(vst_raw)[1] <- "gene"

expr_df <- vst_raw %>%
  tidyr::pivot_longer(cols = -gene, names_to = "sample_id", values_to = "expr") %>%
  tidyr::pivot_wider(names_from = gene, values_from = expr)

plot_df <- labels_df %>% inner_join(expr_df, by = "sample_id")

# ----------------------------
# Targets
# ----------------------------
# IMPORTANT: case-sensitive gene symbols as they appear in vst_matrix.csv
genes_needed <- c("PIEZO1", "YAP1", "TEAD1", "GSDMD", "GDF15")

missing_genes <- setdiff(genes_needed, colnames(plot_df))
if (length(missing_genes) > 0) {
  stop(
    "Missing genes in vst_matrix.csv: ",
    paste(missing_genes, collapse = ", "),
    "\nCheck exact gene symbols (case-sensitive)."
  )
}

pairs <- tibble::tribble(
  ~x,        ~y,        ~pair_name,
  "PIEZO1",  "YAP1",    "PIEZO1_vs_YAP1",
  "PIEZO1",  "TEAD1",   "PIEZO1_vs_TEAD1",
  "TEAD1",   "GSDMD",   "TEAD1_vs_GSDMD",
  "YAP1",    "GSDMD",   "YAP1_vs_GSDMD",
  "PIEZO1",  "GSDMD",   "PIEZO1_vs_GSDMD"
)

# ----------------------------
# Helpers
# ----------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-3) return(format(p, scientific = TRUE, digits = 2))
  format(round(p, 4), nsmall = 4)
}

corr_annot_text <- function(df, x, y) {
  xx <- df[[x]]; yy <- df[[y]]
  ok <- is.finite(xx) & is.finite(yy)
  if (sum(ok) < 3) return("Not enough points")

  sp <- suppressWarnings(cor.test(xx[ok], yy[ok], method = "spearman"))
  pe <- suppressWarnings(cor.test(xx[ok], yy[ok], method = "pearson"))

  paste0(
    "Spearman \u03c1 = ", round(unname(sp$estimate), 3), ", p = ", fmt_p(sp$p.value), "\n",
    "Pearson  r = ", round(unname(pe$estimate), 3), ", p = ", fmt_p(pe$p.value)
  )
}

# A) all-samples scatter
plot_scatter_all <- function(df, x, y) {
  ann <- corr_annot_text(df, x, y)

  ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
             label = ann, size = 5) +
    labs(
      title = paste0(x, " vs ", y, " (all samples)"),
      x = paste0(x, " expression (DESeq2 vst)"),
      y = paste0(y, " expression (DESeq2 vst)")
    ) +
    theme_bw(base_size = 18) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# B) by-cluster facets
plot_scatter_by_cluster <- function(df, x, y) {
  ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(size = 2.0, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~ cluster, ncol = 2) +
    labs(
      title = paste0(x, " vs ", y, " (by cluster)"),
      x = paste0(x, " (vst)"),
      y = paste0(y, " (vst)")
    ) +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# C) cluster-wise Spearman rho summary for all pairs
clusterwise_spearman <- function(df, pairs_tbl) {
  df %>%
    tidyr::crossing(pairs_tbl) %>%
    group_by(cluster, pair_name, x, y) %>%
    summarise(
      rho = {
        xx <- .data[[x]]; yy <- .data[[y]]
        ok <- is.finite(xx) & is.finite(yy)
        if (sum(ok) < 3) NA_real_ else suppressWarnings(cor(xx[ok], yy[ok], method = "spearman"))
      },
      .groups = "drop"
    )
}

plot_corr_heatmap <- function(cw_df) {
  # y-axis as cluster, x-axis as pair
  ggplot(cw_df, aes(x = pair_name, y = cluster, fill = rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(rho), "NA", round(rho, 3))), size = 6) +
    labs(title = "Cluster-specific Spearman correlations", x = NULL, y = "Cluster", fill = "Spearman \u03c1") +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
}

# D) GDF15 boxplot + jitter + Kruskal-Wallis
plot_box_jitter_kw <- function(df, gene) {
  kw <- suppressWarnings(kruskal.test(df[[gene]] ~ df$cluster))
  kw_txt <- paste0("Kruskal-Wallis, p = ", fmt_p(kw$p.value))

  ggplot(df, aes(x = cluster, y = .data[[gene]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.75) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
             label = kw_txt, size = 5) +
    labs(
      title = gene,
      x = "Cluster",
      y = "Expression (DESeq2 vst)"
    ) +
    theme_bw(base_size = 18) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# ----------------------------
# Generate outputs
# ----------------------------
# 1) per-pair scatter plots (all + by cluster)
for (i in seq_len(nrow(pairs))) {
  x <- pairs$x[i]; y <- pairs$y[i]; nm <- pairs$pair_name[i]

  p_all <- plot_scatter_all(plot_df, x, y)
  p_by  <- plot_scatter_by_cluster(plot_df, x, y)

  ggsave(file.path(out_dir, paste0("Figure2_", nm, "_all.png")), p_all, width = 8, height = 7, dpi = 400)
  ggsave(file.path(out_dir, paste0("Figure2_", nm, "_by_cluster.png")), p_by, width = 10, height = 8, dpi = 400)
}

# 2) heatmap summary (cluster-wise Spearman for all pairs)
cw <- clusterwise_spearman(plot_df, pairs)
p_heat <- plot_corr_heatmap(cw)
ggsave(file.path(out_dir, "Figure2_clusterwise_corr_heatmap_5pairs.png"), p_heat, width = 12, height = 7, dpi = 400)

# 3) GDF15 boxplot
p_gdf15 <- plot_box_jitter_kw(plot_df, "GDF15")
ggsave(file.path(out_dir, "FigureX_GDF15_cluster_boxplot.png"), p_gdf15, width = 8, height = 7, dpi = 400)

cat("Saved all figures to: ", out_dir, "\n")
