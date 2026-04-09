#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage:\n")
  cat("  Rscript R/figures/cluster_phenotype_figure1_plus.R <RUN_DIR> [SEED_TAG]\n\n")
  cat("Reads:\n")
  cat("  <RUN_DIR>/labels/cluster_results_RNA_<SEED_TAG>.csv\n")
  cat("  <RUN_DIR>/degs_deseq2/vst_matrix.csv\n\n")
  cat("Writes:\n")
  cat("  <RUN_DIR>/plots/Figure1plus_cluster_level_phenotype.png\n")
  cat("  <RUN_DIR>/plots/Figure1plus_cluster_level_phenotype.pdf\n\n")
  cat("Example:\n")
  cat("  Rscript R/figures/cluster_phenotype_figure1_plus.R results/clustering/RNA_only/run_20251225-180901__0c54699caa seed42\n")
  quit(status = 1)
}

run_dir  <- args[1]
seed_tag <- ifelse(length(args) >= 2, args[2], "seed42")

label_file <- file.path(run_dir, "labels", paste0("cluster_results_RNA_", seed_tag, ".csv"))
vst_file   <- file.path(run_dir, "degs_deseq2", "vst_matrix.csv")
out_dir    <- file.path(run_dir, "plots")
out_png    <- file.path(out_dir, "Figure1plus_cluster_level_phenotype.png")
out_pdf    <- file.path(out_dir, "Figure1plus_cluster_level_phenotype.pdf")

if (!file.exists(label_file)) stop("Missing label_file: ", label_file)
if (!file.exists(vst_file))   stop("Missing vst_file: ", vst_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# Load labels (your confirmed schema: Sample_ID, Cluster)
# -----------------------------------------------------------------------------
labels_df <- read_csv(label_file, show_col_types = FALSE) %>%
  rename(sample_id = Sample_ID, cluster = Cluster) %>%
  mutate(
    sample_id = trimws(sample_id),
    cluster   = factor(cluster)
  ) %>%
  filter(sample_id != "Unnamed: 1")

# -----------------------------------------------------------------------------
# Load vst matrix (gene x sample) -> sample x gene
# -----------------------------------------------------------------------------
vst_mat <- read_csv(vst_file, show_col_types = FALSE)
if (!("gene" %in% colnames(vst_mat))) colnames(vst_mat)[1] <- "gene"

expr_df <- vst_mat %>%
  pivot_longer(cols = -gene, names_to = "sample_id", values_to = "expression") %>%
  pivot_wider(names_from = gene, values_from = expression)

expr_df$sample_id <- trimws(expr_df$sample_id)

# -----------------------------------------------------------------------------
# Merge
# -----------------------------------------------------------------------------
plot_df <- inner_join(labels_df, expr_df, by = "sample_id")
cat("Merged samples:", nrow(plot_df), "\n")

# -----------------------------------------------------------------------------
# Genes for Figure 1+
# (Original 4) + PIEZO2 + TRIM72(MG53)
# -----------------------------------------------------------------------------
genes_of_interest <- c("PIEZO1", "PIEZO2", "YAP1", "TEAD1", "GSDMD", "TRIM72")

missing_genes <- setdiff(genes_of_interest, colnames(plot_df))
if (length(missing_genes) > 0) {
  stop(
    "Missing genes in vst matrix: ", paste(missing_genes, collapse = ", "),
    "\nCheck gene symbols in vst_matrix.csv (case-sensitive)."
  )
}

# -----------------------------------------------------------------------------
# One-panel plot: boxplot + jitter + Kruskal–Wallis p
# -----------------------------------------------------------------------------
plot_gene <- function(df, gene) {
  ggplot(df, aes(x = cluster, y = .data[[gene]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.7) +
    stat_compare_means(method = "kruskal.test", label.y.npc = 0.95, size = 4) +
    labs(
      title = gene,
      x = "Cluster",
      y = "Expression (DESeq2 vst)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 11)
    )
}

plot_list <- lapply(genes_of_interest, function(g) plot_gene(plot_df, g))

# Layout: 2 columns x 3 rows (6 panels)
fig <- cowplot::plot_grid(
  plotlist = plot_list,
  ncol = 2,
  labels = LETTERS[1:length(plot_list)],
  align = "hv"
)

# Save both PNG and PDF (PDF preferred for manuscripts)
ggsave(out_png, fig, width = 10, height = 12, dpi = 600)
ggsave(out_pdf, fig, width = 10, height = 12)
cat("Saved:\n  ", out_png, "\n  ", out_pdf, "\n")
