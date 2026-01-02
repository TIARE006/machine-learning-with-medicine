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
  cat("  Rscript R/cluster_phenotype_figure1.R <RUN_DIR> [SEED_TAG]\n\n")
  cat("Reads:\n")
  cat("  <RUN_DIR>/labels/cluster_results_RNA_<SEED_TAG>.csv\n")
  cat("  <RUN_DIR>/degs_deseq2/vst_matrix.csv\n\n")
  cat("Writes:\n")
  cat("  <RUN_DIR>/plots/Figure1_cluster_level_phenotype.png\n\n")
  cat("Example:\n")
  cat("  Rscript R/cluster_phenotype_figure1.R results/clustering/RNA_only/run_20251225-180901__0c54699caa seed42\n")
  quit(status = 1)
}

run_dir  <- args[1]
seed_tag <- ifelse(length(args) >= 2, args[2], "seed42")

label_file <- file.path(run_dir, "labels", paste0("cluster_results_RNA_", seed_tag, ".csv"))
vst_file   <- file.path(run_dir, "degs_deseq2", "vst_matrix.csv")
out_dir    <- file.path(run_dir, "plots")
out_file   <- file.path(out_dir, "Figure1_cluster_level_phenotype.png")

if (!file.exists(label_file)) stop("Missing label_file: ", label_file)
if (!file.exists(vst_file))   stop("Missing vst_file: ", vst_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# Load cluster labels
# -----------------------------------------------------------------------------
labels_df <- readr::read_csv(label_file, show_col_types = FALSE) %>%
  dplyr::rename(sample_id = Sample_ID, cluster = Cluster) %>%
  dplyr::mutate(
    sample_id = trimws(sample_id),
    cluster   = factor(cluster)
  ) %>%
  # Drop the spurious header-like row
  dplyr::filter(sample_id != "Unnamed: 1")

# -----------------------------------------------------------------------------
# Load vst matrix: gene x sample (first col = gene)
# Convert to sample x gene (one row per sample)
# -----------------------------------------------------------------------------
vst_mat <- read_csv(vst_file, show_col_types = FALSE)

if (!("gene" %in% colnames(vst_mat))) {
  stop("Expected first column named 'gene' in vst_matrix.csv, but got: ",
       paste(colnames(vst_mat)[1:3], collapse = ", "))
}

expr_df <- vst_mat %>%
  pivot_longer(cols = -gene, names_to = "sample_id", values_to = "expression") %>%
  pivot_wider(names_from = gene, values_from = expression)
  expr_df$sample_id <- trimws(expr_df$sample_id)

# -----------------------------------------------------------------------------
# Merge
# -----------------------------------------------------------------------------
plot_df <- inner_join(labels_df, expr_df, by = "sample_id")

# -----------------------------------------------------------------------------
# Genes (main text)
# -----------------------------------------------------------------------------
genes_of_interest <- c("PIEZO1", "YAP1", "TEAD1", "GSDMD")

missing_genes <- setdiff(genes_of_interest, colnames(plot_df))
if (length(missing_genes) > 0) {
  stop("Missing genes in vst matrix: ", paste(missing_genes, collapse = ", "),
       "\nCheck gene symbols in vst_matrix.csv (case-sensitive).")
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

figure1 <- cowplot::plot_grid(
  plotlist = plot_list,
  ncol = 2,
  labels = LETTERS[1:length(plot_list)],
  align = "hv"
)

ggsave(out_file, figure1, width = 10, height = 8, dpi = 600)
cat("Saved:", out_file, "\n")
