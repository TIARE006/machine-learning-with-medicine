# ============================================================
# Bulk pseudotime analysis using Slingshot (NO re-clustering)
# Input:
#   - vst_matrix.csv
#   - cluster_results_RNA_seed42.csv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(matrixStats)
  library(ggplot2)
  library(slingshot)
  library(RColorBrewer)
})

# ============================================================
# 1. Parse arguments
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript bulk_pseudotime_slingshot.R <RUN_DIR>")
}

run_dir <- args[1]
out_dir <- file.path(run_dir, "pseudotime")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

vst_path <- file.path(run_dir, "degs_deseq2", "vst_matrix.csv")
clu_path <- file.path(run_dir, "labels", "cluster_results_RNA_seed42.csv")

cat("Using vst:", vst_path, "\n")
cat("Using labels:", clu_path, "\n")
stopifnot(file.exists(vst_path), file.exists(clu_path))

# ============================================================
# 2. Read VST matrix robustly
# ============================================================
read_vst <- function(path){
  df <- read.csv(path, check.names = FALSE)
  if (!is.numeric(df[[1]][1])) {
    rownames(df) <- df[[1]]
    df[[1]] <- NULL
  }
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  mat
}

vst <- read_vst(vst_path)

# If matrix is genes x samples → transpose
if (nrow(vst) > 1000 && ncol(vst) < 500) {
  vst <- t(vst)
}

# ============================================================
# 3. Read cluster labels
# ============================================================
clu <- read.csv(clu_path, check.names = FALSE)

sample_col  <- names(clu)[str_detect(tolower(names(clu)), "sample|id")][1]
cluster_col <- names(clu)[str_detect(tolower(names(clu)), "cluster")][1]

clu <- clu %>%
  transmute(
    sample  = as.character(.data[[sample_col]]),
    cluster = factor(.data[[cluster_col]])
  )

# ============================================================
# 4. Align samples
# ============================================================
common <- intersect(rownames(vst), clu$sample)
stopifnot(length(common) >= 10)

vst2 <- vst[common, , drop = FALSE]
clu2 <- clu %>%
  filter(sample %in% common) %>%
  arrange(match(sample, common))

# ============================================================
# 5. Select variable genes (bulk pseudotime stability)
# ============================================================
gene_var <- matrixStats::colVars(vst2)
topN <- min(2000, length(gene_var))
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:topN]

X <- scale(vst2[, top_genes, drop = FALSE])

# ============================================================
# 6. PCA
# ============================================================
pca <- prcomp(X, center = FALSE, scale. = FALSE)

pcs <- as.data.frame(pca$x[, 1:5])
pcs$sample <- rownames(pcs)
pcs <- pcs %>% left_join(clu2, by = "sample")

# ============================================================
# 7. Slingshot (NO reclustering)
# ============================================================
ROOT_CLUSTER <- "0"   # early state (can change to 1/2/3)

sce <- slingshot(
  pca$x[, 1:5],
  clusterLabels = pcs$cluster,
  start.clus = ROOT_CLUSTER
)

pt <- slingPseudotime(sce)[, 1]
pt <- (pt - min(pt, na.rm = TRUE)) /
      (max(pt, na.rm = TRUE) - min(pt, na.rm = TRUE))

pcs$pseudotime <- pt

write.csv(
  pcs %>% select(sample, cluster, pseudotime),
  file = file.path(out_dir, "bulk_pseudotime_table.csv"),
  row.names = FALSE
)

# ============================================================
# 8. Plot 1 (CLEAN): PCA + single principal trajectory
# ============================================================

# Extract first lineage principal curve
curve_df <- as.data.frame(slingCurves(sce)[[1]]$s[, 1:2])
colnames(curve_df) <- c("PC1", "PC2")

p1 <- ggplot(pcs, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_path(
    data = curve_df,
    aes(PC1, PC2),
    inherit.aes = FALSE,
    linewidth = 1.2,
    color = "black"
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Bulk pseudotime trajectory (Slingshot)",
    color = "Cluster"
  )

ggsave(
  file.path(out_dir, "Pseudotime_PCA_slingshot.png"),
  p1, width = 8, height = 6, dpi = 300
)

# ============================================================
# 9. Plot 2: pseudotime by cluster
# ============================================================
p2 <- ggplot(pcs, aes(cluster, pseudotime)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
  theme_bw(base_size = 14) +
  labs(title = "Pseudotime distribution by cluster",
       y = "Pseudotime (0–1)")

ggsave(file.path(out_dir, "Pseudotime_by_cluster.png"),
       p2, width = 7, height = 5, dpi = 300)

# ============================================================
# 10. Plot 3: Gene expression vs pseudotime
# ============================================================
genes_to_plot <- c("PIEZO1", "YAP1", "TEAD1", "GSDMD", "GDF15")

expr_long <- vst2 %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  select(any_of(c("sample", genes_to_plot))) %>%
  pivot_longer(-sample, names_to = "gene", values_to = "expr") %>%
  left_join(pcs %>% select(sample, cluster, pseudotime), by = "sample")

p3 <- ggplot(expr_long, aes(pseudotime, expr, color = cluster)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_smooth(se = FALSE, method = "gam",
              formula = y ~ s(x, k = 5)) +
  facet_wrap(~gene, scales = "free_y") +
  theme_bw(base_size = 13) +
  labs(title = "Gene expression dynamics along bulk pseudotime",
       x = "Pseudotime (0–1)",
       y = "Expression (VST)")

ggsave(file.path(out_dir, "Genes_vs_pseudotime.png"),
       p3, width = 11, height = 7, dpi = 300)

