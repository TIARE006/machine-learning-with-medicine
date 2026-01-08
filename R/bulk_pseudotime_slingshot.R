# ============================================================
# Bulk pseudotime analysis using DiffusionMap + Slingshot
# (NO re-clustering)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(matrixStats)
  library(ggplot2)
  library(slingshot)
  library(destiny)
  library(RColorBrewer)
  library(patchwork)
})

# ============================================================
# 1. Parse arguments
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript bulk_pseudotime_slingshot_DM.R <RUN_DIR>")
}

run_dir <- args[1]
out_dir <- file.path(run_dir, "pseudotime")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

vst_path <- file.path(run_dir, "degs_deseq2", "vst_matrix.csv")
clu_path <- file.path(run_dir, "labels", "cluster_results_RNA_seed42.csv")

stopifnot(file.exists(vst_path), file.exists(clu_path))

# ============================================================
# 2. Read VST matrix
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

# genes x samples → transpose
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
# 5. Select variable genes
# ============================================================
gene_var <- matrixStats::colVars(vst2)
topN <- min(2000, length(gene_var))
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:topN]

X <- scale(vst2[, top_genes, drop = FALSE])

# ============================================================
# 6. Diffusion Map (KEY CHANGE)
# ============================================================
dm <- DiffusionMap(X, n_eigs = 5)

emb <- eigenvectors(dm)[, 1:5]
colnames(emb) <- paste0("DC", 1:5)

pcs <- as.data.frame(emb)
pcs$sample <- rownames(emb)
pcs <- pcs %>% left_join(clu2, by = "sample")

# ============================================================
# 7. Slingshot (MULTI-LINEAGE)
# ============================================================
ROOT_CLUSTER <- "0"
END_CLUSTERS <- c("1", "2")

sce <- slingshot(
  emb,
  clusterLabels = pcs$cluster,
  start.clus = ROOT_CLUSTER,
  end.clus   = END_CLUSTERS
)

pt_mat <- slingPseudotime(sce)

pt <- apply(pt_mat, 1, function(x) {
  if (all(is.na(x))) NA else max(x, na.rm = TRUE)
})

pt <- (pt - min(pt, na.rm = TRUE)) /
      (max(pt, na.rm = TRUE) - min(pt, na.rm = TRUE))

pcs$pseudotime <- pt

# ============================================================
# 8. Composite figure (distribution + branching)
# ============================================================

# Left: distribution
p_left <- ggplot(pcs, aes(DC1, DC2)) +
  geom_point(aes(color = cluster), size = 2.3, alpha = 0.85) +
  theme_bw(base_size = 14) +
  labs(
    title = "Protein-coding genes (bulk)",
    x = "Component 1",
    y = "Component 2",
    color = "Cluster"
  )

# Right: branching trajectory
curve_list <- slingCurves(sce)

curve_df <- purrr::map_dfr(
  seq_along(curve_list),
  function(i) {
    crv <- curve_list[[i]]
    df <- as.data.frame(crv$s[, 1:2])
    colnames(df) <- c("DC1", "DC2")
    df$pseudotime <- crv$lambda
    df <- df[order(df$pseudotime), ]

    tibble(
      DC1 = smooth.spline(df$pseudotime, df$DC1, spar = 0.7)$y,
      DC2 = smooth.spline(df$pseudotime, df$DC2, spar = 0.7)$y,
      lineage = paste0("Lineage_", i)
    )
  }
)

p_right <- ggplot(pcs, aes(DC1, DC2)) +
  geom_point(aes(color = cluster), size = 2.3, alpha = 0.85) +
  geom_path(
    data = curve_df,
    aes(DC1, DC2, group = lineage),
    inherit.aes = FALSE,
    linewidth = 1.6,
    color = "black"
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Bulk pseudotime trajectory with branching",
    x = "Component 1",
    y = "Component 2"
  ) +
  theme(legend.position = "none")

p_final <- p_left + p_right + plot_layout(ncol = 2)

ggsave(
  file.path(out_dir, "Figure_bulk_pseudotime_DM_composite.png"),
  p_final,
  width = 13,
  height = 6,
  dpi = 300
)

# ============================================================
# 9. Pseudotime by cluster
# ============================================================
p2 <- ggplot(pcs, aes(cluster, pseudotime)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
  theme_bw(base_size = 14)

ggsave(file.path(out_dir, "Pseudotime_by_cluster.png"),
       p2, width = 7, height = 5, dpi = 300)

# ============================================================
# 10. Gene expression vs pseudotime
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
  geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, k = 5)) +
  facet_wrap(~gene, scales = "free_y") +
  theme_bw(base_size = 13)

ggsave(file.path(out_dir, "Genes_vs_pseudotime.png"),
       p3, width = 11, height = 7, dpi = 300)

