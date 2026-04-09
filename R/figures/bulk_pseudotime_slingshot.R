#!/usr/bin/env Rscript
# ============================================================
# Bulk pseudotime analysis using DiffusionMap + Slingshot
# + Pseudotime heatmap (all genes)
# + Dynamic modules (Early/Mid/Late)
# + TF->Target coupling curves
# + TF enrichment per module (via msigdbr C3:TFT)
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
  library(ComplexHeatmap)
  library(circlize)
  library(msigdbr)
  library(fgsea)
})

# ============================================================
# 1. Parse arguments
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript bulk_pseudotime_slingshot_DM_plusTF.R <RUN_DIR>")
}
run_dir <- args[1]

out_dir <- file.path(run_dir, "pseudotime")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

vst_path <- file.path(run_dir, "degs_deseq2", "vst_matrix.csv")
clu_path <- file.path(run_dir, "labels", "cluster_results_RNA_seed42.csv")
stopifnot(file.exists(vst_path), file.exists(clu_path))

# ============================================================
# 2. Read VST matrix (samples x genes)
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

# if genes x samples, transpose to samples x genes
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
    cluster = as.character(.data[[cluster_col]])
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
# 5. Select variable genes for embedding
# ============================================================
gene_var <- matrixStats::colVars(vst2)
topN <- min(2000, length(gene_var))
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:topN]

X <- scale(vst2[, top_genes, drop = FALSE])

# ============================================================
# 6. Diffusion Map
# ============================================================
dm <- DiffusionMap(X, n_eigs = 5)
emb <- eigenvectors(dm)[, 1:5]
colnames(emb) <- paste0("DC", 1:5)

pcs <- as.data.frame(emb)
pcs$sample <- rownames(emb)
pcs <- pcs %>% left_join(clu2, by = "sample")

# ============================================================
# 7. Slingshot (multi-lineage)
# ============================================================
ROOT_CLUSTER <- "0"
END_CLUSTERS <- c("1", "2")

sce <- slingshot(
  emb,
  clusterLabels = factor(pcs$cluster),
  start.clus = ROOT_CLUSTER,
  end.clus   = END_CLUSTERS
)

pt_mat <- slingPseudotime(sce)

# Combine lineages (bulk): take max pseudotime across lineages
pt <- apply(pt_mat, 1, function(x) {
  if (all(is.na(x))) NA else max(x, na.rm = TRUE)
})

# Normalize to [0,1]
pt <- (pt - min(pt, na.rm = TRUE)) /
      (max(pt, na.rm = TRUE) - min(pt, na.rm = TRUE))

pcs$pseudotime <- pt

# Save pseudotime table
write.csv(pcs %>% select(sample, cluster, pseudotime),
          file.path(out_dir, "Pseudotime_table.csv"),
          row.names = FALSE)

# ============================================================
# 8. Composite figure (distribution + branching)
# ============================================================
p_left <- ggplot(pcs, aes(DC1, DC2)) +
  geom_point(aes(color = cluster), size = 2.3, alpha = 0.85) +
  theme_bw(base_size = 14) +
  labs(
    title = "Bulk diffusion map embedding",
    x = "DC1",
    y = "DC2",
    color = "Cluster"
  )

curve_list <- slingCurves(sce)
curve_df <- purrr::map_dfr(
  seq_along(curve_list),
  function(i) {
    crv <- curve_list[[i]]
    df <- as.data.frame(crv$s[, 1:2])
    colnames(df) <- c("DC1", "DC2")
    df$pseudotime <- crv$lambda
    df <- df[order(df$pseudotime), , drop = FALSE]

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
    x = "DC1",
    y = "DC2"
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
  theme_bw(base_size = 14) +
  labs(title = "Pseudotime distribution by cluster")
ggsave(file.path(out_dir, "Pseudotime_by_cluster.png"),
       p2, width = 7, height = 5, dpi = 300)

# ============================================================
# 10. Gene expression vs "bulk pseudotime" (FIXED)
#     OLD: sample-level scatter + GAM smooth
#     NEW: cluster-level mean ± SEM along enforced order 0 -> 3 -> 1 -> 2
# ============================================================
genes_to_plot <- c("GAPDH","TRIM72","LIF","PIEZO1","PIEZO2","GSDMD","MMP14","HMGCS2","YAP1","TEAD1","GDF15")

# keep the old sample-level plot as a reference
expr_long <- vst2 %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  select(any_of(c("sample", genes_to_plot))) %>%
  pivot_longer(-sample, names_to = "gene", values_to = "expr") %>%
  left_join(pcs %>% select(sample, cluster, pseudotime), by = "sample")

p3_old <- ggplot(expr_long, aes(pseudotime, expr, color = cluster)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, k = 5)) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 13) +
  labs(title = "Selected genes vs pseudotime (sample-level; reference)", x = "Pseudotime", y = "VST")
ggsave(file.path(out_dir, "Genes_vs_pseudotime_samplelevel.png"),
       p3_old, width = 13, height = 8, dpi = 300)

# ---- NEW: cluster-level trajectory plot ----
# Enforced biological order (your requirement)
cluster_order <- c("0","3","1","2")

# Build a discrete pseudo-stage axis from that order
stage_map <- tibble(
  cluster = cluster_order,
  stage = seq_along(cluster_order) - 1
)

# summarize mean ± SEM per (gene, cluster), then map to stage
expr_sum <- expr_long %>%
  filter(!is.na(cluster), cluster %in% cluster_order) %>%
  group_by(gene, cluster) %>%
  summarise(
    mean_expr = mean(expr, na.rm = TRUE),
    sd_expr   = sd(expr, na.rm = TRUE),
    n         = sum(!is.na(expr)),
    sem_expr  = sd_expr / sqrt(pmax(n, 1)),
    .groups = "drop"
  ) %>%
  left_join(stage_map, by = "cluster") %>%
  mutate(
    cluster = factor(cluster, levels = cluster_order)
  ) %>%
  arrange(gene, stage)

# plotting: x = stage (0,1,2,3), label with cluster order
p3_new <- ggplot(expr_sum, aes(x = stage, y = mean_expr, group = 1)) +
  geom_ribbon(aes(ymin = mean_expr - sem_expr, ymax = mean_expr + sem_expr),
              alpha = 0.25, color = NA) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  scale_x_continuous(
    breaks = stage_map$stage,
    labels = paste0("C", stage_map$cluster),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Selected genes along enforced cluster progression (C0 → C3 → C1 → C2)",
    x = "Ordered clusters (discrete pseudo-time)",
    y = "Mean VST (± SEM)"
  )

ggsave(file.path(out_dir, "Genes_vs_pseudotime.png"),
       p3_new, width = 13, height = 8, dpi = 300)

# ============================================================
# 11. Pseudotime heatmap for many genes (paper-like panel c)
#     - choose top variable genes
#     - order samples by pseudotime
#     - z-score per gene
# ============================================================
expr_gs <- t(vst2)   # genes x samples
sample_order <- order(pcs$pseudotime)
samples_ord <- pcs$sample[sample_order]
pt_ord <- pcs$pseudotime[sample_order]
clu_ord <- pcs$cluster[sample_order]

gvar <- matrixStats::rowVars(expr_gs, na.rm = TRUE)
topK <- min(1500, length(gvar))
genes_hm <- names(sort(gvar, decreasing = TRUE))[1:topK]

hm_mat <- expr_gs[genes_hm, samples_ord, drop = FALSE]
hm_mat_z <- t(scale(t(hm_mat)))   # gene-wise z-score
hm_mat_z[is.na(hm_mat_z)] <- 0

peak_idx <- apply(hm_mat_z, 1, which.max)

q <- quantile(peak_idx, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
module <- cut(peak_idx, breaks = q, include.lowest = TRUE, labels = c("Early","Mid","Late"))

row_ord <- order(module, peak_idx)

ha_col <- HeatmapAnnotation(
  Cluster = clu_ord,
  Pseudotime = anno_simple(pt_ord),
  col = list(
    Cluster = structure(
      brewer.pal(max(3, length(unique(pcs$cluster))), "Set2")[seq_along(sort(unique(pcs$cluster)))],
      names = sort(unique(pcs$cluster))
    )
  ),
  annotation_name_side = "left"
)

ha_row <- rowAnnotation(
  Module = module[row_ord],
  col = list(Module = c(Early="#D55E00", Mid="#E69F00", Late="#0072B2"))
)

ht <- Heatmap(
  hm_mat_z[row_ord, , drop = FALSE],
  name = "Z",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha_col,
  left_annotation = ha_row,
  column_title = "Genes ordered by dynamic peak along pseudotime (bulk)",
  use_raster = TRUE
)

png(file.path(out_dir, "Heatmap_pseudotime_allgenes.png"), width = 1400, height = 900, res = 150)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

module_tbl <- tibble(gene = rownames(hm_mat_z), module = module) %>%
  arrange(module)
write.csv(module_tbl, file.path(out_dir, "Gene_modules_EarlyMidLate.csv"), row.names = FALSE)

# ============================================================
# 12. TF enrichment per module (paper-like panel b/f)
# ============================================================
message("Running TF enrichment (MSigDB C3:TFT-like) ...")

m_all_c3 <- msigdbr(
  species = "Homo sapiens",
  collection = "C3"
)

m_tft <- m_all_c3 %>%
  filter(str_detect(gs_name, "^TFT_")) %>%
  select(gs_name, gene_symbol) %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

tft_sets <- split(m_tft$gene_symbol, m_tft$gs_name)

universe <- unique(m_tft$gene_symbol)
universe <- intersect(universe, rownames(expr_gs))

if (length(universe) < 50 || length(tft_sets) < 10) {
  warning("TF enrichment skipped: TFT universe/gene sets too small after intersection.")
  write.csv(
    tibble(module=character(), gene_set=character(), overlap=integer(),
           odds_ratio=double(), pvalue=double(), padj=double()),
    file.path(out_dir, "TF_enrichment_ORA_C3TFT.csv"),
    row.names = FALSE
  )
} else {

  ora_one <- function(gset, module_genes, universe) {
    out <- tibble(
      overlap = NA_integer_,
      odds_ratio = NA_real_,
      pvalue = NA_real_
    )

    gset <- intersect(gset, universe)
    module_genes <- intersect(module_genes, universe)

    if (length(gset) < 10 || length(module_genes) < 10) return(out)

    a <- length(intersect(gset, module_genes))
    b <- length(setdiff(module_genes, gset))
    c <- length(setdiff(gset, module_genes))
    d <- length(setdiff(universe, union(gset, module_genes)))

    if (any(c(a, b, c, d) < 0) || (a + b + c + d) == 0) return(out)

    ft <- tryCatch(
      fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater"),
      error = function(e) NULL
    )
    if (is.null(ft)) return(out)

    out$overlap    <- a
    out$odds_ratio <- unname(ft$estimate)
    out$pvalue     <- ft$p.value
    out
  }

  module_genes_list <- split(module_tbl$gene, module_tbl$module)

  ora_res <- purrr::imap_dfr(
    module_genes_list,
    function(glist, mlabel) {
      purrr::imap_dfr(
        tft_sets,
        function(gset, nm) {
          tmp <- ora_one(gset, glist, universe)
          tmp$gene_set <- nm
          tmp$module   <- as.character(mlabel)
          tmp
        }
      )
    }
  ) %>%
    mutate(
      pvalue = as.numeric(pvalue),
      padj   = p.adjust(pvalue, method = "BH")
    ) %>%
    arrange(module, padj, desc(odds_ratio))

  write.csv(ora_res, file.path(out_dir, "TF_enrichment_ORA_C3TFT.csv"), row.names = FALSE)

  top_plot <- ora_res %>%
    filter(!is.na(padj)) %>%
    group_by(module) %>%
    slice_min(order_by = padj, n = 12, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(label = str_replace(gene_set, "^TFT_", "")) %>%
    mutate(label = str_replace(label, "_TARGETS.*$", ""))

  if (nrow(top_plot) == 0) {
    warning("No TF sets available for plotting after ORA; skipping TF_enrichment_by_module.png")
  } else {
    p_tf <- ggplot(top_plot, aes(x = reorder(label, -log10(padj)), y = -log10(padj))) +
      geom_col() +
      coord_flip() +
      facet_wrap(~module, scales = "free_y") +
      theme_bw(base_size = 12) +
      labs(
        title = "TF enrichment per dynamic module (MSigDB C3 TFT_* ORA)",
        x = "TF (from gene set name)",
        y = "-log10(FDR)"
      )

    ggsave(file.path(out_dir, "TF_enrichment_by_module.png"),
           p_tf, width = 12, height = 7, dpi = 300)
  }
}

# ============================================================
# 13. TF -> Target coupling curves (paper-like panel e)
# ============================================================
tf_targets <- list(
  STAT3 = c("LIF","SOCS3","CEBPD","CXCL8","IL6ST"),
  NFkB  = c("NFKBIA","TNFAIP3","ICAM1","CXCL2","CCL2"),
  TEAD1 = c("CTGF","CYR61","ANKRD1","COL1A1","COL1A2"),
  YAP1  = c("CTGF","CYR61","ANKRD1","AMOTL2","BIRC5"),
  HIF1A = c("SLC2A1","ENO1","LDHA","PGK1","HK2"),
  PPARG = c("CD36","CPT1A","ACOX1","FABP4","LPL")
)

gene_z <- function(g){
  if (!g %in% rownames(expr_gs)) return(rep(NA_real_, length(samples_ord)))
  as.numeric(scale(expr_gs[g, samples_ord]))
}

target_score <- function(targets){
  tg <- intersect(targets, rownames(expr_gs))
  if (length(tg) < 3) {
    return(rep(NA_real_, length(samples_ord)))
  }
  mat <- t(apply(
    expr_gs[tg, samples_ord, drop = FALSE],
    1,
    function(v) as.numeric(scale(v)))
  )
  colMeans(mat, na.rm = TRUE)
}

couple_df <- purrr::imap_dfr(tf_targets, function(tg, tf){
  tibble(
    sample = samples_ord,
    pseudotime = pt_ord,
    cluster = clu_ord,
    TF = tf,
    tf_z = gene_z(tf),
    target_z = target_score(tg)
  )
})

couple_df <- couple_df %>% filter(!is.na(tf_z), !is.na(target_z))

p_couple <- ggplot(couple_df, aes(pseudotime)) +
  geom_smooth(aes(y = tf_z), se = FALSE, method = "gam", formula = y ~ s(x, k = 5)) +
  geom_smooth(aes(y = target_z), se = FALSE, method = "gam", formula = y ~ s(x, k = 5),
              linetype = "dashed") +
  facet_wrap(~TF, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 12) +
  labs(
    title = "TF (solid) vs Target-module score (dashed) along pseudotime",
    x = "Pseudotime", y = "Z-score"
  )

ggsave(file.path(out_dir, "TF_Target_coupling.png"),
       p_couple, width = 12, height = 8, dpi = 300)

# ============================================================
# 14. Save workspace
# ============================================================
save.image(file = file.path(out_dir, "all_results_DM_plusTF.RData"))

message("Done. Outputs written to: ", out_dir)

