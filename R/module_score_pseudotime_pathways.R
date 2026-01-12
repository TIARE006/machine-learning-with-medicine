#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript R/module_score_pseudotime_pathways.R <RUN_DIR>")
}
run_dir <- args[1]

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
vst_path <- file.path(run_dir, "degs_deseq2", "vst_matrix.csv")
lab_dir  <- file.path(run_dir, "labels")
pt_dir   <- file.path(run_dir, "pseudotime")
out_dir  <- file.path(pt_dir, "D_pathways")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

label_files <- list.files(lab_dir, pattern="^cluster_results_.*\\.csv$", full.names=TRUE)
if (length(label_files) < 1) stop("No cluster_results_*.csv found")
clu_path <- label_files[1]

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
vst <- fread(vst_path, data.table = FALSE, check.names = FALSE)
if (!is.numeric(vst[[1]][1])) {
  rownames(vst) <- vst[[1]]
  vst[[1]] <- NULL
}
expr <- as.matrix(vst)
storage.mode(expr) <- "numeric"  # genes x samples

clu <- fread(clu_path, data.table = FALSE)
sample_col  <- names(clu)[grepl("sample|id", tolower(names(clu)))] [1]
cluster_col <- names(clu)[grepl("cluster", tolower(names(clu)))] [1]
clu <- clu %>%
  transmute(sample = as.character(.data[[sample_col]]),
            cluster = as.character(.data[[cluster_col]]))

common <- intersect(colnames(expr), clu$sample)
expr <- expr[, common, drop=FALSE]
clu  <- clu %>% filter(sample %in% common) %>%
  arrange(match(sample, colnames(expr)))

ORDER <- c("0","3","1","2")
clu$cluster <- factor(clu$cluster, levels=ORDER)

# ------------------------------------------------------------
# Define pathways (Hallmark-style, curated)
# ------------------------------------------------------------
PATHWAYS <- list(
  TNFA_NFKB = c("NFKBIA","RELA","TNFAIP3","ICAM1","CXCL8"),
  JAK_STAT  = c("LIF","STAT3","SOCS3","IL6ST"),
  TGF_BETA = c("TGFB1","SMAD3","SMAD2","CTGF","COL1A1"),
  ECM_EMT  = c("MMP14","COL1A1","COL1A2","CTGF"),
  FA_METAB = c("HMGCS2","CPT1A","CD36","ACOX1"),
  GLYCOLYSIS = c("GAPDH","ENO1","PKM","LDHA")
)

# ------------------------------------------------------------
# Compute module scores
# ------------------------------------------------------------
z_expr <- t(scale(t(expr)))  # gene-wise z-score

module_scores <- lapply(PATHWAYS, function(genes){
  g <- intersect(genes, rownames(z_expr))
  if (length(g) < 3) return(rep(NA, ncol(z_expr)))
  colMeans(z_expr[g, , drop=FALSE], na.rm=TRUE)
})

module_scores <- as.data.frame(module_scores)
module_scores$sample <- colnames(expr)

df <- module_scores %>%
  left_join(clu, by="sample") %>%
  mutate(pseudotime = as.numeric(cluster))

write.csv(df, file.path(out_dir, "module_scores_all_samples.csv"), row.names=FALSE)

# ------------------------------------------------------------
# Plot: pathway trend along ordered clusters
# ------------------------------------------------------------
trend <- df %>%
  pivot_longer(cols = names(PATHWAYS), names_to="pathway", values_to="score") %>%
  group_by(pathway, cluster) %>%
  summarise(mean = mean(score, na.rm=TRUE),
            sem  = sd(score, na.rm=TRUE)/sqrt(sum(!is.na(score))),
            .groups="drop")

p <- ggplot(trend, aes(cluster, mean, group=pathway)) +
  geom_line() +
  geom_point(size=2) +
  geom_ribbon(aes(ymin=mean-sem, ymax=mean+sem),
              fill="grey70", alpha=0.4) +
  facet_wrap(~pathway, scales="free_y", ncol=3) +
  theme_bw(base_size=12) +
  labs(
    title="Pathway module scores along ordered clusters (bulk)",
    x="Cluster order (C0→C3→C1→C2)",
    y="Module score (z-score mean)"
  )

ggsave(file.path(out_dir, "Module_scores_trend.png"),
       p, width=13, height=7, dpi=300)

# ------------------------------------------------------------
# Spearman monotonicity (pathway-level)
# ------------------------------------------------------------
time_vec <- 0:3
stats <- trend %>%
  group_by(pathway) %>%
  summarise(
    spearman_rho = suppressWarnings(cor(mean, time_vec, method="spearman")),
    pvalue = suppressWarnings(cor.test(mean, time_vec, method="spearman")$p.value),
    .groups="drop"
  ) %>%
  mutate(padj = p.adjust(pvalue, method="BH")) %>%
  arrange(padj)

write.csv(stats, file.path(out_dir, "Module_scores_monotonicity.csv"),
          row.names=FALSE)

message("Done. Results written to: ", out_dir)

