#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage:\n")
  cat("  Rscript R/rna_marker_heatmap.R <RUN_DIR> [TOPN] [FDR] [LFC] [LABELS_CSV]\n\n")
  cat("Inputs expected:\n")
  cat("  <RUN_DIR>/degs_deseq2/vst_matrix.csv\n")
  cat("  <RUN_DIR>/degs_deseq2/DEG_cluster_*_vs_rest.csv\n\n")
  cat("Labels (optional):\n")
  cat("  default: <RUN_DIR>/labels/cluster_labels.csv\n")
  cat("  CSV columns supported:\n")
  cat("    sample/cluster  OR  Sample_ID/cluster\n")
  quit(status = 1)
}

run_dir  <- args[1]
topn     <- ifelse(length(args) >= 2, as.integer(args[2]), 15)
fdr_cut  <- ifelse(length(args) >= 3, as.numeric(args[3]), 0.05)
lfc_cut  <- ifelse(length(args) >= 4, as.numeric(args[4]), 1.0)

# labels file: user-specified > default
labels_fp <- if (length(args) >= 5) {
  args[5]
} else {
  file.path(run_dir, "labels", "cluster_labels.csv")
}

degs_dir <- file.path(run_dir, "degs_deseq2")
vst_fp   <- file.path(degs_dir, "vst_matrix.csv")
out_dir  <- file.path(run_dir, "heatmaps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(vst_fp)) stop("Missing vst_matrix.csv: ", vst_fp)

deg_files <- list.files(degs_dir, pattern="^DEG_cluster_\\d+_vs_rest\\.csv$", full.names=TRUE)
if (length(deg_files) == 0) stop("No DEG files found in: ", degs_dir)

message(">>> [Heatmap] run_dir: ", run_dir)
message(">>> [Heatmap] vst    : ", vst_fp)
message(">>> [Heatmap] degs   : ", degs_dir)
message(">>> [Heatmap] out    : ", out_dir)
message(">>> [Heatmap] params : TOPN=", topn, " FDR<=", fdr_cut, " |log2FC|>=", lfc_cut)
message(">>> [Heatmap] labels : ", labels_fp)

# ---------- helper: robustly read vst matrix ----------
# Expect: first column = gene (or named gene/symbol); columns = samples
vst_raw <- suppressMessages(read_csv(vst_fp, show_col_types = FALSE))

# detect gene column name
gene_col <- NULL
cand <- c("gene", "Gene", "symbol", "SYMBOL", "feature", "Feature", "id", "ID")
for (cc in cand) {
  if (cc %in% colnames(vst_raw)) { gene_col <- cc; break }
}
if (is.null(gene_col)) gene_col <- colnames(vst_raw)[1]

vst_mat <- vst_raw %>%
  rename(gene = all_of(gene_col)) %>%
  mutate(gene = as.character(gene)) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene, .keep_all = TRUE)

vst_m <- vst_mat %>%
  select(-gene) %>%
  as.data.frame()
rownames(vst_m) <- vst_mat$gene

# ---------- collect top markers per cluster ----------
markers_list <- list()

for (fp in deg_files) {
  bn <- basename(fp)
  k <- str_match(bn, "DEG_cluster_(\\d+)_vs_rest\\.csv")[,2]
  message(">>> [Heatmap] reading DEG: cluster ", k, " | ", bn)

  df <- suppressMessages(read_csv(fp, show_col_types = FALSE))

  # gene column
  if (!("gene" %in% colnames(df))) {
    g2 <- intersect(colnames(df), c("Gene", "SYMBOL", "symbol", "gene_name", "genes"))
    if (length(g2) == 0) stop("No gene column in: ", fp)
    df <- df %>% rename(gene = all_of(g2[1]))
  }

  # FDR / padj column
  if (!("FDR" %in% colnames(df))) {
    f2 <- intersect(colnames(df), c("padj", "p_adj", "adj_p", "qvalue", "Fdr", "fdr"))
    if (length(f2) == 0) stop("No FDR/padj column in: ", fp)
    df <- df %>% rename(FDR = all_of(f2[1]))
  }

  # log2FC column
  if (!("log2FC" %in% colnames(df))) {
    l2 <- intersect(colnames(df), c("log2FoldChange", "logFC", "log2_fc", "lfc", "LFC"))
    if (length(l2) == 0) stop("No log2FC/log2FoldChange column in: ", fp)
    df <- df %>% rename(log2FC = all_of(l2[1]))
  }

  df_sig <- df %>%
    mutate(gene = as.character(gene)) %>%
    filter(!is.na(FDR), !is.na(log2FC)) %>%
    filter(FDR <= fdr_cut, abs(log2FC) >= lfc_cut)

  message("    cluster ", k, ": significant genes = ", nrow(df_sig))
  if (nrow(df_sig) == 0) next

  top_genes <- df_sig %>%
    arrange(desc(abs(log2FC)), FDR) %>%
    slice_head(n = topn) %>%
    pull(gene) %>%
    unique()

  markers_list[[paste0("cluster_", k)]] <- top_genes
}

markers <- unique(unlist(markers_list))
if (length(markers) == 0) stop("No markers found under current thresholds. Try relax FDR/LFC or reduce TOPN.")

# save marker list for reporting
marker_out <- bind_rows(lapply(names(markers_list), function(nm) {
  tibble(cluster = nm, gene = markers_list[[nm]])
}))
marker_csv <- file.path(out_dir, paste0("top_markers_TOP", topn, "_FDR", fdr_cut, "_LFC", lfc_cut, ".csv"))
write_csv(marker_out, marker_csv)
message(">>> [Heatmap] saved markers: ", marker_csv)

# subset VST matrix by markers
markers_in_vst <- intersect(markers, rownames(vst_m))
if (length(markers_in_vst) < 5) {
  stop("Too few markers found in vst_matrix (matched=", length(markers_in_vst),
       "). Check gene naming consistency (SYMBOL vs Ensembl).")
}
vst_sub <- vst_m[markers_in_vst, , drop = FALSE]

# Row z-score
zscore_rows <- function(mat) {
  m <- as.matrix(mat)
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  (m - mu) / sdv
}
zmat <- zscore_rows(vst_sub)

# ---------- sample->cluster labels (optional but recommended) ----------
anno_col <- NULL
gaps_col <- NULL

if (!is.na(labels_fp) && file.exists(labels_fp)) {
  lab <- suppressMessages(read_csv(labels_fp, show_col_types = FALSE))

  # support Sample_ID/cluster OR sample/cluster
  if ("Sample_ID" %in% colnames(lab)) {
    lab <- lab %>% rename(sample = Sample_ID)
  }
  if (!("sample" %in% colnames(lab))) {
    s2 <- intersect(colnames(lab), c("sample", "Sample", "samples", "SampleID", "sample_id"))
    if (length(s2) > 0) lab <- lab %>% rename(sample = all_of(s2[1]))
  }
  if (!("cluster" %in% colnames(lab))) {
    c2 <- intersect(colnames(lab), c("cluster", "Cluster", "label", "Label"))
    if (length(c2) > 0) lab <- lab %>% rename(cluster = all_of(c2[1]))
  }

  if (!all(c("sample","cluster") %in% colnames(lab))) {
    stop("LABELS_CSV must contain columns: sample, cluster (or Sample_ID, cluster). Got: ",
         paste(colnames(lab), collapse = ", "))
  }

  lab <- lab %>%
    mutate(
      sample  = as.character(sample),
      cluster = str_trim(as.character(cluster))
    ) %>%
    filter(sample %in% colnames(zmat))

  if (nrow(lab) == 0) stop("No overlapping samples between labels and vst_matrix columns.")

  # ---- force 4 blocks: cluster order 0 -> 1 -> 2 -> 3 ----
  lab <- lab %>%
    mutate(cluster = factor(cluster, levels = c("0","1","2","3"))) %>%
    arrange(cluster, sample)

  # re-order columns by cluster blocks
  zmat <- zmat[, lab$sample, drop = FALSE]

  # annotation
  anno_col <- data.frame(cluster = lab$cluster)
  rownames(anno_col) <- lab$sample

  # gaps between cluster blocks
  cluster_sizes <- as.integer(table(lab$cluster))
  gaps_col <- cumsum(cluster_sizes)
  gaps_col <- gaps_col[-length(gaps_col)]

  message(">>> [Heatmap] using labels (n=", nrow(lab), "): ", labels_fp)
  message(">>> [Heatmap] cluster sizes: ",
          paste(names(table(lab$cluster)), as.integer(table(lab$cluster)), sep="=", collapse=" | "))
} else {
  message(">>> [Heatmap] no labels file found; heatmap will not be cluster-annotated.")
}

# ---------- plot ----------
tag <- paste0("TOP", topn, "_FDR", fdr_cut, "_LFC", lfc_cut)
out_png <- file.path(out_dir, paste0("marker_heatmap_", tag, ".png"))
out_pdf <- file.path(out_dir, paste0("marker_heatmap_", tag, ".pdf"))

png(out_png, width = 2400, height = 1400, res = 220)
pheatmap(
  zmat,
  annotation_col = anno_col,
  cluster_cols = FALSE,     # keep your block order
  gaps_col = gaps_col,      # draw 4 blocks
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  border_color = NA,
  clustering_method = "complete",
  main = paste0("Top markers heatmap (Top ", topn, "/cluster; FDR<=", fdr_cut, " |log2FC|>=", lfc_cut, ")")
)
dev.off()

pdf(out_pdf, width = 11, height = 7)
pheatmap(
  zmat,
  annotation_col = anno_col,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  border_color = NA,
  clustering_method = "complete",
  main = paste0("Top markers heatmap (Top ", topn, "/cluster; FDR<=", fdr_cut, " |log2FC|>=", lfc_cut, ")")
)
dev.off()

message(">>> [Heatmap] saved: ", out_png)
message(">>> [Heatmap] saved: ", out_pdf)
message(">>> [Heatmap] done.")
