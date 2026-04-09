#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage:\n  Rscript R/differential_expression/merge_deg_tables.R <RUN_DIR> [OUT_NAME]\n\n")
  cat("Reads:\n  <RUN_DIR>/degs_deseq2/DEG_cluster_*_vs_rest.csv\n")
  cat("Writes:\n  <RUN_DIR>/degs_deseq2/<OUT_NAME>.csv (Excel-readable)\n")
  cat("  <RUN_DIR>/degs_deseq2/<OUT_NAME>.xlsx (if openxlsx available)\n")
  quit(status = 1)
}

run_dir <- args[1]
out_name <- ifelse(length(args) >= 2, args[2], "DEG_all_clusters_merged")

degs_dir <- file.path(run_dir, "degs_deseq2")
if (!dir.exists(degs_dir)) stop("Missing degs_dir: ", degs_dir)

files <- list.files(degs_dir, pattern = "^DEG_cluster_\\d+_vs_rest\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No DEG_cluster_*_vs_rest.csv found in: ", degs_dir)

message(">>> [MERGE] run_dir : ", run_dir)
message(">>> [MERGE] degs_dir: ", degs_dir)
message(">>> [MERGE] files  : ", length(files))

merged <- lapply(files, function(fp) {
  k <- str_match(basename(fp), "^DEG_cluster_(\\d+)_vs_rest\\.csv$")[,2]
  df <- suppressMessages(read_csv(fp, show_col_types = FALSE))

  # 标准化列名（兼容你之前各种脚本输出）
  # 期望列：gene, log2FC, p_value, FDR, baseMean, lfcSE, stat
  if (!("gene" %in% colnames(df))) {
    g2 <- intersect(colnames(df), c("Gene", "SYMBOL", "symbol"))
    if (length(g2) > 0) df <- df %>% rename(gene = all_of(g2[1]))
  }
  if (!("log2FC" %in% colnames(df))) {
    l2 <- intersect(colnames(df), c("log2FoldChange", "logFC"))
    if (length(l2) > 0) df <- df %>% rename(log2FC = all_of(l2[1]))
  }
  if (!("FDR" %in% colnames(df))) {
    f2 <- intersect(colnames(df), c("padj", "p_adj", "qvalue"))
    if (length(f2) > 0) df <- df %>% rename(FDR = all_of(f2[1]))
  }
  if (!("p_value" %in% colnames(df))) {
    p2 <- intersect(colnames(df), c("pvalue", "p_value", "p.val"))
    if (length(p2) > 0) df <- df %>% rename(p_value = all_of(p2[1]))
  }

  df %>%
    mutate(
      Cluster = as.integer(k),
      gene = as.character(gene)
    ) %>%
    relocate(Cluster, gene)
}) %>% bind_rows()

# 排序：先按 Cluster，再按 FDR（小到大），再按 |log2FC|（大到小）
if ("FDR" %in% colnames(merged)) {
  merged <- merged %>% arrange(Cluster, FDR, desc(abs(log2FC)))
} else {
  merged <- merged %>% arrange(Cluster, desc(abs(log2FC)))
}

out_csv  <- file.path(degs_dir, paste0(out_name, ".csv"))
write_csv(merged, out_csv)
message(">>> [MERGE] saved CSV : ", out_csv, " (rows=", nrow(merged), ")")

# 可选输出 xlsx（如果 openxlsx 存在）
if (requireNamespace("openxlsx", quietly = TRUE)) {
  out_xlsx <- file.path(degs_dir, paste0(out_name, ".xlsx"))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "DEG_merged")
  openxlsx::writeData(wb, "DEG_merged", merged)
  openxlsx::freezePane(wb, "DEG_merged", firstRow = TRUE, firstCol = 2)
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  message(">>> [MERGE] saved XLSX: ", out_xlsx)
} else {
  message(">>> [MERGE] openxlsx not installed; skip XLSX. (CSV is Excel-readable)")
  message(">>> [MERGE] To enable XLSX: R -q -e 'install.packages(\"openxlsx\")'")
}

message(">>> [MERGE] done.")
