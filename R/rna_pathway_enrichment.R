#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)      # 如果不是人，改成 org.Mm.eg.db
  library(enrichplot)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
cat("ARGS:\n"); print(args)

if (length(args) < 1) {
  cat("Usage:\n  Rscript R/rna_pathway_enrichment.R <RUN_DIR> [FDR] [LFC]\n")
  quit(status = 1)
}

run_dir <- args[1]
fdr_cut <- ifelse(length(args) >= 2, as.numeric(args[2]), 0.05)
lfc_cut <- ifelse(length(args) >= 3, as.numeric(args[3]), 1.0)

degs_dir <- file.path(run_dir, "degs_deseq2")
out_dir  <- file.path(run_dir, "pathways")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message(">>> [Pathway] run_dir: ", run_dir)
message(">>> [Pathway] degs_dir: ", degs_dir)
message(">>> [Pathway] out_dir : ", out_dir)
message(">>> [Pathway] thresholds: FDR<=", fdr_cut, " |log2FC|>=", lfc_cut)

if (!dir.exists(degs_dir)) {
  stop("degs_dir not found: ", degs_dir, "\nDid you run DESeq2 first?")
}

deg_files <- list.files(degs_dir, pattern="^DEG_cluster_.*_vs_rest\\.csv$", full.names=TRUE)
if (length(deg_files) == 0) stop("No DEG files found in: ", degs_dir)

# 过滤一些常见“非标准 SYMBOL”的基因名（会导致 bitr 映射率很低）
is_good_symbol <- function(x) {
  x <- as.character(x)
  ok <- !is.na(x) & x != ""
  ok <- ok & !str_detect(x, "^(AC\\d+|AL\\d+|FP\\d+|LOC\\d+|LINC\\d+)")
  ok <- ok & str_detect(x, "^[A-Za-z0-9][A-Za-z0-9\\-\\.]*$")  # 基本字符集
  ok
}

# helper: symbol -> ENTREZ (带保护，不会因为 bitr 崩掉)
sym2entrez <- function(symbols) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[is_good_symbol(symbols)]
  if (length(symbols) == 0) return(character(0))

  mp <- tryCatch(
    clusterProfiler::bitr(
      symbols,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    ),
    error = function(e) {
      message("    [warn] bitr mapping failed: ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(mp) || nrow(mp) == 0) return(character(0))
  unique(mp$ENTREZID)
}

for (fp in deg_files) {
  bn <- basename(fp)
  m  <- str_match(bn, "DEG_cluster_(\\d+)_vs_rest\\.csv")
  k  <- m[,2]

  if (is.na(k) || k == "") {
    message("\n>>> [Pathway] [skip] cannot parse cluster id from filename: ", bn)
    next
  }

  message("\n>>> [Pathway] Cluster ", k, " file: ", bn)

  df <- read_csv(fp, show_col_types = FALSE) %>%
    mutate(gene = as.character(gene)) %>%
    filter(!is.na(FDR), !is.na(log2FC))

  df_sig <- df %>% filter(FDR <= fdr_cut, abs(log2FC) >= lfc_cut)
  message("    significant genes (raw): ", nrow(df_sig))

  if (nrow(df_sig) < 20) {
    message("    [skip] too few genes for enrichment.")
    next
  }

  entrez <- sym2entrez(df_sig$gene)
  message("    mapped ENTREZ (after filtering): ", length(entrez))

  if (length(entrez) < 20) {
    message("    [skip] too few mapped ENTREZ IDs.")
    next
  }

  ego <- tryCatch(
    enrichGO(
      gene          = entrez,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    ),
    error = function(e) {
      message("    [warn] enrichGO failed: ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("    [warn] no GO terms enriched.")
    next
  }

  out_csv <- file.path(out_dir, paste0("GO_BP_cluster_", k, ".csv"))
  write_csv(as.data.frame(ego), out_csv)
  message("    saved: ", out_csv)

  p <- dotplot(ego, showCategory = 15, font.size = 10) +
    ggtitle(paste0("GO BP enrichment (Cluster ", k, " vs rest)")) +
    theme(plot.title = element_text(size = 14, face = "bold"))

  out_png <- file.path(out_dir, paste0("GO_BP_cluster_", k, "_dotplot.png"))
  ggsave(out_png, p, width = 9, height = 6, dpi = 400)
  message("    saved: ", out_png)
}

message("\n>>> [Pathway] done.")
