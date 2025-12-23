#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)      # 如果不是人，改成 org.Mm.eg.db 等
  library(enrichplot)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
cat("ARGS:\n"); print(args)

if (length(args) < 1) {
  cat("Usage:\n")
  cat("  Rscript R/rna_pathway_enrichment_all.R <RUN_DIR> [FDR] [LFC] [TOPN]\n\n")
  cat("Outputs (per cluster):\n")
  cat("  <RUN_DIR>/pathways/GO_BP_cluster_<k>.csv + _dotplot.png\n")
  cat("  <RUN_DIR>/pathways/GO_MF_cluster_<k>.csv + _dotplot.png\n")
  cat("  <RUN_DIR>/pathways/GO_CC_cluster_<k>.csv + _dotplot.png\n")
  cat("  <RUN_DIR>/pathways/KEGG_cluster_<k>.csv  + _dotplot.png\n\n")
  cat("Outputs (compare):\n")
  cat("  <RUN_DIR>/pathways/GO_BP_cluster_compare.png\n")
  cat("  <RUN_DIR>/pathways/GO_MF_cluster_compare.png\n")
  cat("  <RUN_DIR>/pathways/GO_CC_cluster_compare.png\n")
  cat("  <RUN_DIR>/pathways/KEGG_cluster_compare.png\n")
  quit(status = 1)
}

run_dir <- args[1]
fdr_cut <- ifelse(length(args) >= 2, as.numeric(args[2]), 0.05)
lfc_cut <- ifelse(length(args) >= 3, as.numeric(args[3]), 1.0)
topn    <- ifelse(length(args) >= 4, as.integer(args[4]), 15)

degs_dir <- file.path(run_dir, "degs_deseq2")
out_dir  <- file.path(run_dir, "pathways")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message(">>> [Enrich] run_dir: ", run_dir)
message(">>> [Enrich] degs_dir: ", degs_dir)
message(">>> [Enrich] out_dir : ", out_dir)
message(">>> [Enrich] thresholds: FDR<=", fdr_cut, " |log2FC|>=", lfc_cut)
message(">>> [Enrich] TOPN: ", topn)

if (!dir.exists(degs_dir)) {
  stop("degs_dir not found: ", degs_dir, "\nDid you run DESeq2 first?")
}

deg_files <- list.files(degs_dir, pattern="^DEG_cluster_\\d+_vs_rest\\.csv$", full.names=TRUE)
if (length(deg_files) == 0) stop("No DEG files found in: ", degs_dir)

# ---- filter non-standard symbols (optional but improves mapping) ----
is_good_symbol <- function(x) {
  x <- as.character(x)
  ok <- !is.na(x) & x != ""
  ok <- ok & !str_detect(x, "^(AC\\d+|AL\\d+|FP\\d+|LOC\\d+|LINC\\d+)")
  ok <- ok & str_detect(x, "^[A-Za-z0-9][A-Za-z0-9\\-\\.]*$")
  ok
}

# ---- SYMBOL -> ENTREZID (safe) ----
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

# ---- helper: write + dotplot ----
save_enrich <- function(res, prefix, k, title_txt) {
  if (is.null(res)) return(FALSE)
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(FALSE)

  out_csv <- file.path(out_dir, paste0(prefix, "_cluster_", k, ".csv"))
  write_csv(df, out_csv)

  p <- dotplot(res, showCategory = topn, font.size = 10) +
    ggtitle(title_txt) +
    theme(plot.title = element_text(size = 14, face = "bold"))

  out_png <- file.path(out_dir, paste0(prefix, "_cluster_", k, "_dotplot.png"))
  ggsave(out_png, p, width = 9, height = 6, dpi = 400)

  message("    saved: ", out_csv)
  message("    saved: ", out_png)
  TRUE
}

# -------------------------
# A) per-cluster enrichment (GO BP/MF/CC + KEGG)
# -------------------------
for (fp in deg_files) {
  bn <- basename(fp)
  m  <- str_match(bn, "DEG_cluster_(\\d+)_vs_rest\\.csv")
  k  <- m[,2]
  if (is.na(k) || k == "") {
    message("\n>>> [Enrich] [skip] cannot parse cluster id from filename: ", bn)
    next
  }

  message("\n>>> [Enrich] Cluster ", k, " file: ", bn)

  df <- suppressMessages(read_csv(fp, show_col_types = FALSE))

  # robust column mapping
  if (!("gene" %in% colnames(df))) {
    g2 <- intersect(colnames(df), c("Gene", "SYMBOL", "symbol"))
    if (length(g2) == 0) stop("No gene column in: ", fp)
    df <- df %>% rename(gene = all_of(g2[1]))
  }
  if (!("FDR" %in% colnames(df))) {
    f2 <- intersect(colnames(df), c("padj", "p_adj", "adj_p", "qvalue"))
    if (length(f2) == 0) stop("No FDR/padj column in: ", fp)
    df <- df %>% rename(FDR = all_of(f2[1]))
  }
  if (!("log2FC" %in% colnames(df))) {
    l2 <- intersect(colnames(df), c("log2FoldChange", "logFC", "log2_fc"))
    if (length(l2) == 0) stop("No log2FC/log2FoldChange column in: ", fp)
    df <- df %>% rename(log2FC = all_of(l2[1]))
  }

  df <- df %>%
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

  # ---- GO BP/MF/CC ----
  for (ont in c("BP","MF","CC")) {
    ego <- tryCatch(
      enrichGO(
        gene          = entrez,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        readable      = TRUE
      ),
      error = function(e) {
        message("    [warn] enrichGO(", ont, ") failed: ", conditionMessage(e))
        return(NULL)
      }
    )

    prefix <- paste0("GO_", ont)
    ok <- save_enrich(
      ego, prefix, k,
      paste0("GO ", ont, " enrichment (Cluster ", k, " vs rest)")
    )
    if (!ok) message("    [warn] no GO ", ont, " terms enriched for cluster ", k)
  }

  # ---- KEGG (human = 'hsa'; mouse = 'mmu') ----
  ekegg <- tryCatch(
    enrichKEGG(
      gene          = entrez,
      organism      = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    ),
    error = function(e) {
      message("    [warn] enrichKEGG failed: ", conditionMessage(e))
      return(NULL)
    }
  )

  # (optional) convert KEGG gene IDs to readable SYMBOL (needs OrgDb)
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg <- tryCatch(setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                     error = function(e) ekegg)
  }

  okk <- save_enrich(
    ekegg, "KEGG", k,
    paste0("KEGG pathway enrichment (Cluster ", k, " vs rest)")
  )
  if (!okk) message("    [warn] no KEGG pathways enriched for cluster ", k)
}

# -------------------------
# B) compare plot for each prefix
# -------------------------
make_compare <- function(prefix, title_main) {
  csvs <- list.files(out_dir, pattern = paste0("^", prefix, "_cluster_\\d+\\.csv$"), full.names = TRUE)
  if (length(csvs) < 2) {
    message(">>> [Compare] skip ", prefix, ": need >=2 cluster csv files.")
    return(invisible(FALSE))
  }

  cmp <- lapply(csvs, function(fp) {
    k <- str_match(basename(fp), paste0("^", prefix, "_cluster_(\\d+)\\.csv$"))[,2]
    df <- suppressMessages(read_csv(fp, show_col_types = FALSE))

    # clusterProfiler output columns
    if (!("Description" %in% colnames(df))) {
      if ("Term" %in% colnames(df)) df <- df %>% rename(Description = Term)
    }
    if (!("p.adjust" %in% colnames(df))) {
      if ("p.adjust" %in% colnames(df)) { } else stop("Missing p.adjust in: ", fp)
    }
    if (!("Count" %in% colnames(df))) {
      # some outputs use 'GeneRatio'/'BgRatio' but usually Count exists
      if (!("Count" %in% colnames(df))) stop("Missing Count in: ", fp)
    }

    df %>%
      transmute(
        Cluster = paste0("C", k),
        Description = as.character(Description),
        p.adjust = as.numeric(p.adjust),
        Count = as.numeric(Count)
      ) %>%
      filter(!is.na(p.adjust)) %>%
      arrange(p.adjust) %>%
      slice_head(n = topn)
  }) %>% bind_rows()

  if (nrow(cmp) == 0) {
    message(">>> [Compare] skip ", prefix, ": empty after topN filtering.")
    return(invisible(FALSE))
  }

  term_order <- cmp %>%
    group_by(Description) %>%
    summarise(best = min(p.adjust, na.rm = TRUE), .groups="drop") %>%
    arrange(best) %>%
    pull(Description)

  cmp <- cmp %>%
    mutate(
      Description = factor(Description, levels = rev(term_order)),
      neglog10FDR = -log10(p.adjust + 1e-300)
    )

  out_png <- file.path(out_dir, paste0(prefix, "_cluster_compare.png"))

  p_cmp <- ggplot(cmp, aes(x = Cluster, y = Description)) +
    geom_point(aes(size = Count, color = neglog10FDR), alpha = 0.9) +
    scale_size_continuous(name = "Count") +
    scale_color_gradient(name = "-log10(FDR)") +
    labs(
      title = paste0(title_main, " (top ", topn, " per cluster)"),
      subtitle = paste0("DEG thresholds: FDR<=", fdr_cut, " |log2FC|>=", lfc_cut),
      x = "Cluster",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9),
      panel.grid.minor = element_blank()
    )

  ggsave(out_png, p_cmp,
         width = 10.5,
         height = max(6, 0.25 * length(unique(cmp$Description))),
         dpi = 400)
  message(">>> [Compare] saved: ", out_png)
  TRUE
}

message("\n>>> [Compare] generating comparison plots...")
make_compare("GO_BP",  "GO BP enrichment comparison")
make_compare("GO_MF",  "GO MF enrichment comparison")
make_compare("GO_CC",  "GO CC enrichment comparison")
make_compare("KEGG",   "KEGG pathway enrichment comparison")

message("\n>>> [Enrich] done.")

