#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
})

# ============ 路径 & 参数 ============
base_dir <- "/home/tiare/Desktop/machine learning with medicine"
seed     <- 42  # 跟 Python RANDOM_STATE 保持一致

# SNF 聚类结果（已经跑过 run_snf_v2_pipeline 生成的）
cluster_file <- file.path(
  base_dir,
  "data", "integrated_results",
  sprintf("cluster_results_SNF_RNA_miRNA_seed%d.csv", seed)
)

cat(">>> [R-DEG] Load SNF clusters from:", cluster_file, "\n")
clust <- read_csv(cluster_file, show_col_types = FALSE)

# 列名是 Sample_ID + Cluster_SNF，这里统一成 Sample_ID / Cluster
clust <- clust %>%
  rename(
    Sample_ID = Sample_ID,
    Cluster   = Cluster_SNF
  )

# 只做 Cluster 0 vs 1（跟你 Python 现在的 group1=0, group2=1 一致）
clust_01 <- clust %>% filter(Cluster %in% c(0, 1))
clust_01$Sample_ID <- as.character(clust_01$Sample_ID)
clust_01$Cluster   <- factor(clust_01$Cluster, levels = c(0, 1))

cat(">>> [R-DEG] SNF cluster counts (0 vs 1):\n")
print(table(clust_01$Cluster))


# ============ 简化样本名：和 Python 逻辑一致 ============
# Python 日志里简化后变成 2066, 3005 这种，所以这里用相同规则：
# 取最后一段数字（在最后一个点后面），可选尾巴 "_R1.bam"
simplify_id <- function(x) {
  # e.g. "NS.xxx.2066_R1.bam" or "NS.xxx.3016" 或者已经是 "2066"
  x <- as.character(x)
  out <- sub(".*\\.(\\d+)(?:_R1\\.bam)?$", "\\1", x, perl = TRUE)
  # 如果没匹配上（比如原本就是 "2066"），sub 不会改，直接返回原字符串
  out
}


# ============ 通用函数：对一种 omics 跑 DESeq2 ============
run_deseq2_for_omics <- function(
  data_type = c("RNA", "smallRNA", "lncRNA")
) {
  data_type <- match.arg(data_type)

  if (data_type == "RNA") {
    expr_file <- file.path(
      base_dir, "data", "RNA-seq", "raw",
      "GSE254877_raw_counts_expression.csv"
    )
    deg_dir <- file.path(base_dir, "data", "RNA-seq", "deg_SNF")
  } else if (data_type == "smallRNA") {
    expr_file <- file.path(
      base_dir, "data", "small RNA-seq", "raw",
      "GSE254878_smallRNAs_raw_counts_expression.csv"
    )
    deg_dir <- file.path(base_dir, "data", "small RNA-seq", "deg_SNF")
  } else if (data_type == "lncRNA") {
    expr_file <- file.path(
      base_dir, "data", "lncRNA-seq", "raw",
      "GSE254877_lncRNA_raw_counts_expression.csv"
    )
    deg_dir <- file.path(base_dir, "data", "lncRNA-seq", "deg_SNF")
  }

  if (!dir.exists(deg_dir)) dir.create(deg_dir, recursive = TRUE)

  cat("\n==============================\n")
  cat(">>> [R-DEG]", data_type, "\n")
  cat("    Expression file:", expr_file, "\n")

  if (!file.exists(expr_file)) {
    stop("Expression file not found: ", expr_file)
  }

  # ---- 读取表达矩阵 ----
  df <- read_csv(expr_file, show_col_types = FALSE)

  # 按你之前的说明：第一行是注释，删掉
  df <- df[-1, ]

  # 第一列是基因 / miRNA / lncRNA ID
  gene_col <- colnames(df)[1]
  gene_ids <- df[[gene_col]]

  df <- df %>%
    select(-all_of(gene_col)) %>%
    as.data.frame()

  # smallRNA 有 "type" 列，删掉
  if ("type" %in% colnames(df)) {
    df$type <- NULL
  }

  # 全部转 numeric
  df[] <- lapply(df, function(x) as.numeric(x))

  # 行 = gene, 列 = sample（DESeq2 标准形状）
  rownames(df) <- gene_ids
  raw_counts <- as.matrix(df)

  # ---- 简化样本名，与 SNF Sample_ID 对齐 ----
  raw_sample_ids <- colnames(raw_counts)
  simp_ids       <- simplify_id(raw_sample_ids)

  colnames(raw_counts) <- simp_ids

  # 和 SNF 的 0/1 样本取交集
  snf_ids    <- clust_01$Sample_ID
  common_ids <- intersect(simp_ids, snf_ids)

  cat("    #raw samples:", length(simp_ids),
      ", #common with SNF (cluster 0/1):", length(common_ids), "\n")

  if (length(common_ids) < 4) {
    stop("[R-DEG] Too few common samples for ", data_type,
         " after matching with SNF 0/1.")
  }

  # 只保留这些公共样本的列
  counts_use <- raw_counts[, common_ids, drop = FALSE]

  # ---- 构造 colData（样本分组） ----
  clust_sub <- clust_01[match(common_ids, clust_01$Sample_ID), ]
  # 检查顺序完全一致
  stopifnot(all(as.character(clust_sub$Sample_ID) == colnames(counts_use)))

  coldata <- data.frame(
    row.names = common_ids,
    group = factor(clust_sub$Cluster, levels = c(0, 1))
  )

  cat("    Group table:\n")
  print(table(coldata$group))

  # ---- 构建 DESeq2 对象 ----
  dds <- DESeqDataSetFromMatrix(
    countData = counts_use,
    colData   = coldata,
    design    = ~ group
  )

  # 过滤低表达基因（行和 >= 10）
  keep <- rowSums(counts(dds)) >= 10
  dds  <- dds[keep, ]
  cat("    #Genes after filtering (rowSums >= 10):", nrow(dds), "\n")

  # ---- 跑 DESeq2 ----
  dds <- DESeq(dds)

  # 对比：group 0 vs group 1（和 Python 的 0 vs 1 一致）
  res <- results(dds, contrast = c("group", "0", "1"))

  # 收缩 log2FC（用 normal，避免 apeglm 额外依赖）
  res <- lfcShrink(dds,
                   contrast = c("group", "0", "1"),
                   res      = res,
                   type     = "normal")

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)

  # 去掉 padj 为 NA 的行，按 FDR 排序
  res_df <- res_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj)

  # 输出成 Python 友好的列名
  out_df <- res_df %>%
    transmute(
      gene        = gene,
      log2FC      = log2FoldChange,
      p_value     = pvalue,
      FDR         = padj,
      baseMean    = baseMean,
      lfcSE       = lfcSE,
      stat        = stat
    )

  # 文件名跟 Python 的 save_deg_tables 对齐：
  # DEG_full_RNA_SNF_seed42.csv 这种格式
  out_file <- file.path(
    deg_dir,
    sprintf("DEG_full_%s_SNF_seed%d.csv", data_type, seed)
  )

  write_csv(out_df, out_file)
  cat(">>> [R-DEG] Save:", out_file, "\n")
  cat("    #rows:", nrow(out_df), "\n")
}


# ============ 顺序跑三种组学 ============
run_deseq2_for_omics("RNA")
run_deseq2_for_omics("smallRNA")
run_deseq2_for_omics("lncRNA")

cat("\n>>> [R-DEG] All omics finished.\n")

