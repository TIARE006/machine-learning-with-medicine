# ============================================================
# DESeq2 Differential Expression Analysis
# Dataset: GSE254877
# Comparison: Pancreas vs Colorectal
# Files:
#   data/RNA-seq/GSE254877_raw_counts_expression.csv
#   data/RNA-seq/GSE254877_pheno.csv
# ============================================================

suppressPackageStartupMessages(library(DESeq2))

setwd("/home/tiare/Desktop/machine learning with medicine")

counts_file <- "data/RNA-seq/GSE254877_raw_counts_expression.csv"
pheno_file  <- "data/RNA-seq/GSE254877_pheno.csv"

## 1. 读入 counts，清洗成整数矩阵 --------------------------
counts <- read.csv(counts_file,
                   header = TRUE,
                   row.names = 1,
                   check.names = FALSE)

# 去掉多余的 Unnamed 列
if (any(grepl("^Unnamed", colnames(counts)))) {
  counts <- counts[ , !grepl("^Unnamed", colnames(counts)), drop = FALSE]
}

# 转数值 + 填 NA = 0
counts <- as.matrix(counts)
mode(counts) <- "numeric"
counts[is.na(counts)] <- 0
storage.mode(counts) <- "integer"

cat("counts dim (genes x samples):", paste(dim(counts), collapse = " x "), "\n")

## 2. 读入 pheno，并按顺序对齐 -----------------------------
pheno <- read.csv(pheno_file,
                  row.names = 1,
                  check.names = FALSE)

if (ncol(counts) != nrow(pheno)) {
  stop("counts 列数和 pheno 行数不一致，无法按顺序对齐。")
}

pheno <- pheno[seq_len(ncol(counts)), , drop = FALSE]
rownames(pheno) <- colnames(counts)

## 3. 从 title 中提取分组信息（Colorectal / Pancreas） ------
title_vec <- as.character(pheno$title)
disease   <- trimws(sub(".*,", "", title_vec))  # 取逗号后部分

cat("Disease table:\n")
print(table(disease))

condition <- factor(disease,
                    levels = c("Colorectal", "Pancreas"))

cat("Group distribution:\n")
print(table(condition))

colData <- data.frame(
  row.names = colnames(counts),
  condition = condition
)

## 4. 构建 DESeq2 对象并运行 -------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = colData,
                              design    = ~ condition)

# 过滤低表达、
keep <- rowSums(counts(dds) >= 5) >= 1
dds  <- dds[keep, ]

cat("Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Pancreas", "Colorectal"))
res_df <- as.data.frame(res)

cat("DESeq2 summary:\n")
print(summary(res))

## 5. 严格阈值筛选差异基因（可以根据需要改） ---------------
resSig <- subset(res_df,
                 !is.na(padj) &
                 padj < 0.1 &
                 abs(log2FoldChange) > 0.5)

resSig <- resSig[order(resSig$log2FoldChange, decreasing = TRUE), ]

up_genes   <- sum(resSig$log2FoldChange >  1)
down_genes <- sum(resSig$log2FoldChange < -1)

cat("Significant DEGs:", nrow(resSig), "\n")
cat("  Upregulated:", up_genes, "\n")
cat("  Downregulated:", down_genes, "\n")
cat("Top genes:\n")
print(head(resSig, 10))

# ========= 创建差异结果目录 =========
out_dir <- "data/DEG_Pancreas_vs_Colorectal_DESeq2"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ========= 导出结果 =========
write.csv(resSig,
          file.path(out_dir,
                    "DEG_Pancreas_vs_Colorectal_DESeq2.csv"))

write.csv(counts(dds, normalized = TRUE),
          file.path(out_dir,
                    "Normalized_counts_Pancreas_vs_Colorectal.csv"))

cat("DONE. Results saved in:\n", out_dir, "\n")

