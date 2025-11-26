# ============================================================
# DESeq2: GSE254878 (small RNA)
# Comparison: Pancreas vs Colorectal
# ============================================================

suppressPackageStartupMessages(library(DESeq2))

# 设置工作目录
setwd("/home/tiare/Desktop/machine learning with medicine")

counts_file <- "data/small RNA-seq/GSE254878_smallRNAs_raw_counts_expression.csv"
pheno_file  <- "data/small RNA-seq/GSE254878_pheno.csv"

cat("counts file : ", counts_file, "\n")
cat("pheno  file : ", pheno_file,  "\n")

## 1. 读入 counts -------------------------------------------
counts <- read.csv(counts_file,
                   header = TRUE,
                   row.names = 1,
                   check.names = FALSE)

# 去掉 Unnamed 列（如果有）
if (any(grepl("^Unnamed", colnames(counts)))) {
  counts <- counts[, !grepl("^Unnamed", colnames(counts)), drop = FALSE]
}

counts <- as.matrix(counts)
mode(counts) <- "numeric"
counts[is.na(counts)] <- 0      # NA 填 0
storage.mode(counts) <- "integer"

cat("counts dim (miRNA x columns):",
    paste(dim(counts), collapse = " x "), "\n")

## 2. 读入 pheno ---------------------------------------------
pheno <- read.csv(pheno_file,
                  row.names = 1,
                  check.names = FALSE)

cat("pheno rows:", nrow(pheno), "\n")

## 3. 处理 85 vs 84 的问题：如果多 1 列，就裁掉 1 列 ----------
if (ncol(counts) == nrow(pheno) + 1) {
  cat("⚠ counts 比 pheno 多 1 列，自动只保留前 ", nrow(pheno), " 列用于分析。\n")
  counts <- counts[, 1:nrow(pheno), drop = FALSE]
}

if (ncol(counts) != nrow(pheno)) {
  stop("裁剪后 counts 列数(", ncol(counts),
       ") 和 pheno 行数(", nrow(pheno),
       ") 仍然不一致，请检查数据。")
}

## 4. 按“顺序”对齐 -----------------------------------------
# 假设导出时 pheno 和 counts 的顺序是一致的
pheno_aligned <- pheno[seq_len(ncol(counts)), , drop = FALSE]
rownames(pheno_aligned) <- colnames(counts)

cat("对齐后 counts 维度:",
    paste(dim(counts), collapse = " x "), "\n")
cat("对齐后 pheno_aligned 行数:", nrow(pheno_aligned), "\n")

## 5. 从 title 中提取疾病类型 -------------------------------
title_vec <- as.character(pheno_aligned$title)
disease   <- trimws(sub(".*,", "", title_vec))  # 取逗号后部分: Colorectal / Pancreas

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

## 6. 构建 DESeq2 对象并运行 -------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = colData,
                              design    = ~ condition)

# small RNA 稍微宽松一点的基因过滤
keep <- rowSums(counts(dds) >= 5) >= 2
dds  <- dds[keep, ]

cat("miRNAs after filtering:", nrow(dds), "\n")

# 检查每个样本的总 reads 数（防止某个样本全是 0）
libsize <- colSums(counts(dds))
cat("Library sizes (per sample):\n")
print(libsize)

if (any(libsize == 0)) {
  cat("⚠ 发现总 counts 为 0 的样本，将其从分析中删除:\n")
  print(names(libsize[libsize == 0]))
  dds <- dds[, libsize > 0]
}

# 关键：使用 poscounts 方式估计 size factor，适合稀疏 small RNA
dds <- DESeq(dds, sfType = "poscounts")

res <- results(dds, contrast = c("condition", "Pancreas", "Colorectal"))
res_df <- as.data.frame(res)


## 7. 筛差异 miRNA -----------------------------------------
resSig <- subset(res_df,
                 !is.na(padj) &
                 padj < 0.05 &
                 abs(log2FoldChange) > 1)

resSig <- resSig[order(resSig$log2FoldChange, decreasing = TRUE), ]

up_genes   <- sum(resSig$log2FoldChange >  1)
down_genes <- sum(resSig$log2FoldChange < -1)

cat("Significant DE miRNAs:", nrow(resSig), "\n")
cat("  Upregulated:", up_genes, "\n")
cat("  Downregulated:", down_genes, "\n")
cat("Top miRNAs:\n")
print(head(resSig, 10))

## 8. 导出结果 ----------------------------------------------
out_dir <- "data/DEG_smRNA_Pancreas_vs_Colorectal_DESeq2"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write.csv(resSig,
          file.path(out_dir,
                    "DEG_GSE254878_smRNA_DESeq2_Pancreas_vs_Colorectal.csv"))

write.csv(counts(dds, normalized = TRUE),
          file.path(out_dir,
                    "Normalized_counts_GSE254878_smRNA.csv"))

cat("DONE. Results saved in:\n", out_dir, "\n")


