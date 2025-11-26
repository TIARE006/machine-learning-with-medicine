# know.R  —— 从 GSE254878 读取样本注释，看看能不能做恶病质

suppressPackageStartupMessages(library(GEOquery))

# 根目录（注意这里是整个项目目录，不是 small RNA 那一层）
setwd("/home/tiare/Desktop/machine learning with medicine")

cat("=== 通过 GEO 直接获取 GSE254878 ===\n")
gse_list <- getGEO("GSE254878", GSEMatrix = TRUE)

# 通常只有一个表达矩阵
eset  <- gse_list[[1]]
pheno <- pData(eset)

cat("=== pheno 的列名 ===\n")
print(colnames(pheno))

cat("\n=== title 前 20 行 ===\n")
print(head(pheno$title, 20))

cat("\n=== 包含 'characteristics' 的列 ===\n")
char_cols <- grep("characteristics", colnames(pheno), value = TRUE)
print(char_cols)

if (length(char_cols) > 0) {
  cat("\n=== characteristics 列前 20 行 ===\n")
  print(head(pheno[, char_cols, drop = FALSE], 20))
} else {
  cat("\n没有发现任何 'characteristics' 相关列。\n")
}

cat("\n=== 粗暴搜一下是否有 cachexia 相关字样 ===\n")
has_cachexia <- apply(pheno, 2, function(col) {
  any(grepl("cachexia|wasting|weight", col, ignore.case = TRUE))
})
print(has_cachexia)



# setwd("/home/tiare/Desktop/machine learning with medicine")

# pheno <- read.csv("data/small RNA-seq/GSE254878_pheno.csv",
#                   row.names = 1,
#                   check.names = FALSE)

# cat("=== characteristics_ch1 ===\n")
# print(table(pheno$characteristics_ch1))

# cat("\n=== characteristics_ch1.1 ===\n")
# print(table(pheno$characteristics_ch1.1))

# cat("\n=== title 前10行 ===\n")
# print(head(pheno$title, 10))
