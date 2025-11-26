setwd("/home/tiare/Desktop/machine learning with medicine")

pheno <- read.csv("data/small RNA-seq/GSE254878_pheno.csv",
                  row.names = 1,
                  check.names = FALSE)

cat("=== characteristics_ch1 ===\n")
print(table(pheno$characteristics_ch1))

cat("\n=== characteristics_ch1.1 ===\n")
print(table(pheno$characteristics_ch1.1))

cat("\n=== title 前10行 ===\n")
print(head(pheno$title, 10))
