import os
import pandas as pd

BASE_DIR = "/home/tiare/Desktop/machine learning with medicine"

# ---- 文件路径 ----
cluster_file = os.path.join(
    BASE_DIR,
    "data", "integrated_results",
    "cluster_results_multiomics_RNA_miRNA_seed42.csv"
)

pheno_file = os.path.join(
    BASE_DIR,
    "data", "RNA-seq",
    "GSE254877_pheno.csv"   # 如果名称不同，根据你的文件改一下
)

print("===== Load Files =====")
print("Cluster file:", cluster_file)
print("Pheno file  :", pheno_file)

# ---- 读取文件 ----
clus = pd.read_csv(cluster_file)
pheno = pd.read_csv(pheno_file)

print("\nCluster file columns:", clus.columns.tolist())
print("Pheno file columns  :", pheno.columns.tolist())

# 你 pheno 表的病人 ID 列可能叫别的名字（比如 Sample_ID）
# 我猜你现在的 cluster 表有 patient_id，所以 pheno 可能也要改成相同列名
# 如果 pheno 里 patient ID 列不是 patient_id，请替换下面这一行

# 统一列名（按你的真实情况调整）
if "patient_id" not in pheno.columns:
    # 尝试找到可能的列名
    for col in pheno.columns:
        if "ID" in col or "sample" in col.lower():
            print(f"[INFO] 自动检测到可能的 patient ID 列: {col}")
            pheno = pheno.rename(columns={col: "patient_id"})
            break

# 合并
merged = pheno.merge(clus, on="patient_id", how="inner")

print("\n===== Merged Data =====")
print("Merged rows:", merged.shape[0])
print(merged.head())

# ---- 每个 cluster 的样本数 ----
print("\n===== Cluster Size =====")
print(merged["Cluster"].value_counts().sort_index())

# ---- 每个 cluster 的 cachexia / non-cachexia 统计 ----
# 假设 pheno 表里有一列 cachexia_status（如: Cachexia / Control）
possible_cols = [c for c in merged.columns if "cache" in c.lower()]

if possible_cols:
    status_col = possible_cols[0]
    print(f"\nDetected phenotype column: {status_col}")

    print("\n===== Cachexia Distribution per Cluster =====")
    print(merged.groupby("Cluster")[status_col].value_counts())
else:
    print("\n[WARN] 你的 pheno 表中没有检测到 cachexia 状态相关列。")
    print("请确认 pheno 表中 cachexia 列叫什么名字，并手动修改脚本。")

# ---- 生成论文方法部分英文描述 ----
n_samples = merged.shape[0]
n_clusters = merged["Cluster"].nunique()

method_text = f"""
We performed early-fusion multi-omics clustering using matched mRNA and miRNA 
expression profiles. After harmonizing sample identifiers between datasets, 
a total of {n_samples} subjects with both mRNA and miRNA data were retained.

For each omics layer, highly variable features were selected (top 2000 for mRNA 
and top 500 for miRNA), followed by z-score normalization and PCA dimensionality 
reduction (50 PCs for mRNA and 20 PCs for miRNA). The reduced features were 
concatenated to form a unified representation, and k-means clustering was applied 
across K = 3–5. The optimal number of clusters was determined using silhouette scores, 
resulting in {n_clusters} molecular subtypes.

This multi-omics subtype structure was subsequently used for downstream interpretation, 
including phenotype association, differential expression, pathway analysis, and ceRNA 
network integration.
"""

print("\n===== English Method Text (Copy into your report) =====")
print(method_text)

