import os
import numpy as np
import pandas as pd

from scipy.stats import ttest_ind
import gseapy as gp   # pip install gseapy


# =========================
# 配置区（按需要改）
# =========================
DATA_TYPE = "RNA"          # "smallRNA" 或 "RNA"
RANDOM_STATE = 42          # 跟聚类脚本一致即可

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# 原始表达矩阵 & 聚类结果文件
if DATA_TYPE == "smallRNA":
    expr_file = os.path.join(
        BASE_DIR,
        "data",
        "small RNA-seq",
        "raw",   # ← 修正路径
        "GSE254878_smallRNAs_raw_counts_expression.csv"
    )
    cluster_file = os.path.join(
        BASE_DIR,
        "data",
        "small RNA-seq",
        "clustering",    # ← 聚类结果现在在 clustering 目录
        f"cluster_results_smallRNA_seed{RANDOM_STATE}.csv"
    )
elif DATA_TYPE == "RNA":
    expr_file = os.path.join(
        BASE_DIR,
        "data",
        "RNA-seq",
        "raw",   # ← 修正路径
        "GSE254877_raw_counts_expression.csv"
    )
    cluster_file = os.path.join(
        BASE_DIR,
        "data",
        "RNA-seq",
        "clustering",    # ← 聚类结果现在在 clustering 目录
        f"cluster_results_RNA_seed{RANDOM_STATE}.csv"
    )
else:
    raise ValueError("DATA_TYPE 必须是 'smallRNA' 或 'RNA'")



# =========================
# 1. 读取表达矩阵 + 聚类结果
# =========================
print(f"[INFO] Using expression file: {expr_file}")
print(f"[INFO] Using cluster file   : {cluster_file}")

df = pd.read_csv(expr_file, low_memory=False)

# 删除第一行描述行
df = df.drop(index=0)

# 第一列作为基因 ID
gene_col = df.columns[0]
df = df.set_index(gene_col)

# smallRNA 有 type 列，删掉
if "type" in df.columns:
    df = df.drop(columns=["type"])

# 转数值并清理
df = df.apply(pd.to_numeric, errors="coerce")
df = df.dropna(how="all")
df = df.fillna(0)

# 转置：行=样本，列=基因
X = df.T
print(f"[INFO] Expression matrix shape: {X.shape}")

# 读取聚类结果
clusters = pd.read_csv(cluster_file)
# 要求列名：Sample_ID, Cluster
clusters = clusters.set_index("Sample_ID")

# 与表达矩阵对齐
clusters = clusters.loc[X.index]
cluster_labels = clusters["Cluster"].values

print(f"[INFO] Loaded cluster labels, k = {len(np.unique(cluster_labels))}")
print(f"[INFO] Cluster counts:\n{clusters['Cluster'].value_counts()}")


# =========================
# 2. Benjamini-Hochberg FDR 校正
# =========================
def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """
    BH FDR 校正，返回与 pvals 同顺序的 FDR 数组
    """
    pvals = np.asarray(pvals)
    n = pvals.size
    order = np.argsort(pvals)
    ranked_pvals = pvals[order]

    fdr = ranked_pvals * n / (np.arange(1, n + 1))
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr = np.clip(fdr, 0, 1.0)

    fdr_corrected = np.empty_like(fdr)
    fdr_corrected[order] = fdr
    return fdr_corrected


# =========================
# 3. 差异基因分析（Cluster 0 vs 1）
# =========================
def compute_deg(X_df: pd.DataFrame,
                labels: np.ndarray,
                group1: int = 0,
                group2: int = 1) -> pd.DataFrame:
    """
    X_df: 行=样本，列=基因
    labels: 每个样本的 cluster 标签
    group1 vs group2 做差异
    """
    mask1 = (labels == group1)
    mask2 = (labels == group2)

    X1 = X_df.loc[mask1]
    X2 = X_df.loc[mask2]

    print(f"[INFO] DEG: comparing Cluster {group1} (n={mask1.sum()}) "
          f"vs Cluster {group2} (n={mask2.sum()})")

    p_vals = []
    log2fc = []

    for gene in X_df.columns:
        g1 = X1[gene].values
        g2 = X2[gene].values

        # Welch t-test
        stat, p = ttest_ind(g1, g2, equal_var=False)
        p_vals.append(p)

        # log2 fold change (group1 / group2)
        fc = (g1.mean() + 1e-9) / (g2.mean() + 1e-9)
        log2fc.append(np.log2(fc))

    p_vals = np.array(p_vals)
    fdr = benjamini_hochberg(p_vals)

    deg = pd.DataFrame({
        "gene": X_df.columns,
        "log2FC": log2fc,
        "p_value": p_vals,
        "FDR": fdr
    })

    deg = deg.sort_values("FDR")
    return deg


deg = compute_deg(X, cluster_labels, group1=0, group2=1)

# =========================
# 4. 输出目录（按功能分类）
# =========================
if DATA_TYPE == "smallRNA":
    base_out = os.path.join(BASE_DIR, "data", "small RNA-seq")
else:
    base_out = os.path.join(BASE_DIR, "data", "RNA-seq")

deg_dir = os.path.join(base_out, "deg")
pathway_dir = os.path.join(base_out, "pathway")

os.makedirs(deg_dir, exist_ok=True)
os.makedirs(pathway_dir, exist_ok=True)

# ---- 4.1 保存完整 DEG 表 ----
deg_file = os.path.join(deg_dir, f"DEG_full_{DATA_TYPE}_seed{RANDOM_STATE}.csv")
deg.to_csv(deg_file, index=False)
print(f"[SAVE] Full DEG table -> {deg_file}")

# ---- 4.2 按阈值筛选显著差异基因 ----
FDR_CUT = 0.05
LOG2FC_CUT = 1.0

deg_sig = deg[(deg["FDR"] < FDR_CUT) & (deg["log2FC"].abs() > LOG2FC_CUT)]
deg_sig_file = os.path.join(
    deg_dir,
    f"DEG_sig_{DATA_TYPE}_FDR{FDR_CUT}_log2FC{LOG2FC_CUT}_seed{RANDOM_STATE}.csv"
)
deg_sig.to_csv(deg_sig_file, index=False)
print(f"[SAVE] Significant DEGs -> {deg_sig_file}")
print(f"[INFO] #Significant DEGs: {deg_sig.shape[0]}")

# 上调 / 下调列表（给通路富集用）
up_genes = deg_sig[deg_sig["log2FC"] > 0]["gene"].tolist()
down_genes = deg_sig[deg_sig["log2FC"] < 0]["gene"].tolist()

up_file = os.path.join(deg_dir, f"DEG_up_genes_{DATA_TYPE}_seed{RANDOM_STATE}.txt")
down_file = os.path.join(deg_dir, f"DEG_down_genes_{DATA_TYPE}_seed{RANDOM_STATE}.txt")

with open(up_file, "w") as f:
    for g in up_genes:
        f.write(str(g) + "\n")

with open(down_file, "w") as f:
    for g in down_genes:
        f.write(str(g) + "\n")

print(f"[SAVE] Up-regulated genes   -> {up_file} (n={len(up_genes)})")
print(f"[SAVE] Down-regulated genes -> {down_file} (n={len(down_genes)})")


# =========================
# 5. 通路富集分析（只保留 GO-BP）
# =========================
def run_enrichr(gene_list, label, out_prefix):
    if len(gene_list) == 0:
        print(f"[WARN] No genes in list for {label}, skip enrichment.")
        return

    print(f"[INFO] Running GO enrichment for {label}, #genes={len(gene_list)}")

    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets=['GO_Biological_Process_2021'],   # ✅ 只用 GO-BP
        cutoff=0.05
    )

    res = enr.results
    out_file = os.path.join(pathway_dir, f"{out_prefix}_{DATA_TYPE}_seed{RANDOM_STATE}.csv")
    res.to_csv(out_file, index=False)
    print(f"[SAVE] Enrichment result ({label}) -> {out_file}")


run_enrichr(up_genes,   label="Up-regulated",   out_prefix="Pathway_up")
run_enrichr(down_genes, label="Down-regulated", out_prefix="Pathway_down")

print("[DONE] DEG + GO-BP Pathway pipeline finished.")

