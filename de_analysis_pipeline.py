import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import gseapy as gp


# =========================
# 0. 全局配置
# =========================
DATA_TYPE     = "RNA"       # "smallRNA" 或 "RNA"
RANDOM_STATE  = 42               # 跟聚类脚本中的 seed 一致
TOP_N         = 200              # 导出 p 最小的前 N 个
STRICT_FDR    = 0.05
STRICT_LOG2FC = 1.0
RELAX_FDR     = 0.10
RELAX_LOG2FC  = 0.50

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

MIRNA_TARGET_FILE = os.path.join(
    BASE_DIR, "data", "reference", "mirna_target_human.csv"
)


# =========================
# 工具函数
# =========================
def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR 校正"""
    pvals = np.asarray(pvals, dtype=float)
    n = pvals.size
    order = np.argsort(pvals)
    ranked = pvals[order]

    fdr = ranked * n / (np.arange(1, n + 1))
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr = np.clip(fdr, 0, 1.0)

    out = np.empty_like(fdr)
    out[order] = fdr
    return out


def load_expression_and_clusters(data_type: str, seed: int):
    """读取表达矩阵和聚类结果，返回 X(样本×特征)、cluster_labels、输出目录们"""
    if data_type == "smallRNA":
        expr_file = os.path.join(
            BASE_DIR, "data", "small RNA-seq", "raw",
            "GSE254878_smallRNAs_raw_counts_expression.csv"
        )
        cluster_file = os.path.join(
            BASE_DIR, "data", "small RNA-seq", "clustering",
            f"cluster_results_smallRNA_seed{seed}.csv"
        )
        base_out = os.path.join(BASE_DIR, "data", "small RNA-seq")
    elif data_type == "RNA":
        expr_file = os.path.join(
            BASE_DIR, "data", "RNA-seq", "raw",
            "GSE254877_raw_counts_expression.csv"
        )
        cluster_file = os.path.join(
            BASE_DIR, "data", "RNA-seq", "clustering",
            f"cluster_results_RNA_seed{seed}.csv"
        )
        base_out = os.path.join(BASE_DIR, "data", "RNA-seq")
    else:
        raise ValueError("DATA_TYPE 必须是 'smallRNA' 或 'RNA'")

    print(f"[INFO] Using expression file: {expr_file}")
    print(f"[INFO] Using cluster file   : {cluster_file}")

    # ---- 表达矩阵 ----
    df = pd.read_csv(expr_file, low_memory=False)
    df = df.drop(index=0)                 # 第一行注释
    gene_col = df.columns[0]             # 第一列作为 ID
    df = df.set_index(gene_col)

    if "type" in df.columns:             # smallRNA 的 type 列
        df = df.drop(columns=["type"])

    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all").fillna(0)

    X = df.T  # 行=样本，列=基因/miRNA
    print(f"[INFO] Expression matrix shape: {X.shape}")

    # ---- 聚类标签 ----
    clust = pd.read_csv(cluster_file).set_index("Sample_ID")
    clust = clust.loc[X.index]           # 对齐样本顺序
    labels = clust["Cluster"].values

    print(f"[INFO] Loaded cluster labels, k = {len(np.unique(labels))}")
    print(f"[INFO] Cluster counts:\n{clust['Cluster'].value_counts()}")

    # 输出目录
    deg_dir     = os.path.join(base_out, "deg")
    pathway_dir = os.path.join(base_out, "pathway")
    targets_dir = os.path.join(base_out, "targets")
    os.makedirs(deg_dir, exist_ok=True)
    os.makedirs(pathway_dir, exist_ok=True)
    os.makedirs(targets_dir, exist_ok=True)

    return X, labels, deg_dir, pathway_dir, targets_dir


def compute_deg(X_df: pd.DataFrame,
                labels: np.ndarray,
                group1: int = 0,
                group2: int = 1) -> pd.DataFrame:
    """Cluster group1 vs group2，做差异分析，返回 DEG 表"""
    mask1 = (labels == group1)
    mask2 = (labels == group2)

    X1 = X_df.loc[mask1]
    X2 = X_df.loc[mask2]

    print(f"[INFO] DEG: Cluster {group1} (n={mask1.sum()}) "
          f"vs Cluster {group2} (n={mask2.sum()})")

    p_vals = []
    log2fc = []

    for gene in X_df.columns:
        g1 = X1[gene].values
        g2 = X2[gene].values

        _, p = ttest_ind(g1, g2, equal_var=False)
        p_vals.append(p)

        fc = (g1.mean() + 1e-9) / (g2.mean() + 1e-9)
        log2fc.append(np.log2(fc))

    p_vals = np.array(p_vals)
    fdr = bh_fdr(p_vals)

    deg = pd.DataFrame({
        "gene":   X_df.columns,
        "log2FC": log2fc,
        "p_value": p_vals,
        "FDR":    fdr
    }).sort_values("FDR")

    return deg


def save_deg_tables(deg: pd.DataFrame,
                    deg_dir: str,
                    data_type: str,
                    seed: int):
    """保存各种 DEG 结果，并返回宽松阈值下的 up/down 列表"""
    # 全表
    full_path = os.path.join(deg_dir, f"DEG_full_{data_type}_seed{seed}.csv")
    deg.to_csv(full_path, index=False)
    print(f"[SAVE] Full DEG table -> {full_path}")

    # 全局统计
    print("========== DEG summary (global) ==========")
    print("Total genes:", deg.shape[0])
    print("FDR < 0.05:", (deg["FDR"] < 0.05).sum())
    print("FDR < 0.10:", (deg["FDR"] < 0.10).sum())
    print("|log2FC| > 1.0:", (deg["log2FC"].abs() > 1.0).sum())
    print("|log2FC| > 0.5:", (deg["log2FC"].abs() > 0.5).sum())
    print("===========================================")

    # 严格版
    strict = deg[(deg["FDR"] < STRICT_FDR) &
                 (deg["log2FC"].abs() > STRICT_LOG2FC)]
    strict_path = os.path.join(
        deg_dir,
        f"DEG_sig_strict_{data_type}_FDR{STRICT_FDR}_log2FC{STRICT_LOG2FC}_seed{seed}.csv"
    )
    strict.to_csv(strict_path, index=False)
    print(f"[SAVE] Strict significant DEGs -> {strict_path}")
    print(f"[INFO] #Strict DEGs: {strict.shape[0]}")

    # Top N
    top = deg.sort_values("p_value").head(TOP_N)
    top_path = os.path.join(
        deg_dir, f"DEG_top{TOP_N}_{data_type}_seed{seed}.csv"
    )
    top.to_csv(top_path, index=False)
    print(f"[SAVE] Top {TOP_N} DEG (by p-value) -> {top_path}")

    # 宽松版（给富集、target 映射用）
    print(f"[INFO] Using relaxed cutoff: FDR<{RELAX_FDR}, |log2FC|>{RELAX_LOG2FC}")
    relaxed = deg[(deg["FDR"] < RELAX_FDR) &
                  (deg["log2FC"].abs() > RELAX_LOG2FC)]
    relaxed_path = os.path.join(
        deg_dir,
        f"DEG_sig_relaxed_{data_type}_FDR{RELAX_FDR}_log2FC{RELAX_LOG2FC}_seed{seed}.csv"
    )
    relaxed.to_csv(relaxed_path, index=False)
    print(f"[SAVE] Relaxed significant DEGs -> {relaxed_path}")
    print(f"[INFO] #Relaxed DEGs: {relaxed.shape[0]}")

    up_genes = relaxed[relaxed["log2FC"] > 0]["gene"].tolist()
    down_genes = relaxed[relaxed["log2FC"] < 0]["gene"].tolist()

    up_file = os.path.join(deg_dir, f"DEG_up_genes_{data_type}_seed{seed}.txt")
    down_file = os.path.join(deg_dir, f"DEG_down_genes_{data_type}_seed{seed}.txt")

    with open(up_file, "w") as f:
        f.write("\n".join(map(str, up_genes)))
    with open(down_file, "w") as f:
        f.write("\n".join(map(str, down_genes)))

    print(f"[SAVE] Up-regulated genes   -> {up_file} (n={len(up_genes)})")
    print(f"[SAVE] Down-regulated genes -> {down_file} (n={len(down_genes)})")

    return up_genes, down_genes


def run_enrichr(gene_list, label, out_prefix, pathway_dir, data_type, seed):
    """调用 gseapy.enrichr 做 GO-BP 富集并保存结果"""
    if not gene_list:
        print(f"[WARN] No genes for {label}, skip enrichment.")
        return

    print(f"[INFO] Running GO enrichment for {label}, #genes={len(gene_list)}")

    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets=["GO_Biological_Process_2021"],
        cutoff=0.05
    )

    res = enr.results
    out_file = os.path.join(
        pathway_dir, f"{out_prefix}_{data_type}_seed{seed}.csv"
    )
    res.to_csv(out_file, index=False)
    print(f"[SAVE] Enrichment result ({label}) -> {out_file}")


def smallrna_to_targets(up_mirna, down_mirna, targets_dir):
    """根据 mirna_target_human.csv 做 smallRNA → mRNA 映射"""
    print(f"[INFO] Loading miRNA target mapping from: {MIRNA_TARGET_FILE}")
    if not os.path.exists(MIRNA_TARGET_FILE):
        raise FileNotFoundError(
            f"miRNA–target mapping file not found: {MIRNA_TARGET_FILE}\n"
            "请准备一个包含列: miRNA, TargetGene 的 CSV."
        )

    db = pd.read_csv(MIRNA_TARGET_FILE)
    db["miRNA"] = db["miRNA"].astype(str).str.strip()
    db["TargetGene"] = db["TargetGene"].astype(str).str.strip()

    targets_up = (
        db[db["miRNA"].isin(up_mirna)]["TargetGene"]
        .dropna().unique().tolist()
    )
    targets_down = (
        db[db["miRNA"].isin(down_mirna)]["TargetGene"]
        .dropna().unique().tolist()
    )

    print(f"[INFO] #Up miRNA: {len(up_mirna)}, mapped target genes: {len(targets_up)}")
    print(f"[INFO] #Down miRNA: {len(down_mirna)}, mapped target genes: {len(targets_down)}")

    up_path = os.path.join(targets_dir, "Targets_up_from_smallRNA_seed42.txt")
    down_path = os.path.join(targets_dir, "Targets_down_from_smallRNA_seed42.txt")

    with open(up_path, "w") as f:
        f.write("\n".join(map(str, targets_up)))
    with open(down_path, "w") as f:
        f.write("\n".join(map(str, targets_down)))

    print(f"[SAVE] Up miRNA target genes   -> {up_path}")
    print(f"[SAVE] Down miRNA target genes -> {down_path}")

    return targets_up, targets_down


# =========================
# main 流程
# =========================
def main():
    # 1. 读数据 + 聚类
    X, labels, deg_dir, pathway_dir, targets_dir = load_expression_and_clusters(
        DATA_TYPE, RANDOM_STATE
    )

    # 2. DEG
    deg = compute_deg(X, labels, group1=0, group2=1)

    # 3. 保存 DEG 各种表，拿到 up/down 基因或 miRNA
    up_genes, down_genes = save_deg_tables(
        deg, deg_dir, DATA_TYPE, RANDOM_STATE
    )

    # 4. 富集分析
    if DATA_TYPE == "RNA":
        # 直接对 diff 基因做 GO-BP 富集
        run_enrichr(up_genes,
                    label="Up-regulated genes (RNA, relaxed)",
                    out_prefix="Pathway_up",
                    pathway_dir=pathway_dir,
                    data_type=DATA_TYPE,
                    seed=RANDOM_STATE)

        run_enrichr(down_genes,
                    label="Down-regulated genes (RNA, relaxed)",
                    out_prefix="Pathway_down",
                    pathway_dir=pathway_dir,
                    data_type=DATA_TYPE,
                    seed=RANDOM_STATE)

    elif DATA_TYPE == "smallRNA":
        # 先 miRNA → mRNA，再对靶基因富集
        targets_up, targets_down = smallrna_to_targets(
            up_genes, down_genes, targets_dir
        )

        run_enrichr(targets_up,
                    label="Targets of up-regulated smallRNA (relaxed)",
                    out_prefix="Pathway_targets_up",
                    pathway_dir=pathway_dir,
                    data_type=DATA_TYPE,
                    seed=RANDOM_STATE)

        run_enrichr(targets_down,
                    label="Targets of down-regulated smallRNA (relaxed)",
                    out_prefix="Pathway_targets_down",
                    pathway_dir=pathway_dir,
                    data_type=DATA_TYPE,
                    seed=RANDOM_STATE)

    print("[DONE] DEG + (smallRNA→mRNA) + GO-BP Pathway pipeline finished.")


if __name__ == "__main__":
    main()



