import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import gseapy as gp

# =====================================================
# 0. 全局配置
# =====================================================
RANDOM_STATE  = 42
TOP_N         = 200
STRICT_FDR    = 0.05
STRICT_LOG2FC = 1.0
RELAX_FDR     = 0.10
RELAX_LOG2FC  = 0.50

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# miRNA -> mRNA
MIRNA_TARGET_FILE = os.path.join(
    BASE_DIR, "data", "reference", "mirna_target_human.csv"
)

# lncRNA -> miRNA
LNCRNA_MIRNA_FILE = os.path.join(
    BASE_DIR, "data", "reference", "lncrna_mirna_human.csv"
)

# 集中管理整合结果目录
INTEGRATED_DIR = os.path.join(BASE_DIR, "data", "integrated_results")
OVERLAP_DIR    = os.path.join(INTEGRATED_DIR, "overlap")
TRIPLETS2_DIR  = os.path.join(INTEGRATED_DIR, "triplets")   # miRNA–mRNA 两层
CERNA3_DIR     = os.path.join(INTEGRATED_DIR, "ceRNA3")     # lncRNA–miRNA–mRNA 三层
CERNA3_TRIP    = os.path.join(CERNA3_DIR, "triplets")
CERNA3_ENRICH  = os.path.join(CERNA3_DIR, "enrichment")

os.makedirs(OVERLAP_DIR, exist_ok=True)
os.makedirs(TRIPLETS2_DIR, exist_ok=True)
os.makedirs(CERNA3_TRIP, exist_ok=True)
os.makedirs(CERNA3_ENRICH, exist_ok=True)


# =====================================================
# 工具函数
# =====================================================
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
    """
    读取表达矩阵和聚类结果
    返回:
        X: (样本 × 特征) DataFrame
        labels: 聚类标签
        deg_dir, pathway_dir, targets_dir (targets_dir 仅 smallRNA 用)
    """
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

    elif data_type == "lncRNA":
        # 建议你把 lncRNA 矩阵放在这个位置
        expr_file = os.path.join(
            BASE_DIR, "data", "lncRNA-seq", "raw",
            "GSE254877_lncRNA_raw_counts_expression.csv"
        )
        # 聚类可以直接复用 RNA 的结果（样本一一对应）
        cluster_file = os.path.join(
            BASE_DIR, "data", "RNA-seq", "clustering",
            f"cluster_results_RNA_seed{seed}.csv"
        )
        base_out = os.path.join(BASE_DIR, "data", "lncRNA-seq")

    else:
        raise ValueError("data_type 必须是 'smallRNA'、'RNA' 或 'lncRNA'")

    print(f"\n===== [{data_type}] Load data =====")
    print(f"[INFO] Using expression file: {expr_file}")
    print(f"[INFO] Using cluster file   : {cluster_file}")

    if not os.path.exists(expr_file):
        raise FileNotFoundError(f"表达矩阵不存在: {expr_file}")
    if not os.path.exists(cluster_file):
        raise FileNotFoundError(f"聚类结果不存在: {cluster_file}")

    # ---- 表达矩阵 ----
    df = pd.read_csv(expr_file, low_memory=False)
    df = df.drop(index=0)                 # 第一行注释
    gene_col = df.columns[0]              # 第一列作为 ID
    df = df.set_index(gene_col)

    if "type" in df.columns:              # smallRNA 的 type 列
        df = df.drop(columns=["type"])

    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all").fillna(0)

    X = df.T  # 行=样本，列=基因/miRNA/lncRNA
    print(f"[INFO] Expression matrix shape: {X.shape}")

    # ---- 聚类标签 ----
    clust = pd.read_csv(cluster_file).set_index("Sample_ID")
    clust = clust.loc[X.index]           # 对齐样本
    labels = clust["Cluster"].values

    print(f"[INFO] Loaded cluster labels, k = {len(np.unique(labels))}")
    print(f"[INFO] Cluster counts:\n{clust['Cluster'].value_counts()}")

    deg_dir     = os.path.join(base_out, "deg")
    pathway_dir = os.path.join(base_out, "pathway")
    targets_dir = os.path.join(base_out, "targets")  # smallRNA 用

    os.makedirs(deg_dir, exist_ok=True)
    os.makedirs(pathway_dir, exist_ok=True)
    os.makedirs(targets_dir, exist_ok=True)

    return X, labels, deg_dir, pathway_dir, targets_dir


def compute_deg(X_df: pd.DataFrame,
                labels: np.ndarray,
                group1: int = 0,
                group2: int = 1) -> pd.DataFrame:
    """Cluster group1 vs group2 差异分析"""
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
        "gene":    X_df.columns,
        "log2FC":  log2fc,
        "p_value": p_vals,
        "FDR":     fdr
    }).sort_values("FDR")

    return deg


def save_deg_tables(deg: pd.DataFrame,
                    deg_dir: str,
                    data_type: str,
                    seed: int):
    """
    保存:
      - Full
      - Strict
      - Top N
      - Relaxed (后续富集 & ceRNA 用)
    返回:
      up_genes_relaxed, down_genes_relaxed, relaxed_df
    """
    full_path = os.path.join(deg_dir, f"DEG_full_{data_type}_seed{seed}.csv")
    deg.to_csv(full_path, index=False)
    print(f"[SAVE] Full DEG table -> {full_path}")

    print("---------- DEG summary ----------")
    print("Total genes:", deg.shape[0])
    print("FDR < 0.05:", (deg["FDR"] < 0.05).sum())
    print("FDR < 0.10:", (deg["FDR"] < 0.10).sum())
    print("|log2FC| > 1.0:", (deg["log2FC"].abs() > 1.0).sum())
    print("|log2FC| > 0.5:", (deg["log2FC"].abs() > 0.5).sum())
    print("---------------------------------")

    strict = deg[(deg["FDR"] < STRICT_FDR) &
                 (deg["log2FC"].abs() > STRICT_LOG2FC)]
    strict_path = os.path.join(
        deg_dir,
        f"DEG_sig_strict_{data_type}_FDR{STRICT_FDR}_log2FC{STRICT_LOG2FC}_seed{seed}.csv"
    )
    strict.to_csv(strict_path, index=False)
    print(f"[SAVE] Strict significant DEGs -> {strict_path}")
    print(f"[INFO] #Strict DEGs: {strict.shape[0]}")

    top = deg.sort_values("p_value").head(TOP_N)
    top_path = os.path.join(
        deg_dir, f"DEG_top{TOP_N}_{data_type}_seed{seed}.csv"
    )
    top.to_csv(top_path, index=False)
    print(f"[SAVE] Top {TOP_N} DEG (by p-value) -> {top_path}")

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

    return up_genes, down_genes, relaxed


def run_enrichr(gene_list, label, out_prefix, pathway_dir, data_type, seed):
    """GO_Biological_Process_2021 富集"""
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


def smallrna_to_targets(up_mirna, down_mirna, targets_dir, seed):
    """miRNA -> mRNA 靶基因映射"""
    print(f"[INFO] Loading miRNA target mapping from: {MIRNA_TARGET_FILE}")
    if not os.path.exists(MIRNA_TARGET_FILE):
        raise FileNotFoundError(
            f"miRNA–target mapping file not found: {MIRNA_TARGET_FILE}"
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

    up_path = os.path.join(targets_dir, f"Targets_up_from_smallRNA_seed{seed}.txt")
    down_path = os.path.join(targets_dir, f"Targets_down_from_smallRNA_seed{seed}.txt")

    with open(up_path, "w") as f:
        f.write("\n".join(map(str, targets_up)))
    with open(down_path, "w") as f:
        f.write("\n".join(map(str, targets_down)))

    print(f"[SAVE] Up miRNA target genes   -> {up_path}")
    print(f"[SAVE] Down miRNA target genes -> {down_path}")

    return targets_up, targets_down, db


def build_smallrna_rna_overlap(rna_up, rna_down,
                               targets_up, targets_down,
                               small_up, small_down,
                               mirna_db,
                               seed: int):
    """两层：smallRNA–mRNA 的 overlap + triplets（你之前已经用过）"""
    rna_up_set    = set(rna_up)
    rna_down_set  = set(rna_down)
    t_up_set      = set(targets_up)
    t_down_set    = set(targets_down)

    overlap_up_up      = sorted(t_up_set   & rna_up_set)
    overlap_up_down    = sorted(t_up_set   & rna_down_set)
    overlap_down_up    = sorted(t_down_set & rna_up_set)
    overlap_down_down  = sorted(t_down_set & rna_down_set)

    def save_list(lst, name):
        path = os.path.join(OVERLAP_DIR, name)
        with open(path, "w") as f:
            f.write("\n".join(lst))
        print(f"[SAVE] Overlap list -> {path} (n={len(lst)})")

    save_list(overlap_up_up,     f"overlap_targetsUp_RNAup_seed{seed}.txt")
    save_list(overlap_up_down,   f"overlap_targetsUp_RNAdown_seed{seed}.txt")
    save_list(overlap_down_up,   f"overlap_targetsDown_RNAup_seed{seed}.txt")
    save_list(overlap_down_down, f"overlap_targetsDown_RNAdown_seed{seed}.txt")

    # 构建两层 triplets（miRNA–mRNA）
    all_rna_deg = set(rna_up) | set(rna_down)
    all_small   = set(small_up) | set(small_down)

    trip = mirna_db[
        mirna_db["miRNA"].isin(all_small) &
        mirna_db["TargetGene"].isin(all_rna_deg)
    ].copy()

    small_reg = {g: "up" for g in small_up}
    small_reg.update({g: "down" for g in small_down})

    rna_reg = {g: "up" for g in rna_up}
    rna_reg.update({g: "down" for g in rna_down})

    trip["miRNA_regulation"] = trip["miRNA"].map(small_reg).fillna("NA")
    trip["RNA_regulation"]   = trip["TargetGene"].map(rna_reg).fillna("NA")
    trip = trip[trip["RNA_regulation"] != "NA"]

    trip_out = os.path.join(
        TRIPLETS2_DIR, f"miRNA_target_RNA_triplets_seed{seed}.csv"
    )
    trip.to_csv(trip_out, index=False)
    print(f"[SAVE] miRNA–target–RNA triplets -> {trip_out} (n={trip.shape[0]})")


def build_ceRNA3_triplets(lnc_up, lnc_down,
                          mir_up, mir_down,
                          rna_up, rna_down,
                          mirna_db, lncrna_db,
                          seed: int):
    """
    三层：lncRNA–miRNA–mRNA ceRNA 网络
    规则：
      Pattern1: lnc↑, miR↓, mRNA↑
      Pattern2: lnc↓, miR↑, mRNA↓
    """
    print("\n===== Build 3-layer ceRNA triplets (lncRNA–miRNA–mRNA) =====")

    # clean lncRNA–miRNA mapping
    lncrna_db["lncRNA"] = lncrna_db["lncRNA"].astype(str).str.strip()
    lncrna_db["miRNA"]  = lncrna_db["miRNA"].astype(str).str.strip()


    # merge on miRNA
    trip3 = lncrna_db.merge(mirna_db, on="miRNA", how="inner")
    trip3.rename(columns={"TargetGene": "mRNA"}, inplace=True)

    # 构建上下调字典
    lnc_reg = {g: "up" for g in lnc_up}
    lnc_reg.update({g: "down" for g in lnc_down})

    mir_reg = {g: "up" for g in mir_up}
    mir_reg.update({g: "down" for g in mir_down})

    rna_reg = {g: "up" for g in rna_up}
    rna_reg.update({g: "down" for g in rna_down})

    trip3["lncRNA_reg"] = trip3["lncRNA"].map(lnc_reg).fillna("NA")
    trip3["miRNA_reg"]  = trip3["miRNA"].map(mir_reg).fillna("NA")
    trip3["mRNA_reg"]   = trip3["mRNA"].map(rna_reg).fillna("NA")

    # 只保留三者都有方向信息的
    trip3 = trip3[
        (trip3["lncRNA_reg"].isin(["up", "down"])) &
        (trip3["miRNA_reg"].isin(["up", "down"])) &
        (trip3["mRNA_reg"].isin(["up", "down"]))
    ].copy()

    all_path = os.path.join(
        CERNA3_TRIP, f"lncRNA_miRNA_mRNA_triplets_all_seed{seed}.csv"
    )
    trip3.to_csv(all_path, index=False)
    print(f"[SAVE] All ceRNA triplets -> {all_path} (n={trip3.shape[0]})")

    # Pattern1: lnc↑, miR↓, mRNA↑
    p1 = trip3[
        (trip3["lncRNA_reg"] == "up") &
        (trip3["miRNA_reg"] == "down") &
        (trip3["mRNA_reg"] == "up")
    ].copy()

    # Pattern2: lnc↓, miR↑, mRNA↓
    p2 = trip3[
        (trip3["lncRNA_reg"] == "down") &
        (trip3["miRNA_reg"] == "up") &
        (trip3["mRNA_reg"] == "down")
    ].copy()

    p1_path = os.path.join(
        CERNA3_TRIP, f"ceRNA_pattern1_lncUp_miRDown_mRNAUp_seed{seed}.csv"
    )
    p2_path = os.path.join(
        CERNA3_TRIP, f"ceRNA_pattern2_lncDown_miRUp_mRNADown_seed{seed}.csv"
    )

    p1.to_csv(p1_path, index=False)
    p2.to_csv(p2_path, index=False)

    print(f"[SAVE] Pattern1 triplets -> {p1_path} (n={p1.shape[0]})")
    print(f"[SAVE] Pattern2 triplets -> {p2_path} (n={p2.shape[0]})")

    # 对每种模式的 mRNA 做 GO 富集（看对应功能）
    genes_p1 = sorted(p1["mRNA"].unique())
    genes_p2 = sorted(p2["mRNA"].unique())

    def enrich_for_ceRNA(genes, name):
        if not genes:
            print(f"[WARN] No genes for {name}, skip enrichment.")
            return
        print(f"[INFO] Enrichr for {name}, #genes={len(genes)}")
        enr = gp.enrichr(
            gene_list=genes,
            gene_sets=["GO_Biological_Process_2021"],
            cutoff=0.05
        )
        res = enr.results
        out_file = os.path.join(
            CERNA3_ENRICH, f"{name}_GO_BP_seed{seed}.csv"
        )
        res.to_csv(out_file, index=False)
        print(f"[SAVE] GO-BP result -> {out_file}")

    enrich_for_ceRNA(genes_p1, "ceRNA_pattern1_mRNA")
    enrich_for_ceRNA(genes_p2, "ceRNA_pattern2_mRNA")


# =====================================================
# main
# =====================================================
def main():
    # ---------- RNA ----------
    X_rna, lab_rna, deg_dir_rna, path_dir_rna, _ = load_expression_and_clusters(
        "RNA", RANDOM_STATE
    )
    deg_rna = compute_deg(X_rna, lab_rna, group1=0, group2=1)
    rna_up, rna_down, _ = save_deg_tables(
        deg_rna, deg_dir_rna, "RNA", RANDOM_STATE
    )

    run_enrichr(rna_up,   "Up-regulated genes (RNA, relaxed)",
                "Pathway_up",   path_dir_rna, "RNA", RANDOM_STATE)
    run_enrichr(rna_down, "Down-regulated genes (RNA, relaxed)",
                "Pathway_down", path_dir_rna, "RNA", RANDOM_STATE)

    # ---------- smallRNA ----------
    X_s, lab_s, deg_dir_s, path_dir_s, targets_dir_s = load_expression_and_clusters(
        "smallRNA", RANDOM_STATE
    )
    deg_s = compute_deg(X_s, lab_s, group1=0, group2=1)
    mir_up, mir_down, _ = save_deg_tables(
        deg_s, deg_dir_s, "smallRNA", RANDOM_STATE
    )

    targets_up, targets_down, mirna_db = smallrna_to_targets(
        mir_up, mir_down, targets_dir_s, RANDOM_STATE
    )

    run_enrichr(targets_up,
                "Targets of up-regulated smallRNA (relaxed)",
                "Pathway_targets_up", path_dir_s, "smallRNA", RANDOM_STATE)
    run_enrichr(targets_down,
                "Targets of down-regulated smallRNA (relaxed)",
                "Pathway_targets_down", path_dir_s, "smallRNA", RANDOM_STATE)

    # 两层 smallRNA–mRNA 的 overlap + triplets（沿用之前逻辑）
    build_smallrna_rna_overlap(
        rna_up=rna_up,
        rna_down=rna_down,
        targets_up=targets_up,
        targets_down=targets_down,
        small_up=mir_up,
        small_down=mir_down,
        mirna_db=mirna_db,
        seed=RANDOM_STATE
    )

    # ---------- lncRNA ----------
    X_lnc, lab_lnc, deg_dir_lnc, path_dir_lnc, _ = load_expression_and_clusters(
        "lncRNA", RANDOM_STATE
    )
    deg_lnc = compute_deg(X_lnc, lab_lnc, group1=0, group2=1)
    lnc_up, lnc_down, _ = save_deg_tables(
        deg_lnc, deg_dir_lnc, "lncRNA", RANDOM_STATE
    )

    # （可选）对 lncRNA 自己做 GO 富集
    run_enrichr(lnc_up,
                "Up-regulated lncRNA (relaxed)",
                "Pathway_up", path_dir_lnc, "lncRNA", RANDOM_STATE)
    run_enrichr(lnc_down,
                "Down-regulated lncRNA (relaxed)",
                "Pathway_down", path_dir_lnc, "lncRNA", RANDOM_STATE)

    # ---------- 读取 lncRNA–miRNA 数据，构建三层 ceRNA ----------
    print(f"[INFO] Loading lncRNA–miRNA mapping from: {LNCRNA_MIRNA_FILE}")
    if not os.path.exists(LNCRNA_MIRNA_FILE):
        raise FileNotFoundError(
            f"lncRNA–miRNA mapping file not found: {LNCRNA_MIRNA_FILE}\n"
            "需要一个包含列 lncRNA, miRNA 的 CSV."
        )
    lncrna_db = pd.read_csv(LNCRNA_MIRNA_FILE)

    build_ceRNA3_triplets(
        lnc_up=lnc_up,
        lnc_down=lnc_down,
        mir_up=mir_up,
        mir_down=mir_down,
        rna_up=rna_up,
        rna_down=rna_down,
        mirna_db=mirna_db,
        lncrna_db=lncrna_db,
        seed=RANDOM_STATE
    )

    print("\n[DONE] RNA + smallRNA + lncRNA + smallRNA→mRNA + lncRNA→miRNA "
          "+ 2-layer & 3-layer ceRNA + GO-BP 全流程完成")


if __name__ == "__main__":
    main()




