import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import gseapy as gp

# =====================================================
# 0. 全局配置
# =====================================================
RANDOM_STATE = 42               # 和聚类脚本里的 seed 一致
TOP_N        = 200              # 导出 p 最小的前 N 个
STRICT_FDR   = 0.05
STRICT_LOG2FC = 1.0
RELAX_FDR     = 0.10
RELAX_LOG2FC  = 0.50

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# smallRNA -> mRNA 映射表
MIRNA_TARGET_FILE = os.path.join(
    BASE_DIR, "data", "reference", "mirna_target_human.csv"
)

# lncRNA -> miRNA 映射表（NPInter 构建）
LNCRNA_MIRNA_FILE = os.path.join(
    BASE_DIR, "data", "reference", "lncrna_mirna_human.csv"
)


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
        labels: 聚类标签 (numpy array)
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
        expr_file = os.path.join(
            BASE_DIR, "data", "lncRNA-seq", "raw",
            "GSE254877_lncRNA_raw_counts_expression.csv"
        )
        # lncRNA 用的是 RNA 的聚类标签（同一批样本）
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

    # 行=样本，列=基因/miRNA
    X = df.T
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
    targets_dir = os.path.join(base_out, "targets")  # smallRNA / lncRNA 用
    os.makedirs(deg_dir, exist_ok=True)
    os.makedirs(pathway_dir, exist_ok=True)
    os.makedirs(targets_dir, exist_ok=True)

    return X, labels, deg_dir, pathway_dir, targets_dir


def compute_deg(X_df: pd.DataFrame,
                labels: np.ndarray,
                group1: int = 0,
                group2: int = 1) -> pd.DataFrame:
    """
    Cluster group1 vs group2，做差异分析，返回 DEG 表
    log2FC = mean(group1) / mean(group2)
    """
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
    保存各种 DEG 结果：
      - Full
      - Strict
      - Top N
      - Relaxed (富集 & overlap 用)
    返回：
      up_genes_relaxed, down_genes_relaxed, relaxed_deg_df
    """
    # 全表
    full_path = os.path.join(deg_dir, f"DEG_full_{data_type}_seed{seed}.csv")
    deg.to_csv(full_path, index=False)
    print(f"[SAVE] Full DEG table -> {full_path}")

    # 全局统计
    print("---------- DEG summary ----------")
    print("Total genes:", deg.shape[0])
    print("FDR < 0.05:", (deg["FDR"] < 0.05).sum())
    print("FDR < 0.10:", (deg["FDR"] < 0.10).sum())
    print("|log2FC| > 1.0:", (deg["log2FC"].abs() > 1.0).sum())
    print("|log2FC| > 0.5:", (deg["log2FC"].abs() > 0.5).sum())
    print("---------------------------------")

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

    # 宽松版（给富集 & target 映射 & overlap / ceRNA 用）
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


def smallrna_to_targets(up_mirna, down_mirna, targets_dir, seed):
    """
    根据 mirna_target_human.csv 做 smallRNA → mRNA 映射
    返回:
        targets_up, targets_down, mirna_db
    """
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
    """
    smallRNA–target–RNA overlap 整合分析

    输出到:
      data/integrated_results/
        - overlap/overlap_targetsUp_RNAup_seedXX.txt / ...
        - triplets/miRNA_target_RNA_triplets_seedXX.csv
    """
    integrated_dir = os.path.join(BASE_DIR, "data", "integrated_results")
    overlap_dir    = os.path.join(integrated_dir, "overlap")
    triplet_dir    = os.path.join(integrated_dir, "triplets")
    os.makedirs(overlap_dir, exist_ok=True)
    os.makedirs(triplet_dir, exist_ok=True)

    rna_up_set    = set(rna_up)
    rna_down_set  = set(rna_down)
    t_up_set      = set(targets_up)
    t_down_set    = set(targets_down)

    # 不同组合的基因 overlap
    overlap_up_up      = sorted(t_up_set   & rna_up_set)
    overlap_up_down    = sorted(t_up_set   & rna_down_set)
    overlap_down_up    = sorted(t_down_set & rna_up_set)
    overlap_down_down  = sorted(t_down_set & rna_down_set)

    def save_list(lst, name):
        path = os.path.join(overlap_dir, name)
        with open(path, "w") as f:
            f.write("\n".join(lst))
        print(f"[SAVE] Overlap list -> {path} (n={len(lst)})")

    # 保存 4 类 overlap
    save_list(overlap_up_up,     f"overlap_targetsUp_RNAup_seed{seed}.txt")
    save_list(overlap_up_down,   f"overlap_targetsUp_RNAdown_seed{seed}.txt")
    save_list(overlap_down_up,   f"overlap_targetsDown_RNAup_seed{seed}.txt")
    save_list(overlap_down_down, f"overlap_targetsDown_RNAdown_seed{seed}.txt")

    # ===== miRNA–target–RNA triplet 表 =====
    # 只保留：miRNA 是差异 smallRNA，TargetGene 是差异 RNA（up 或 down）
    all_rna_deg = set(rna_up) | set(rna_down)
    all_small   = set(small_up) | set(small_down)

    trip = mirna_db[
        mirna_db["miRNA"].isin(all_small) &
        mirna_db["TargetGene"].isin(all_rna_deg)
    ].copy()

    # 加上 smallRNA / RNA 的上下调方向
    small_reg = {}
    for g in small_up:
        small_reg[g] = "up"
    for g in small_down:
        small_reg[g] = "down"

    rna_reg = {}
    for g in rna_up:
        rna_reg[g] = "up"
    for g in rna_down:
        rna_reg[g] = "down"

    trip["miRNA_regulation"] = trip["miRNA"].map(small_reg).fillna("NA")
    trip["RNA_regulation"]   = trip["TargetGene"].map(rna_reg).fillna("NA")

    # 只保留有方向信息的行（真正参与交集的）
    trip = trip[trip["RNA_regulation"] != "NA"]

    trip_out = os.path.join(
        triplet_dir, f"miRNA_target_RNA_triplets_seed{seed}.csv"
    )
    trip.to_csv(trip_out, index=False)
    print(f"[SAVE] miRNA–target–RNA triplets -> {trip_out} (n={trip.shape[0]})")


# =====================================================
# lncRNA -> miRNA 映射 + ceRNA（三层）相关
# =====================================================
def lncrna_to_mirna(lnc_up, lnc_down, targets_dir, seed):
    """
    根据 lncrna_mirna_human.csv 做 lncRNA → miRNA 映射
    返回:
        mir_up, mir_down, lncrna_db
    """
    print(f"[INFO] Loading lncRNA–miRNA mapping from: {LNCRNA_MIRNA_FILE}")
    if not os.path.exists(LNCRNA_MIRNA_FILE):
        raise FileNotFoundError(
            f"lncRNA–miRNA mapping file not found: {LNCRNA_MIRNA_FILE}\n"
            "需要一个包含列 lncRNA, miRNA 的 CSV."
        )

    db = pd.read_csv(LNCRNA_MIRNA_FILE)
    db["lncRNA"] = db["lncRNA"].astype(str).str.strip()
    db["miRNA"]  = db["miRNA"].astype(str).str.strip()

    mir_up = (
        db[db["lncRNA"].isin(lnc_up)]["miRNA"]
        .dropna().unique().tolist()
    )
    mir_down = (
        db[db["lncRNA"].isin(lnc_down)]["miRNA"]
        .dropna().unique().tolist()
    )

    print(f"[INFO] #Up lncRNA: {len(lnc_up)}, mapped miRNAs: {len(mir_up)}")
    print(f"[INFO] #Down lncRNA: {len(lnc_down)}, mapped miRNAs: {len(mir_down)}")

    up_path = os.path.join(targets_dir, f"miRNA_from_up_lncRNA_seed{seed}.txt")
    down_path = os.path.join(targets_dir, f"miRNA_from_down_lncRNA_seed{seed}.txt")

    with open(up_path, "w") as f:
        f.write("\n".join(map(str, mir_up)))
    with open(down_path, "w") as f:
        f.write("\n".join(map(str, mir_down)))

    print(f"[SAVE] miRNA from up-lncRNA   -> {up_path}")
    print(f"[SAVE] miRNA from down-lncRNA -> {down_path}")

    return mir_up, mir_down, db


def build_ceRNA3_triplets(lnc_up, lnc_down,
                          mirna_up, mirna_down,
                          mirna_rna_triplets_file: str,
                          lncrna_mirna_db: pd.DataFrame,
                          seed: int):
    """
    构建三层 ceRNA 关系:
        lncRNA – miRNA – mRNA

    利用:
      - lncRNA–miRNA 映射 (lncrna_mirna_db)
      - miRNA–mRNA triplets (miRNA_target_RNA_triplets_seedXX.csv)
    """
    integrated_dir = os.path.join(BASE_DIR, "data", "integrated_results")
    cerna_dir = os.path.join(integrated_dir, "ceRNA3")
    triplet_dir = os.path.join(cerna_dir, "triplets")
    enrich_dir = os.path.join(cerna_dir, "enrichment")
    os.makedirs(triplet_dir, exist_ok=True)
    os.makedirs(enrich_dir, exist_ok=True)

    print("\n===== Build 3-layer ceRNA triplets (lncRNA–miRNA–mRNA) =====")
    print(f"[INFO] Loading miRNA–mRNA triplets from: {mirna_rna_triplets_file}")

    if not os.path.exists(mirna_rna_triplets_file):
        print("[WARN] miRNA–mRNA triplets file 不存在，跳过 ceRNA3 构建.")
        return

    mirna_rna = pd.read_csv(mirna_rna_triplets_file)

    # 只保留方向明确的 miRNA–mRNA 关系
    mirna_rna = mirna_rna[
        (mirna_rna["miRNA_regulation"].isin(["up", "down"])) &
        (mirna_rna["RNA_regulation"].isin(["up", "down"]))
    ].copy()

    # lncRNA–miRNA 映射里只取在差异 lncRNA & 差异 miRNA 里的
    diff_lnc = set(lnc_up) | set(lnc_down)
    diff_mir = set(mirna_up) | set(mirna_down)

    lm = lncrna_mirna_db[
        lncrna_mirna_db["lncRNA"].isin(diff_lnc) &
        lncrna_mirna_db["miRNA"].isin(diff_mir)
    ].copy()

    print(f"[INFO] ceRNA3: #lncRNA–miRNA pairs after diff filter: {lm.shape[0]}")
    print(f"[INFO] ceRNA3: #miRNA–mRNA edges (diff)             : {mirna_rna.shape[0]}")

    # 给 lncRNA 加方向
    lnc_reg = {}
    for g in lnc_up:
        lnc_reg[g] = "up"
    for g in lnc_down:
        lnc_reg[g] = "down"

    lm["lncRNA_regulation"] = lm["lncRNA"].map(lnc_reg).fillna("NA")

    # miRNA 方向信息已有: mirna_up / mirna_down
    mir_reg = {}
    for g in mirna_up:
        mir_reg[g] = "up"
    for g in mirna_down:
        mir_reg[g] = "down"

    lm["miRNA_regulation"] = lm["miRNA"].map(mir_reg).fillna("NA")

    # 和 miRNA–mRNA 的表做 merge -> 得到三层 triplet
    trip = lm.merge(
        mirna_rna,
        left_on="miRNA",
        right_on="miRNA",
        suffixes=("_lnc", "_mRNA")
    )

    # 合并后会有 miRNA_regulation_lnc 和 miRNA_regulation_mRNA
    # 我们用来自 miRNA–mRNA 表的那个（_mRNA），并统一重命名回 miRNA_regulation
    if "miRNA_regulation_mRNA" in trip.columns:
        trip = trip.rename(columns={"miRNA_regulation_mRNA": "miRNA_regulation"})
    elif "miRNA_regulation" not in trip.columns:
        # 如果两个名字都没有，说明上游数据结构变了
        raise KeyError(
            "No miRNA_regulation column found after merge; "
            "check mirna_rna columns."
        )

    # 保留需要的列
    trip = trip[[
        "lncRNA", "miRNA", "TargetGene",
        "lncRNA_regulation",
        "miRNA_regulation",
        "RNA_regulation"
    ]].copy()


    out_file = os.path.join(
        triplet_dir,
        f"lncRNA_miRNA_mRNA_triplets_all_seed{seed}.csv"
    )
    trip.to_csv(out_file, index=False)
    print(f"[SAVE] 3-layer ceRNA triplets -> {out_file} (n={trip.shape[0]})")

    # 你可以在这里继续做 pattern 过滤 / 富集分析（比如 lncUp–miRDown–mRNAUp）
    # 这里先简单做两个常见 pattern 的子集输出:
    pattern1 = trip[
        (trip["lncRNA_regulation"] == "up") &
        (trip["miRNA_regulation"] == "down") &
        (trip["RNA_regulation"] == "up")
    ]
    pattern2 = trip[
        (trip["lncRNA_regulation"] == "down") &
        (trip["miRNA_regulation"] == "up") &
        (trip["RNA_regulation"] == "down")
    ]

    p1_file = os.path.join(
        triplet_dir,
        f"ceRNA_pattern1_lncUp_miRDown_mRNAUp_seed{seed}.csv"
    )
    p2_file = os.path.join(
        triplet_dir,
        f"ceRNA_pattern2_lncDown_miRUp_mRNADown_seed{seed}.csv"
    )
    pattern1.to_csv(p1_file, index=False)
    pattern2.to_csv(p2_file, index=False)

    print(f"[SAVE] ceRNA pattern1 (lncUp–miRDown–mRNAUp)   -> {p1_file} (n={pattern1.shape[0]})")
    print(f"[SAVE] ceRNA pattern2 (lncDown–miRUp–mRNADown) -> {p2_file} (n={pattern2.shape[0]})")


def postprocess_mirna_mrna_triplets(seed: int = RANDOM_STATE):
    """对 miRNA–mRNA triplets 做方向过滤 + GO-BP 富集

    读取:
        data/integrated_results/triplets/miRNA_target_RNA_triplets_seedXX.csv

    输出:
        data/integrated_results/triplets/triplets_miRup_RNAdown_seedXX.csv
        data/integrated_results/triplets/triplets_miRdown_RNAup_seedXX.csv
        data/integrated_results/enrichment/GO_miRup_RNAdown_targets_GO_BP_seedXX.csv
        data/integrated_results/enrichment/GO_miRdown_RNAup_targets_GO_BP_seedXX.csv
    """
    integrated_dir = os.path.join(BASE_DIR, "data", "integrated_results")
    triplet_dir = os.path.join(integrated_dir, "triplets")
    enrich_dir = os.path.join(integrated_dir, "enrichment")
    os.makedirs(triplet_dir, exist_ok=True)
    os.makedirs(enrich_dir, exist_ok=True)

    triplet_file = os.path.join(
        triplet_dir,
        f"miRNA_target_RNA_triplets_seed{seed}.csv"
    )

    if not os.path.exists(triplet_file):
        print(f"[WARN] Triplet file not found: {triplet_file}, skip post-processing.")
        return

    print("\n===== [Postprocess miRNA–mRNA triplets] =====")
    print(f"[INFO] Loading triplets: {triplet_file}")
    df = pd.read_csv(triplet_file)

    # 只保留真正参与调控的（两边都有方向信息）
    df = df[
        (df["miRNA_regulation"].isin(["up", "down"])) &
        (df["RNA_regulation"].isin(["up", "down"]))
    ].copy()

    print(f"[INFO] Total triplets with direction info: {df.shape[0]}")

    # 1) miRNA 上调，mRNA 下调（典型抑制型）
    up_miR_down_RNA = df[
        (df["miRNA_regulation"] == "up") &
        (df["RNA_regulation"] == "down")
    ].copy()

    # 2) miRNA 下调，mRNA 上调（释放抑制）
    down_miR_up_RNA = df[
        (df["miRNA_regulation"] == "down") &
        (df["RNA_regulation"] == "up")
    ].copy()

    # 保存详细表到 triplets/
    up_down_file = os.path.join(
        triplet_dir,
        f"triplets_miRup_RNAdown_seed{seed}.csv"
    )
    down_up_file = os.path.join(
        triplet_dir,
        f"triplets_miRdown_RNAup_seed{seed}.csv"
    )

    up_miR_down_RNA.to_csv(up_down_file, index=False)
    down_miR_up_RNA.to_csv(down_up_file, index=False)

    print(f"[SAVE] miR↑–RNA↓ triplets -> {up_down_file} (n={up_miR_down_RNA.shape[0]})")
    print(f"[SAVE] miR↓–RNA↑ triplets -> {down_up_file} (n={down_miR_up_RNA.shape[0]})")

    # 针对各自的 mRNA 做 GO-BP 富集，存到 enrichment/
    genes_up_miR_down_RNA = sorted(up_miR_down_RNA["TargetGene"].dropna().unique())
    genes_down_miR_up_RNA = sorted(down_miR_up_RNA["TargetGene"].dropna().unique())

    def _run_enrichr_for_triplets(gene_list, label, out_prefix):
        if not gene_list:
            print(f"[WARN] No genes for {label}, skip enrichment.")
            return
        print(f"[INFO] Enrichr for {label}, #genes = {len(gene_list)}")
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=["GO_Biological_Process_2021"],
            cutoff=0.05,
        )
        res = enr.results
        out_path = os.path.join(enrich_dir, f"{out_prefix}_GO_BP_seed{seed}.csv")
        res.to_csv(out_path, index=False)
        print(f"[SAVE] GO-BP result -> {out_path}")

    _run_enrichr_for_triplets(
        genes_up_miR_down_RNA,
        label="miR_up_RNA_down_targets",
        out_prefix="GO_miRup_RNAdown_targets",
    )
    _run_enrichr_for_triplets(
        genes_down_miR_up_RNA,
        label="miR_down_RNA_up_targets",
        out_prefix="GO_miRdown_RNAup_targets",
    )