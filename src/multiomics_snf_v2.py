# src/multiomics_snf_v2.py
# -*- coding: utf-8 -*-

"""
multiomics_snf_v2.py

多组学聚类 V2.0 工具包（真正 SNF + early-fusion baseline 对比 + UMAP 可视化 + RNA DEG/Pathway）：

提供一个入口函数：
    run_snf_v2_pipeline()

由外部脚本调用，例如：
    from src.multiomics_snf_v2 import run_snf_v2_pipeline
    run_snf_v2_pipeline()
"""

import os
import re
from collections import Counter

import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import pairwise_distances

# UMAP 如未安装则自动退回 PCA 2D
try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False

from src.methods_utils import BASE_DIR, RANDOM_STATE


# ============================================================
# 0. 工具函数：样本 ID 简化
# ============================================================
def simplify_sample_id(s: str):
    """
    把原始样本名简化成病人 ID，比如：
    RNA:    NS.1402.003.NEBNext_dual_i7_177---NEBNext_dual_i5_177.2066_R1.bam -> '2066'
    miRNA:  NS.1404.001.NEBNext_S02.2066                                     -> '2066'
    """
    s = str(s)
    if "Unnamed" in s:
        return None  # 明显是奇怪列，直接丢掉

    nums4 = re.findall(r"(\d{4})", s)
    if nums4:
        return nums4[-1]
    return s.strip()


# ============================================================
# 1. 读取表达矩阵
# ============================================================
def load_expression_only(data_type: str) -> pd.DataFrame:
    """
    加载表达矩阵（样本 × 特征）。
    """
    if data_type == "RNA":
        expr_file = os.path.join(
            BASE_DIR, "data", "RNA-seq", "raw",
            "GSE254877_raw_counts_expression.csv"
        )
    elif data_type == "smallRNA":
        expr_file = os.path.join(
            BASE_DIR, "data", "small RNA-seq", "raw",
            "GSE254878_smallRNAs_raw_counts_expression.csv"
        )
    else:
        raise ValueError("data_type 必须是 'RNA' 或 'smallRNA'")

    print(f"\n[LOAD] {data_type} expression from: {expr_file}")
    if not os.path.exists(expr_file):
        raise FileNotFoundError(f"表达矩阵不存在: {expr_file}")

    df = pd.read_csv(expr_file, low_memory=False)
    df = df.drop(index=0)
    gene_col = df.columns[0]
    df = df.set_index(gene_col)

    if "type" in df.columns:
        df = df.drop(columns=["type"])

    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all").fillna(0)

    X = df.T
    X.index = X.index.astype(str).str.strip()

    print(f"[INFO] {data_type} matrix shape: {X.shape} (samples × features)")
    print(f"[DEBUG] {data_type} sample IDs (head): {list(X.index[:5])}")
    return X


# ============================================================
# 2. 高变特征
# ============================================================
def select_top_variable_features(X: pd.DataFrame, top_n: int = 2000) -> pd.DataFrame:
    if X.shape[1] <= top_n:
        return X

    var = X.var(axis=0)
    keep = var.sort_values(ascending=False).head(top_n).index
    X_sub = X[keep]
    print(f"[INFO] Selected top {top_n} variable features: {X_sub.shape}")
    return X_sub


# ============================================================
# 3. PCA 降维
# ============================================================
def pca_transform(X: pd.DataFrame, n_components: int = 50) -> np.ndarray:
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X.values)

    max_dim = min(n_components, X_scaled.shape[0] - 1)
    if max_dim < 2:
        raise ValueError("样本太少，PCA 维度 < 2，无法聚类。")

    pca = PCA(n_components=max_dim, random_state=RANDOM_STATE)
    X_pca = pca.fit_transform(X_scaled)

    print(f"[INFO] PCA -> {X_pca.shape[1]} dims, explained_var_ratio_sum="
          f"{pca.explained_variance_ratio_.sum():.3f}")
    return X_pca


# ============================================================
# 4. SNF: affinity + fusion + clustering
# ============================================================
def construct_affinity_matrix(X: np.ndarray, K: int = 20, mu: float = 0.5) -> np.ndarray:
    print("[SNF] Constructing affinity matrix...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    dist = pairwise_distances(X_scaled, metric="euclidean")
    N = dist.shape[0]

    dist_sort = np.sort(dist, axis=1)
    K_eff = min(K, N - 1)
    sigma = np.mean(dist_sort[:, 1:K_eff+1], axis=1) + 1e-8

    sigma_i = sigma.reshape(-1, 1)
    sigma_j = sigma.reshape(1, -1)
    W = np.exp(- (dist ** 2) / (sigma_i * sigma_j))

    W_knn = np.zeros_like(W)
    for i in range(N):
        idx = np.argsort(dist[i, :])[1:K_eff+1]
        W_knn[i, idx] = W[i, idx]

    W_sym = (W_knn + W_knn.T) / 2.0

    row_sum = W_sym.sum(axis=1, keepdims=True)
    row_sum[row_sum == 0] = 1.0
    W_norm = W_sym / row_sum
    print("[SNF] Affinity matrix constructed.")
    return W_norm


def snf(affinity_list, K: int = 20, t: int = 20) -> np.ndarray:
    print(f"[SNF] Running SNF fusion with {len(affinity_list)} networks, "
          f"K={K}, t={t}...")

    P_list = [W.copy() for W in affinity_list]
    M = len(P_list)
    N = P_list[0].shape[0]

    knn_masks = []
    for P in P_list:
        mask = np.zeros_like(P, dtype=bool)
        for i in range(N):
            idx = np.argsort(P[i, :])[::-1][:K]
            mask[i, idx] = True
        knn_masks.append(mask)

    for it in range(t):
        new_P_list = []
        for m in range(M):
            others = [P_list[j] for j in range(M) if j != m]
            P_bar = sum(others) / (M - 1)

            Pm = P_list[m]
            P_new = Pm @ P_bar @ Pm.T

            P_new[~knn_masks[m]] = 0.0

            row_sum = P_new.sum(axis=1, keepdims=True)
            row_sum[row_sum == 0] = 1.0
            P_new = P_new / row_sum

            new_P_list.append(P_new)

        P_list = new_P_list
        if (it + 1) % 5 == 0:
            print(f"[SNF] Iteration {it+1}/{t} done")

    W_fused = sum(P_list) / M
    print("[SNF] Fusion finished.")
    return W_fused


def snf_clustering(W_fused: np.ndarray,
                   candidate_K=(3, 4, 5),
                   random_state: int = RANDOM_STATE):
    print("\n[SNF] Spectral clustering on fused network...")
    best_k = None
    best_score = -1.0
    best_labels = None

    D = 1.0 - W_fused
    D = (D + D.T) / 2.0
    np.fill_diagonal(D, 0.0)

    for k in candidate_K:
        sc = SpectralClustering(
            n_clusters=k,
            affinity="precomputed",
            random_state=random_state,
            assign_labels="kmeans"
        )
        labels = sc.fit_predict(W_fused)

        if len(set(labels)) <= 1:
            score = -1.0
        else:
            score = silhouette_score(D, labels, metric="precomputed")
        print(f"  [SNF] K={k}: silhouette_score={score:.4f}")

        if score > best_score:
            best_k = k
            best_score = score
            best_labels = labels

    print(f"[SNF] Best K = {best_k}, silhouette = {best_score:.4f}")
    return best_labels, best_k, best_score


# ============================================================
# 5. Early-fusion baseline
# ============================================================
def early_fusion_clustering(X_rna: pd.DataFrame,
                            X_mir: pd.DataFrame,
                            candidate_K=(3, 4, 5),
                            random_state: int = RANDOM_STATE):
    X_rna.index = X_rna.index.astype(str).str.strip()
    X_mir.index = X_mir.index.astype(str).str.strip()

    rna_keys = [simplify_sample_id(s) for s in X_rna.index]
    mir_keys = [simplify_sample_id(s) for s in X_mir.index]

    X_rna.index = rna_keys
    X_mir.index = mir_keys

    X_rna = X_rna[~pd.isna(X_rna.index)]
    X_mir = X_mir[~pd.isna(X_mir.index)]

    X_rna = X_rna[~X_rna.index.duplicated(keep="first")]
    X_mir = X_mir[~X_mir.index.duplicated(keep="first")]

    print(f"[DEBUG] Simplified RNA sample IDs (head): {list(X_rna.index[:5])}")
    print(f"[DEBUG] Simplified miRNA sample IDs (head): {list(X_mir.index[:5])}")

    common_samples = sorted(set(X_rna.index) & set(X_mir.index))
    print(f"[INFO] Common samples (RNA ∩ miRNA): {len(common_samples)}")

    if len(common_samples) == 0:
        print("[ERROR] RNA 和 miRNA 没有共同样本名，无法做多组学聚类。")
        return None, None, None, None, None

    X_rna = X_rna.loc[common_samples]
    X_mir = X_mir.loc[common_samples]

    X_rna_var = select_top_variable_features(X_rna, top_n=2000)
    X_mir_var = select_top_variable_features(X_mir, top_n=500)

    print("\n[PCA] RNA (early fusion)")
    Z_rna = pca_transform(X_rna_var, n_components=50)

    print("\n[PCA] miRNA (early fusion)")
    Z_mir = pca_transform(X_mir_var, n_components=20)

    Z_fused = np.hstack([Z_rna, Z_mir])
    print(f"[INFO] Early-fusion feature shape: {Z_fused.shape}")

    best_k = None
    best_score = -1.0
    best_labels = None

    print("\n[CLUSTER] Early-fusion KMeans, try different K")
    for k in candidate_K:
        km = KMeans(n_clusters=k, random_state=random_state, n_init=20)
        labels = km.fit_predict(Z_fused)

        if len(set(labels)) == 1:
            score = -1.0
        else:
            score = silhouette_score(Z_fused, labels)
        print(f"  [EarlyFusion] K={k}: silhouette_score={score:.4f}")

        if score > best_score:
            best_k = k
            best_score = score
            best_labels = labels

    print(f"[EarlyFusion] Best K = {best_k}, silhouette = {best_score:.4f}")

    return common_samples, Z_fused, best_labels, best_k, best_score


# ============================================================
# 6. UMAP / PCA 可视化
# ============================================================
def plot_embedding(Z_fused: np.ndarray,
                   labels_snf: np.ndarray,
                   labels_early: np.ndarray,
                   out_path: str,
                   random_state: int = RANDOM_STATE):
    import matplotlib.pyplot as plt

    if HAS_UMAP:
        print("[VIS] Using UMAP for 2D embedding...")
        reducer = umap.UMAP(random_state=random_state)
        emb = reducer.fit_transform(Z_fused)
        method_name = "UMAP"
    else:
        print("[VIS] UMAP not installed, fallback to PCA 2D.")
        pca2 = PCA(n_components=2, random_state=random_state)
        emb = pca2.fit_transform(Z_fused)
        method_name = "PCA2D"

    x = emb[:, 0]
    y = emb[:, 1]

    # SNF
    plt.figure(figsize=(6, 5))
    scatter = plt.scatter(x, y, c=labels_snf, cmap="tab10", s=40, alpha=0.8)
    plt.title(f"{method_name} embedding colored by SNF clusters")
    plt.xlabel("Dim 1")
    plt.ylabel("Dim 2")
    plt.colorbar(scatter, label="SNF Cluster")
    plt.tight_layout()
    snf_png = out_path.replace(".png", "_SNF.png")
    plt.savefig(snf_png, dpi=300)
    plt.close()
    print(f"[VIS] Saved SNF cluster embedding -> {snf_png}")

    # Early-fusion
    plt.figure(figsize=(6, 5))
    scatter = plt.scatter(x, y, c=labels_early, cmap="tab10", s=40, alpha=0.8)
    plt.title(f"{method_name} embedding colored by Early-fusion clusters")
    plt.xlabel("Dim 1")
    plt.ylabel("Dim 2")
    plt.colorbar(scatter, label="Early-fusion Cluster")
    plt.tight_layout()
    early_png = out_path.replace(".png", "_EarlyFusion.png")
    plt.savefig(early_png, dpi=300)
    plt.close()
    print(f"[VIS] Saved Early-fusion cluster embedding -> {early_png}")


# ============================================================
# 7. RNA DEG + pathway（基于 EarlyFusion cluster）
# ============================================================
from src.methods_utils import compute_deg, save_deg_tables, run_enrichr


def run_rna_deg_from_earlyfusion(X_rna_aligned: pd.DataFrame,
                                 labels_early: np.ndarray,
                                 sample_ids):
    print("\n[DEG] RNA DEG + pathway using Early-fusion clusters")

    X_rna_aligned = X_rna_aligned.loc[sample_ids]

    unique_clusters = sorted(set(labels_early))
    print(f"[DEG] Unique EarlyFusion clusters: {unique_clusters}")
    print("[DEG] Cluster sizes:")
    counts = Counter(labels_early)
    for c in unique_clusters:
        print(f"  - Cluster {c}: n = {counts[c]}")

    base_deg_dir = os.path.join(BASE_DIR, "data", "RNA-seq", "deg_EarlyFusion")
    base_path_dir = os.path.join(BASE_DIR, "data", "RNA-seq", "pathway_EarlyFusion")
    os.makedirs(base_deg_dir, exist_ok=True)
    os.makedirs(base_path_dir, exist_ok=True)

    for c in unique_clusters:
        print("\n" + "=" * 60)
        print(f"[CLUSTER-DEG] Processing EarlyFusion cluster {c} vs others")

        labels_binary = np.zeros_like(labels_early, dtype=int)
        labels_binary[labels_early == c] = 1

        n_in = int(labels_binary.sum())
        n_out = len(labels_binary) - n_in
        print(f"[INFO] Cluster {c}: n_in = {n_in}, n_out = {n_out}")
        if n_in < 3 or n_out < 3:
            print(f"[WARN] Cluster {c} 或其他组样本数过少（<3），DEG 结果可能不稳定。")

        deg_rna = compute_deg(
            X_rna_aligned,
            labels_binary,
            group1=0,
            group2=1
        )

        deg_dir = os.path.join(base_deg_dir, f"cluster{c}")
        path_dir = os.path.join(base_path_dir, f"cluster{c}")
        os.makedirs(deg_dir, exist_ok=True)
        os.makedirs(path_dir, exist_ok=True)
        print(f"[SAVE] DEG dir:      {deg_dir}")
        print(f"[SAVE] Pathway dir:  {path_dir}")

        data_type_label = f"RNA_EarlyFusion_cluster{c}"
        rna_up, rna_down, rna_relaxed = save_deg_tables(
            deg_rna,
            deg_dir,
            data_type_label,
            RANDOM_STATE
        )

        if len(rna_up) > 0:
            run_enrichr(
                rna_up,
                label=f"Up-regulated genes (EarlyFusion cluster {c}, relaxed)",
                out_prefix="Pathway_up",
                pathway_dir=path_dir,
                data_type=data_type_label,
                seed=RANDOM_STATE
            )
        else:
            print(f"[WARN] Cluster {c}: 没有 up-regulated 基因（按 relaxed 标准）。")

        if len(rna_down) > 0:
            run_enrichr(
                rna_down,
                label=f"Down-regulated genes (EarlyFusion cluster {c}, relaxed)",
                out_prefix="Pathway_down",
                pathway_dir=path_dir,
                data_type=data_type_label,
                seed=RANDOM_STATE
            )
        else:
            print(f"[WARN] Cluster {c}: 没有 down-regulated 基因（按 relaxed 标准）。")

        print(f"[DONE] Cluster {c}: RNA DEG + pathway enrichment finished.")

    print("\n[ALL DONE] 所有 EarlyFusion clusters 的 RNA DEG + pathway 分析完成。")


# ============================================================
# 8. 入口函数：外部脚本调用这个
# ============================================================
def run_snf_v2_pipeline():
    print("===== Multi-omics clustering V2.0 (SNF + Early-Fusion) =====")

    # 1) 加载表达矩阵
    X_rna = load_expression_only("RNA")
    X_mir = load_expression_only("smallRNA")

    # 2) Early-fusion baseline
    (
        common_samples,
        Z_fused,
        labels_early,
        best_k_early,
        score_early
    ) = early_fusion_clustering(
        X_rna,
        X_mir,
        candidate_K=(3, 4, 5),
        random_state=RANDOM_STATE
    )

    if common_samples is None:
        return

    # 3) SNF
    X_rna_aligned = X_rna.loc[common_samples]
    X_mir_aligned = X_mir.loc[common_samples]

    X_rna_var = select_top_variable_features(X_rna_aligned, top_n=2000)
    X_mir_var = select_top_variable_features(X_mir_aligned, top_n=500)

    W_rna = construct_affinity_matrix(X_rna_var.values, K=20, mu=0.5)
    W_mir = construct_affinity_matrix(X_mir_var.values, K=20, mu=0.5)

    W_fused = snf([W_rna, W_mir], K=20, t=20)

    labels_snf, best_k_snf, score_snf = snf_clustering(
        W_fused,
        candidate_K=(3, 4, 5),
        random_state=RANDOM_STATE
    )

    # 4) 保存聚类结果
    out_dir = os.path.join(BASE_DIR, "data", "integrated_results")
    os.makedirs(out_dir, exist_ok=True)

    snf_out_file = os.path.join(
        out_dir,
        f"cluster_results_SNF_RNA_miRNA_seed{RANDOM_STATE}.csv"
    )
    snf_df = pd.DataFrame({
        "Sample_ID": common_samples,
        "Cluster_SNF": labels_snf
    })
    snf_df.to_csv(snf_out_file, index=False)
    print(f"\n[SAVE] SNF cluster results -> {snf_out_file}")

    early_out_file = os.path.join(
        out_dir,
        f"cluster_results_EarlyFusion_RNA_miRNA_seed{RANDOM_STATE}.csv"
    )
    early_df = pd.DataFrame({
        "Sample_ID": common_samples,
        "Cluster_EarlyFusion": labels_early
    })
    early_df.to_csv(early_out_file, index=False)
    print(f"[SAVE] Early-fusion cluster results -> {early_out_file}")

    # 5) 可视化
    emb_png = os.path.join(
        out_dir,
        f"embedding_RNA_miRNA_seed{RANDOM_STATE}.png"
    )
    plot_embedding(
        Z_fused=Z_fused,
        labels_snf=labels_snf,
        labels_early=labels_early,
        out_path=emb_png,
        random_state=RANDOM_STATE
    )

    # 6) silhouette 对比
    print("\n[SUMMARY] Silhouette comparison")
    print(f"  SNF fused network:       K={best_k_snf},   silhouette={score_snf:.4f}")
    print(f"  Early-fusion PCA+KMeans: K={best_k_early}, silhouette={score_early:.4f}")

    summary_txt = os.path.join(
        out_dir,
        f"silhouette_comparison_seed{RANDOM_STATE}.txt"
    )
    with open(summary_txt, "w") as f:
        f.write("Silhouette comparison (SNF vs Early-fusion)\n")
        f.write(f"SNF fused network:       K={best_k_snf},   silhouette={score_snf:.4f}\n")
        f.write(f"Early-fusion PCA+KMeans: K={best_k_early}, silhouette={score_early:.4f}\n")
    print(f"[SAVE] Silhouette summary -> {summary_txt}")

    # 7) RNA DEG + pathway
    run_rna_deg_from_earlyfusion(
        X_rna_aligned=X_rna_aligned,
        labels_early=labels_early,
        sample_ids=common_samples
    )

    print("\n[DONE] SNF + Early-fusion clustering + visualization + RNA DEG/Pathway 完成。")
