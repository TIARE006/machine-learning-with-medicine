import os
import hashlib
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import SpectralClustering, AgglomerativeClustering
from sklearn.metrics import silhouette_score

# ============================================================
# 目标：为 DEG / 通路富集服务的 RNA 聚类
# 关键：logCPM -> HVG -> PCA -> (Spectral/Ward) -> 选K -> 输出结果
#
# 输出严格按你现有格式：
# results/clustering/{RNA_only|smallRNA_only}/run_YYYYMMDD-HHMMSS__hash/
#   - labels/ cluster_results_{DATA_TYPE}_seed{RANDOM_STATE}.csv  (Sample_ID, Cluster)
#   - metrics/ params/ plots/ ... （按你旧结构分目录）
# ============================================================

# -------------------------
# 配置区：按需改
# -------------------------
DATA_TYPE = "RNA"  # "RNA" or "smallRNA"
RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

# 聚类候选K范围
K_MIN = 2
K_MAX = 7

# HVG数量（建议 1000~3000）
N_HVG = 2000

# PCA维度（建议 20~50）
PCA_NCOMP = 30

# 共识聚类参数
CONS_N_ITER = 100
SUBSAMPLE_RATE = 0.85

# 聚类方法： "spectral" 或 "ward"
CLUSTER_METHOD = "spectral"

# Spectral 的邻居数（一般 8~15）
SPECTRAL_N_NEIGHBORS = 10

# 低表达基因过滤（CPM空间）
MIN_CPM_IN_AT_LEAST_N_SAMPLES = 5  # 至少在多少个样本中 CPM > MIN_CPM 才保留
MIN_CPM = 1.0

# 项目根目录（脚本在 scripts/ 下）
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def load_counts_matrix(expr_file: str) -> pd.DataFrame:
    """
    读入 raw counts CSV，返回 基因×样本 的 counts 矩阵（index=gene，columns=samples）
    兼容：第一行异常描述行 + 第一列 gene id + optional type 列
    """
    df = pd.read_csv(expr_file, low_memory=False)
    df = df.drop(index=0)

    gene_col = df.columns[0]
    df = df.set_index(gene_col)

    if "type" in df.columns:
        df = df.drop(columns=["type"])

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    df = df.loc[~(df.sum(axis=1) == 0)]  # 删除全0基因
    return df


def counts_to_logcpm(counts_gxs: pd.DataFrame) -> pd.DataFrame:
    """
    counts -> CPM -> log1p(CPM)
    返回：logCPM（基因×样本）
    """
    counts = counts_gxs.values.astype(float)
    libsize = counts.sum(axis=0, keepdims=True) + 1e-12
    cpm = counts / libsize * 1e6
    logcpm = np.log1p(cpm)
    return pd.DataFrame(logcpm, index=counts_gxs.index, columns=counts_gxs.columns)


def filter_low_expression(logcpm_gxs: pd.DataFrame) -> pd.DataFrame:
    """
    低表达过滤：至少在 N 个样本中 CPM > MIN_CPM 才保留
    在 logCPM 上判断：log1p(CPM) > log1p(MIN_CPM)
    """
    thr = np.log1p(MIN_CPM)
    keep = (logcpm_gxs > thr).sum(axis=1) >= MIN_CPM_IN_AT_LEAST_N_SAMPLES
    return logcpm_gxs.loc[keep]


def select_hvg(logcpm_gxs: pd.DataFrame, n_hvg: int = 2000) -> pd.DataFrame:
    """
    HVG选择：按基因在样本间的方差排序取前 n_hvg
    """
    gene_var = logcpm_gxs.var(axis=1)
    n = min(n_hvg, gene_var.shape[0])
    top_genes = gene_var.sort_values(ascending=False).head(n).index
    return logcpm_gxs.loc[top_genes]


def cluster_fit(X_pca: np.ndarray, k: int, method: str) -> np.ndarray:
    """
    在 PCA 空间做聚类
    """
    if method == "spectral":
        model = SpectralClustering(
            n_clusters=k,
            affinity="nearest_neighbors",
            n_neighbors=min(SPECTRAL_N_NEIGHBORS, max(1, X_pca.shape[0] - 1)),
            random_state=RANDOM_STATE,
            assign_labels="kmeans",
        )
        return model.fit_predict(X_pca)

    if method == "ward":
        model = AgglomerativeClustering(n_clusters=k, linkage="ward")
        return model.fit_predict(X_pca)

    raise ValueError("CLUSTER_METHOD must be 'spectral' or 'ward'")


def consensus_pac_scores(
    X_pca: np.ndarray,
    k_min: int,
    k_max: int,
    n_iter: int,
    subsample_rate: float,
    method: str,
):
    """
    共识聚类（带共同出现归一化）+ PAC 指标（越小越好）
    - same[i,j]：同次抽样且同类计数
    - together[i,j]：同次抽样共同出现计数
    - consensus = same/together (where together>0)
    PAC：consensus 值落在 (0.1, 0.9) 的比例（越小越不模糊）
    """
    n = X_pca.shape[0]
    results = {}

    u, v = 0.1, 0.9

    for k in range(k_min, k_max + 1):
        same = np.zeros((n, n), dtype=float)
        together = np.zeros((n, n), dtype=float)

        for _ in range(n_iter):
            idx = np.random.choice(n, int(np.ceil(subsample_rate * n)), replace=False)
            X_sub = X_pca[idx, :]

            labels = cluster_fit(X_sub, k, method)

            # 更新 together（共同出现）
            for a in range(len(idx)):
                ia = idx[a]
                together[ia, idx] += 1

            # 更新 same（同类）
            for c in np.unique(labels):
                members = idx[labels == c]
                same[np.ix_(members, members)] += 1

        consensus = np.full((n, n), np.nan, dtype=float)
        mask = together > 0
        consensus[mask] = same[mask] / together[mask]

        triu = np.triu_indices(n, k=1)
        vals = consensus[triu]
        vals = vals[~np.isnan(vals)]

        pac = np.mean((vals > u) & (vals < v)) if vals.size > 0 else np.nan
        mean_cons = np.mean(vals) if vals.size > 0 else np.nan

        results[k] = {"PAC": pac, "mean_consensus": mean_cons}

    return results


def make_sample_correlation_heatmap(logcpm_hvg_gxs: pd.DataFrame, labels: np.ndarray, out_png: str):
    """
    用 HVG 的 logCPM 计算样本-样本相关性热图，并按 cluster 排序
    """
    X = logcpm_hvg_gxs.T.values  # samples x genes
    corr = np.corrcoef(X)

    order = np.argsort(labels)
    corr_ord = corr[np.ix_(order, order)]

    plt.figure(figsize=(8, 7))
    plt.imshow(corr_ord, aspect="auto")
    plt.colorbar(label="Pearson correlation")
    plt.title("Sample-Sample Correlation (HVG logCPM), ordered by cluster")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    return order


def make_run_dir(base_dir: str, data_type: str, random_state: int, best_k: int, method: str) -> str:
    """
    results/clustering/{RNA_only|smallRNA_only}/run_YYYYMMDD-HHMMSS__hash/
    """
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    hash_src = f"{data_type}|seed={random_state}|k={best_k}|method={method}"
    run_hash = hashlib.md5(hash_src.encode("utf-8")).hexdigest()[:10]

    run_dir = os.path.join(
        base_dir,
        "results",
        "clustering",
        f"{data_type}_only",
        f"run_{timestamp}__{run_hash}",
    )
    os.makedirs(run_dir, exist_ok=True)
    return run_dir


def main():
    # -------------------------
    # 输入文件路径（保持你原项目结构）
    # -------------------------
    if DATA_TYPE == "smallRNA":
        expr_file = os.path.join(
            BASE_DIR,
            "data",
            "raw", "small_RNA_seq",
            
            "GSE254878_smallRNAs_raw_counts_expression.csv",
        )
    elif DATA_TYPE == "RNA":
        expr_file = os.path.join(
            BASE_DIR,
            "data",
            "raw", "RNA_seq",
            
            "GSE254877_raw_counts_expression.csv",
        )
    else:
        raise ValueError("DATA_TYPE must be 'smallRNA' or 'RNA'")

    print(f"[INFO] Mode: {DATA_TYPE}")
    print(f"[INFO] Input: {expr_file}")

    # -------------------------
    # 1) counts
    # -------------------------
    counts = load_counts_matrix(expr_file)
    print("[INFO] counts shape (genes x samples):", counts.shape)

    # -------------------------
    # 2) logCPM + filter
    # -------------------------
    logcpm = counts_to_logcpm(counts)
    logcpm = filter_low_expression(logcpm)
    print("[INFO] after low-expression filter:", logcpm.shape)

    # -------------------------
    # 3) HVG
    # -------------------------
    logcpm_hvg = select_hvg(logcpm, n_hvg=N_HVG)
    print("[INFO] HVG logcpm shape:", logcpm_hvg.shape)

    # -------------------------
    # 4) PCA
    # -------------------------
    X_df = logcpm_hvg.T  # samples x genes
    X_scaled = StandardScaler().fit_transform(X_df.values)

    ncomp = min(PCA_NCOMP, X_scaled.shape[1], max(1, X_scaled.shape[0] - 1))
    pca = PCA(n_components=ncomp, random_state=RANDOM_STATE)
    X_pca = pca.fit_transform(X_scaled)

    # -------------------------
    # 5) 选K（PAC + silhouette）
    # -------------------------
    print("[INFO] Running consensus clustering for K selection ...")
    cons = consensus_pac_scores(
        X_pca,
        k_min=K_MIN,
        k_max=K_MAX,
        n_iter=CONS_N_ITER,
        subsample_rate=SUBSAMPLE_RATE,
        method=CLUSTER_METHOD,
    )

    k_list, pac_list, mean_cons_list, sil_list = [], [], [], []
    for k in range(K_MIN, K_MAX + 1):
        labels_k = cluster_fit(X_pca, k, CLUSTER_METHOD)
        try:
            sil = silhouette_score(X_pca, labels_k, metric="euclidean")
        except Exception:
            sil = np.nan

        k_list.append(k)
        pac_list.append(cons[k]["PAC"])
        mean_cons_list.append(cons[k]["mean_consensus"])
        sil_list.append(sil)

    metrics = pd.DataFrame(
        {
            "K": k_list,
            "PAC_lower_better": pac_list,
            "mean_consensus": mean_cons_list,
            "silhouette_higher_better": sil_list,
        }
    )
    print("[INFO] K-selection metrics:\n", metrics)

    best_candidates = metrics.dropna(subset=["PAC_lower_better"]).copy()
    min_pac = best_candidates["PAC_lower_better"].min()
    sub = best_candidates[best_candidates["PAC_lower_better"] == min_pac].copy()

    if sub["silhouette_higher_better"].notna().any():
        best_k = sub.sort_values(["silhouette_higher_better", "K"], ascending=[False, True]).iloc[0]["K"]
    else:
        best_k = sub.sort_values(["K"], ascending=[True]).iloc[0]["K"]
    best_k = int(best_k)

    print(f"[INFO] Selected best K = {best_k} (method={CLUSTER_METHOD})")

    # -------------------------
    # 6) 最终聚类
    # -------------------------
    cluster_labels = cluster_fit(X_pca, best_k, CLUSTER_METHOD)

    # -------------------------
    # 7) run 目录 + 子目录（对齐你旧结构）
    # -------------------------
    run_dir = make_run_dir(BASE_DIR, DATA_TYPE, RANDOM_STATE, best_k, CLUSTER_METHOD)

    labels_dir   = os.path.join(run_dir, "labels")
    metrics_dir  = os.path.join(run_dir, "metrics")
    params_dir   = os.path.join(run_dir, "params")
    plots_dir    = os.path.join(run_dir, "plots")

    degs_dir     = os.path.join(run_dir, "degs_deseq2")
    pathways_dir = os.path.join(run_dir, "pathways")
    heatmaps_dir = os.path.join(run_dir, "heatmaps")

    for d in [labels_dir, metrics_dir, params_dir, plots_dir, degs_dir, pathways_dir, heatmaps_dir]:
        os.makedirs(d, exist_ok=True)

    print(f"[INFO] Output run dir: {run_dir}")

    # -------------------------
    # 8) labels：核心接口（Sample_ID, Cluster）
    # -------------------------
    result = pd.DataFrame({"Sample_ID": X_df.index, "Cluster": cluster_labels})
    result_path = os.path.join(labels_dir, f"cluster_results_{DATA_TYPE}_seed{RANDOM_STATE}.csv")
    result.to_csv(result_path, index=False)
    print(f"Cluster results saved to: {result_path}")

    # ordered samples（给热图排序用）
    order = make_sample_correlation_heatmap(
        logcpm_hvg_gxs=logcpm_hvg,
        labels=cluster_labels,
        out_png=os.path.join(plots_dir, "sample_correlation_heatmap.png"),
    )
    ordered_samples = [X_df.index[i] for i in order]
    pd.DataFrame({"Sample_ID": ordered_samples}).to_csv(
        os.path.join(labels_dir, "ordered_samples_for_heatmap.csv"), index=False
    )

    # -------------------------
    # 9) params
    # -------------------------
    logcpm_hvg.index.to_series().to_csv(os.path.join(params_dir, "HVG_genes.txt"), index=False, header=False)
    pd.Series(pca.explained_variance_ratio_, name="explained_variance_ratio").to_csv(
        os.path.join(params_dir, "pca_explained_variance_ratio.csv"), index=False
    )

    # -------------------------
    # 10) metrics
    # -------------------------
    metrics.to_csv(os.path.join(metrics_dir, "K_selection_metrics.csv"), index=False)

    plt.figure()
    plt.plot(metrics["K"], metrics["PAC_lower_better"], marker="o")
    plt.xlabel("K")
    plt.ylabel("PAC (lower is better)")
    plt.title(f"Consensus PAC vs K ({DATA_TYPE}, {CLUSTER_METHOD})")
    plt.tight_layout()
    plt.savefig(os.path.join(metrics_dir, "consensus_PAC_vs_K.png"), dpi=200)
    plt.close()

    plt.figure()
    plt.plot(metrics["K"], metrics["silhouette_higher_better"], marker="o")
    plt.xlabel("K")
    plt.ylabel("Silhouette (higher is better)")
    plt.title(f"Silhouette vs K ({DATA_TYPE}, {CLUSTER_METHOD})")
    plt.tight_layout()
    plt.savefig(os.path.join(metrics_dir, "silhouette_vs_K.png"), dpi=200)
    plt.close()

    # -------------------------
    # 11) plots
    # -------------------------
    plt.figure()
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=cluster_labels)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(f"PCA (HVG logCPM) - {DATA_TYPE} - K={best_k} - {CLUSTER_METHOD}")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "PCA_clusters.png"), dpi=200)
    plt.close()

    pca_coords = (
        pd.DataFrame(X_pca[:, :2], columns=["PC1", "PC2"], index=X_df.index)
        .reset_index()
        .rename(columns={"index": "Sample_ID"})
    )
    pca_coords["Cluster"] = cluster_labels
    pca_coords.to_csv(os.path.join(plots_dir, "pca_coords_PC1_PC2.csv"), index=False)

    print("[DONE] Key outputs:")
    print(f" - labels/{os.path.basename(result_path)} (用于 DESeq2 分组，接口不变)")
    print(" - labels/ordered_samples_for_heatmap.csv")
    print(" - params/HVG_genes.txt")
    print(" - metrics/K_selection_metrics.csv")
    print(" - plots/PCA_clusters.png")
    print(" - plots/sample_correlation_heatmap.png")


if __name__ == "__main__":
    main()

