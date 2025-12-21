# src/rna_plotting.py
# -*- coding: utf-8 -*-

import re
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind

try:
    from upsetplot import UpSet, from_contents
    HAS_UPSETPLOT = True
except ImportError:
    HAS_UPSETPLOT = False


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def compute_degs_ttest_log1p(
    X_counts: pd.DataFrame,  # samples × genes
    labels: np.ndarray,
    min_cluster_size: int = 3,
    top_m: int = 50,
    fdr_thresh: float = 0.05,
    lfc_thresh: float = 1.0,
):
    """
    轻量 DEG：log1p(counts) + Welch t-test + BH-FDR
    返回:
      deg_sets: dict cluster -> set(genes)
      marker_genes: list (union top markers)
      deg_tables: dict cluster -> df (gene, log2FC, p, q)
    """
    X = np.log1p(X_counts.astype(float))
    clusters = np.unique(labels)

    deg_sets, deg_tables = {}, {}
    marker_union = []

    for c in clusters:
        idx_in = labels == c
        idx_out = labels != c

        n_in, n_out = int(idx_in.sum()), int(idx_out.sum())
        if n_in < min_cluster_size or n_out < min_cluster_size:
            deg_sets[int(c)] = set()
            deg_tables[int(c)] = pd.DataFrame(columns=["gene", "log2FC", "p", "q"])
            continue

        Xin = X.loc[idx_in]
        Xout = X.loc[idx_out]

        mu_in = Xin.mean(axis=0)
        mu_out = Xout.mean(axis=0)

        eps = 1e-8
        log2fc = np.log2((mu_in + eps) / (mu_out + eps))

        _, pvals = ttest_ind(
            Xin.values, Xout.values,
            axis=0, equal_var=False, nan_policy="omit"
        )
        qvals = _bh_fdr(pvals)

        df = pd.DataFrame({
            "gene": X.columns,
            "log2FC": log2fc.values,
            "p": pvals,
            "q": qvals
        }).sort_values(["q", "p"], ascending=True)

        df_deg = df[(df["q"] <= fdr_thresh) & (np.abs(df["log2FC"]) >= lfc_thresh)].copy()
        deg_sets[int(c)] = set(df_deg["gene"].tolist())
        deg_tables[int(c)] = df_deg

        df_top = df_deg.sort_values(["q", "p", "log2FC"], ascending=[True, True, False]).head(top_m)
        marker_union.extend(df_top["gene"].tolist())

    # union markers, de-duplicate keep order
    marker_genes = list(dict.fromkeys(marker_union))
    return deg_sets, marker_genes, deg_tables


def plot_expression_heatmap(
    X_counts: pd.DataFrame,   # samples × genes
    labels: np.ndarray,
    genes: list,
    out_png: str,
    max_genes: int = 80,          # 论文图别太多
    show_gene_labels: int = 30,   # 只显示前30个基因名
):
    if len(genes) == 0:
        print("[WARN] No marker genes; skip heatmap.")
        return

    # 选基因（最多 max_genes）
    genes = genes[:max_genes]

    X = np.log1p(X_counts.loc[:, genes].astype(float))

    # z-score per gene across samples
    Xz = pd.DataFrame(
        StandardScaler().fit_transform(X.values),
        index=X.index,
        columns=X.columns
    )

    # 按 cluster 排序样本
    order = np.argsort(labels)
    Xz = Xz.iloc[order]
    labels_ord = labels[order]

    data = Xz.T  # genes × samples

    # 顶部 cluster 色条
    uniq = sorted(np.unique(labels_ord))
    palette = sns.color_palette("Set2", n_colors=len(uniq))
    cluster_to_color = {c: palette[i] for i, c in enumerate(uniq)}
    col_colors = [cluster_to_color[int(c)] for c in labels_ord]

    # 控制 y tick：只显示前 show_gene_labels 个基因名
    yticklabels = [g if i < show_gene_labels else "" for i, g in enumerate(data.index)]

    sns.set_theme(context="paper", style="white", font_scale=1.1)

    g = sns.clustermap(
        data,
        col_cluster=False,
        row_cluster=True,
        col_colors=col_colors,
        xticklabels=False,
        yticklabels=yticklabels,
        cmap="vlag",
        figsize=(12, 8),
        dendrogram_ratio=(0.12, 0.02),   # 缩小树状图占比
        cbar_pos=(0.02, 0.8, 0.03, 0.18) # 把色条放左上角紧凑位置
    )

    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.fig.suptitle("RNA marker-gene heatmap (log1p + z-score)", y=1.02, fontsize=14)

    g.savefig(out_png, dpi=400, bbox_inches="tight")
    plt.close("all")



def plot_deg_upset(
    deg_sets: dict,
    out_png: str,
    min_subset_size: int = 10,
    max_intersections: int = 12,
):
    if not HAS_UPSETPLOT:
        print("[WARN] upsetplot not installed; skip upset plot.")
        return

    contents = {f"Cluster {k}": v for k, v in deg_sets.items() if len(v) > 0}
    if len(contents) < 2:
        print("[WARN] Not enough DEG sets for upset; skip.")
        return

    upset_data = from_contents(contents)

    sns.set_theme(context="paper", style="white", font_scale=1.1)
    plt.figure(figsize=(10, 5))

    up = UpSet(
        upset_data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",          # 交集大小排序
        min_subset_size=min_subset_size,
    )
    up.plot()

    plt.suptitle("DEG overlap across clusters (UpSet)", y=1.02, fontsize=13)
    plt.savefig(out_png, dpi=400, bbox_inches="tight")
    plt.close("all")

def _is_clinical_gene(g: str) -> bool:
    # 过滤难以临床验证/注释弱的ID（可按你需要调整）
    g = str(g)
    return not re.match(r"^(AC\d+|AL\d+|FP\d+|LOC\d+|LINC\d+)\b", g)

def load_deseq2_deg_sets_and_top_genes(
    degs_dir: str,
    top_per_cluster: int = 12,
    fdr: float = 0.05,
    lfc: float = 1.0,
    require_clinical_names: bool = True,
    upset_top_n: int = 800,          # 新增：UpSet 每个 cluster 只取 Top N 显著基因
    upset_require_clinical: bool = True,  # 新增：UpSet 也只保留“临床友好”基因名
):
    """
    读取 R 输出的 DESeq2 one-vs-rest 结果：
      - deg_sets: cluster -> set(genes) （用于 UpSet）
      - top_genes: 用于 heatmap 的 gene list（每个 cluster Top N，去重保序）
    """
    degs_dir = Path(degs_dir)
    deg_sets = {}
    top_genes = []

    for f in sorted(degs_dir.glob("DEG_cluster_*_vs_rest.csv")):
        # 文件名: DEG_cluster_<k>_vs_rest.csv
        k = int(f.name.split("_")[2])
        df = pd.read_csv(f)

        # 清理
        df = df.dropna(subset=["FDR", "log2FC"])
        df_sig = df[(df["FDR"] <= fdr) & (df["log2FC"].abs() >= lfc)].copy()

        # --- UpSet 用：只取 Top N（避免 13080 这种巨量吞掉信息） ---
        df_up = df_sig.sort_values(["FDR", "p_value"], ascending=True)
        if upset_top_n is not None and upset_top_n > 0:
            df_up = df_up.head(int(upset_top_n))

        genes_up = df_up["gene"].astype(str).tolist()
        if upset_require_clinical:
            genes_up = [g for g in genes_up if _is_clinical_gene(g)]
        deg_sets[k] = set(genes_up)

        # --- Heatmap marker：每类取 Top N ---
        df_top = df_sig.sort_values(["FDR", "p_value"], ascending=True).head(top_per_cluster)
        genes = df_top["gene"].astype(str).tolist()
        if require_clinical_names:
            genes = [g for g in genes if _is_clinical_gene(g)]
        top_genes.extend(genes)

    # 去重保序
    top_genes = list(dict.fromkeys(top_genes))
    return deg_sets, top_genes


def plot_heatmap_from_vst(
    vst_csv: str,
    labels_df: pd.DataFrame,  # columns: Sample_ID, Cluster
    genes: list,
    out_png: str,
    max_genes: int = 40,
    show_gene_labels: int = 30,
):
    """
    用 DESeq2 的 VST matrix 画 heatmap（临床/论文更认可）
    vst_matrix.csv 格式：gene × samples（第一列 gene）
    """
    vst = pd.read_csv(vst_csv)
    vst = vst.set_index("gene")
    # 转成 samples × genes
    X = vst.T
    X.index = X.index.astype(str)

    # 对齐样本顺序
    lab = labels_df.copy()
    lab["Sample_ID"] = lab["Sample_ID"].astype(str)
    lab = lab[lab["Sample_ID"].isin(X.index)]
    X = X.loc[lab["Sample_ID"]]

    # 选基因
    genes_use = [g for g in genes if g in X.columns][:max_genes]
    if len(genes_use) == 0:
        print("[WARN] No genes matched VST columns; skip heatmap.")
        return

    Xsub = X[genes_use]

    # z-score per gene（跨样本）让图更清晰
    Xz = pd.DataFrame(
        StandardScaler().fit_transform(Xsub.values),
        index=Xsub.index,
        columns=Xsub.columns
    )

    # 按 cluster 排序样本
    order = np.argsort(lab["Cluster"].values)
    Xz = Xz.iloc[order]
    lab_ord = lab.iloc[order].reset_index(drop=True)

    data = Xz.T  # genes × samples

    uniq = sorted(lab_ord["Cluster"].unique())
    palette = sns.color_palette("Set2", n_colors=len(uniq))
    cluster_to_color = {c: palette[i] for i, c in enumerate(uniq)}
    col_colors = [cluster_to_color[int(c)] for c in lab_ord["Cluster"].values]

    yticklabels = [g if i < show_gene_labels else "" for i, g in enumerate(data.index)]

    sns.set_theme(context="paper", style="white", font_scale=1.1)
    g = sns.clustermap(
        data,
        col_cluster=False,
        row_cluster=True,
        col_colors=col_colors,
        xticklabels=False,
        yticklabels=yticklabels,
        cmap="vlag",
        figsize=(12, 8),
        dendrogram_ratio=(0.12, 0.02),
        cbar_pos=(0.02, 0.8, 0.03, 0.18)
    )
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.fig.suptitle("RNA marker heatmap (DESeq2 one-vs-rest, VST + z-score)", y=1.02, fontsize=14)
    g.savefig(out_png, dpi=400, bbox_inches="tight")
    plt.close("all")


def make_rna_cluster_plots(
    run_dir: str,
    load_expression_only_fn,
    top_m: int = 50,     # 这里不再用于 t-test；保留参数兼容
    fdr: float = 0.05,
    lfc: float = 1.0,
    max_genes_heatmap: int = 40,   # 临床热图建议 30–50
    top_per_cluster: int = 12,     # 每个 cluster 取多少个 marker（临床建议 10–15）
):
    """
    论文/临床级绘图（Option A）：
      - DEG：使用 run_dir/degs_deseq2/ 下的 DESeq2 one-vs-rest 结果
      - Heatmap：使用 DESeq2 的 VST 矩阵 vst_matrix.csv
      - UpSet：使用 DESeq2 的显著 DEG set（FDR + |log2FC|）
    """
    run_dir = Path(run_dir)
    labels_csv = run_dir / "labels" / "cluster_labels.csv"
    plots_dir = run_dir / "plots"
    degs_dir = run_dir / "degs_deseq2"   # 你 R 脚本输出目录
    plots_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1) Load labels ----
    lab = pd.read_csv(labels_csv)
    lab["Sample_ID"] = lab["Sample_ID"].astype(str)

    # ---- 2) Check DESeq2 outputs exist ----
    vst_csv = degs_dir / "vst_matrix.csv"
    if not vst_csv.exists():
        raise FileNotFoundError(f"Missing VST matrix: {vst_csv}. Run the R DESeq2 script first.")

    deg_files = list(degs_dir.glob("DEG_cluster_*_vs_rest.csv"))
    if len(deg_files) == 0:
        raise FileNotFoundError(f"No DESeq2 DEG files found in: {degs_dir}. Run the R DESeq2 script first.")

    # ---- 3) Load DESeq2 DEG sets + select top marker genes ----
    deg_sets, top_genes = load_deseq2_deg_sets_and_top_genes(
        degs_dir=str(degs_dir),
        top_per_cluster=top_per_cluster,
        fdr=fdr,
        lfc=lfc,
        require_clinical_names=True,
        upset_top_n=800,               
        upset_require_clinical=True
    )

    # ================================
    # 只保留 cluster-specific DEG
    # ================================
    # specific_deg_sets = {}
    # for k, genes in deg_sets.items():
    #     other_genes = set().union(
    #         *[v for kk, v in deg_sets.items() if kk != k]
    #     )
    #     specific = genes - other_genes
    #     specific_deg_sets[k] = specific

    # 用于 UpSet（展示交集）
    deg_sets_all = deg_sets

    # 用于 heatmap（展示特异 marker）
    specific_deg_sets = {}
    for k, genes in deg_sets.items():
        other_genes = set().union(
            *[v for kk, v in deg_sets.items() if kk != k]
        )
        specific_deg_sets[k] = genes - other_genes


    # 可选：打印每个 cluster 的 specific DEG 数量（调试用）
    for k, v in specific_deg_sets.items():
        print(f"[INFO] Cluster {k} specific DEG:", len(v))


    # ---- 4) Heatmap from VST ----
    plot_heatmap_from_vst(
        vst_csv=str(vst_csv),
        labels_df=lab,
        genes=top_genes,
        out_png=str(plots_dir / "heatmap_marker_genes_DESeq2.png"),
        max_genes=max_genes_heatmap,
        show_gene_labels=min(30, max_genes_heatmap)
    )

    # ---- 5) UpSet from DESeq2 DEG sets ----
    plot_deg_upset(
        deg_sets=deg_sets_all,
        out_png=str(plots_dir / "upset_DEG_overlap_DESeq2_specific.png"),
        min_subset_size=5
    )

    print(f"[SAVE] plots -> {plots_dir}")
    print(f"[INFO] DESeq2 degs dir -> {degs_dir}")
    return str(plots_dir), str(degs_dir)
