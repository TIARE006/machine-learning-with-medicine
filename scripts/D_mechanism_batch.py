#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import gseapy as gp
import os, sys

# ===============================
# USAGE
# ===============================
if len(sys.argv) < 2:
    print("Usage: python D_mechanism_batch.py <RUN_DIR>")
    sys.exit(1)

RUN_DIR = sys.argv[1]

VST_PATH     = os.path.join(RUN_DIR, "degs_deseq2", "vst_matrix.csv")
CLUSTER_PATH = os.path.join(RUN_DIR, "labels", "cluster_results_RNA_seed42.csv")
DEG_DIR      = os.path.join(RUN_DIR, "degs_deseq2")

BASE_OUTDIR  = os.path.join(RUN_DIR, "pseudotime", "D_mechanism")
os.makedirs(BASE_OUTDIR, exist_ok=True)

# ===============================
# D MECHANISM CONFIG
# ===============================
D_CONFIG = {
    "0": {
        "TFs": ["NR3C2", "ESRRA"],
        "MARKERS": ["GAPDH", "TRIM72"],
        "MODULES": {
            "Energy": ["GAPDH", "PGK1", "ENO1"],
            "Muscle_repair": ["TRIM72", "DYSF"]
        }
    },
    "1": {
        "TFs": ["LIF", "STAT3"],
        "MARKERS": ["LIF"],
        "MODULES": {
            "JAK_STAT": ["LIF", "STAT3", "SOCS3", "IL6ST"],
            "Inflammation": ["NFKBIA", "CXCL8"]
        }
    },
    "2": {
        "TFs": ["ID2", "ID3"],
        "MARKERS": ["PIEZO1", "PIEZO2", "GSDMD", "MMP14"],
        "MODULES": {
            "Mechanical": ["PIEZO1", "PIEZO2"],
            "Inflammation": ["GSDMD", "IL6ST", "SOCS3"],
            "ECM": ["MMP14", "COL1A1", "CTGF"]
        }
    },
    "3": {
        "TFs": ["HMGCS2"],
        "MARKERS": ["HMGCS2"],
        "MODULES": {
            "FAO_Ketogenesis": ["HMGCS2", "CPT1A", "ACOX1"]
        }
    }
}

# ===============================
# LOAD DATA
# ===============================
vst = pd.read_csv(VST_PATH, index_col=0)
clusters = pd.read_csv(CLUSTER_PATH)

sample_col  = [c for c in clusters.columns if "sample" in c.lower() or "id" in c.lower()][0]
cluster_col = [c for c in clusters.columns if "cluster" in c.lower()][0]

clusters = clusters[[sample_col, cluster_col]].rename(
    columns={sample_col: "sample", cluster_col: "cluster"}
)
clusters["cluster"] = clusters["cluster"].astype(str)

expr = vst.T.merge(clusters, left_index=True, right_on="sample")

trrust = gp.get_library(
    name="TRRUST_Transcription_Factors_2019",
    organism="Human"
)

# ===============================
# MAIN LOOP
# ===============================
for cl, cfg in D_CONFIG.items():
    print(f"\n[INFO] Running D mechanism for cluster {cl}")
    OUT_DIR = os.path.join(BASE_OUTDIR, f"cluster{cl}")
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---------------------------
    # D1 | Module scores
    # ---------------------------
    module_names = []
    for name, genes in cfg["MODULES"].items():
        genes = [g for g in genes if g in expr.columns]
        if len(genes) < 2:
            continue
        z = (expr[genes] - expr[genes].mean()) / expr[genes].std()
        expr[name] = z.mean(axis=1)
        module_names.append(name)

    plot_df = expr.melt(
        id_vars=["cluster"],
        value_vars=module_names,
        var_name="Module",
        value_name="Score"
    )

    plt.figure(figsize=(7, 5))
    sns.boxplot(data=plot_df, x="cluster", y="Score", hue="Module")
    plt.title(f"D1 | Cluster {cl} pathway modules")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "D1_module_scores.png"), dpi=300)
    plt.close()

    # ---------------------------
    # D2 | TF activity
    # ---------------------------
    DEG_PATH = os.path.join(DEG_DIR, f"DEG_cluster_{cl}_vs_rest.csv")
    deg = pd.read_csv(DEG_PATH)

    # ---- robust column detection ----
    padj_candidates = ["padj", "FDR", "adj.P.Val", "p_adj"]
    lfc_candidates  = ["log2FoldChange", "log2FC", "logFC", "avg_log2FC", "LFC"]

    padj_col = next((c for c in padj_candidates if c in deg.columns), None)
    lfc_col  = next((c for c in lfc_candidates  if c in deg.columns), None)

    if padj_col is None or lfc_col is None:
        raise RuntimeError(
            f"[ERROR] DEG columns not found for cluster {cl}. "
            f"Columns = {deg.columns.tolist()}"
        )

    deg_up = deg.loc[
        (deg[padj_col] < 0.05) & (deg[lfc_col] > 0),
        "gene"
    ].dropna().tolist()

    tf_cols = []
    for tf in cfg["TFs"]:
        targets = set(trrust.get(tf, []))
        module_genes = list((targets | set(deg_up)) & set(expr.columns))
        if len(module_genes) < 5:
            continue
        z = (expr[module_genes] - expr[module_genes].mean()) / expr[module_genes].std()
        expr[f"{tf}_activity"] = z.mean(axis=1)
        tf_cols.append(f"{tf}_activity")

    plot_df = expr.melt(
        id_vars=["cluster"],
        value_vars=tf_cols,
        var_name="TF",
        value_name="Activity"
    )

    plt.figure(figsize=(6, 4))
    sns.boxplot(data=plot_df, x="cluster", y="Activity", hue="TF")
    plt.title(f"D2 | Cluster {cl} TF activity")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "D2_TF_activity.png"), dpi=300)
    plt.close()

    # ---------------------------
    # D3 | Targets + marker heatmap
    # ---------------------------
    expr_c = expr[expr["cluster"] == cl]

    genes = set(cfg["MARKERS"])
    for tf in cfg["TFs"]:
        genes |= set(trrust.get(tf, []))
    genes |= set(deg_up)

    genes = [g for g in genes if g in expr_c.columns]

    mat = expr_c[genes]
    mat_z = (mat - mat.mean()) / mat.std()

    sns.clustermap(
        mat_z.T,
        cmap="vlag",
        col_cluster=True,
        yticklabels=True,
        figsize=(10, 6)
    )
    plt.savefig(os.path.join(OUT_DIR, "D3_targets_markers_heatmap.png"), dpi=300)
    plt.close()

print("\n[D] All cluster mechanisms completed.")
