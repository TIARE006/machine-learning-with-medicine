#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os, sys

# ===============================
# USAGE
# ===============================
if len(sys.argv) < 2:
    print("Usage: python C_final_markers.py <RUN_DIR>")
    sys.exit(1)

RUN_DIR = sys.argv[1]
OUT_DIR = os.path.join(RUN_DIR, "pseudotime", "C_markers")
os.makedirs(OUT_DIR, exist_ok=True)

VST_PATH     = os.path.join(RUN_DIR, "degs_deseq2", "vst_matrix.csv")
CLUSTER_PATH = os.path.join(RUN_DIR, "labels", "cluster_results_RNA_seed42.csv")

ORDER = ["0", "3", "1", "2"]

MARKERS = {
    "GAPDH": "Cluster0",
    "TRIM72": "Cluster0",
    "LIF": "Cluster1",
    "PIEZO1": "Cluster2",
    "PIEZO2": "Cluster2",
    "GSDMD": "Cluster2",
    "MMP14": "Cluster2",
    "HMGCS2": "Cluster3"
}

# ===============================
# Load data
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

# ===============================
# Plot per gene
# ===============================
for gene, cl in MARKERS.items():
    if gene not in expr.columns:
        print(f"[WARN] {gene} not found, skipped")
        continue

    df = expr[["cluster", gene]].copy()
    df["cluster"] = pd.Categorical(df["cluster"], ORDER)

    # ---------------------------
    # C1 | Mean ± SEM (ordered)
    # ---------------------------
    stats = df.groupby("cluster")[gene].agg(
        mean="mean",
        sem=lambda x: x.std() / np.sqrt(len(x))
    ).loc[ORDER]

    plt.figure(figsize=(4, 3))
    plt.errorbar(
        x=ORDER,
        y=stats["mean"],
        yerr=stats["sem"],
        fmt="-o",
        capsize=4
    )
    plt.title(f"{gene} | ordered progression")
    plt.ylabel("Expression (VST)")
    plt.tight_layout()
    plt.savefig(
        os.path.join(OUT_DIR, f"C1_{gene}_meanSEM.png"),
        dpi=300
    )
    plt.close()

    # ---------------------------
    # C2 | Boxplot by cluster
    # ---------------------------
    plt.figure(figsize=(4, 3))
    sns.boxplot(data=df, x="cluster", y=gene, order=ORDER)
    sns.stripplot(
        data=df, x="cluster", y=gene,
        order=ORDER, color="black", size=3, alpha=0.6
    )
    plt.title(f"{gene} | by cluster")
    plt.tight_layout()
    plt.savefig(
        os.path.join(OUT_DIR, f"C2_{gene}_boxplot.png"),
        dpi=300
    )
    plt.close()

print("[C] Marker gene figures completed.")
print("Saved to:", OUT_DIR)
