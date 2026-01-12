#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys

# ===============================
# USAGE
# ===============================
if len(sys.argv) < 2:
    print("Usage: python A2_branching_signature.py <RUN_DIR>")
    sys.exit(1)

RUN_DIR = sys.argv[1]

VST_PATH = os.path.join(RUN_DIR, "degs_deseq2", "vst_matrix.csv")
LABELS   = os.path.join(RUN_DIR, "labels", "cluster_results_RNA_seed42.csv")

DEG1 = os.path.join(RUN_DIR, "degs_deseq2", "DEG_cluster_1_vs_rest.csv")
DEG2 = os.path.join(RUN_DIR, "degs_deseq2", "DEG_cluster_2_vs_rest.csv")

OUT_DIR = os.path.join(RUN_DIR, "pseudotime", "A2_branching")
os.makedirs(OUT_DIR, exist_ok=True)

ORDER = ["0", "3", "1", "2"]

# ===============================
# Load expression & labels
# ===============================
vst = pd.read_csv(VST_PATH, index_col=0)
labels = pd.read_csv(LABELS)

sample_col  = [c for c in labels.columns if "sample" in c.lower() or "id" in c.lower()][0]
cluster_col = [c for c in labels.columns if "cluster" in c.lower()][0]

labels = labels[[sample_col, cluster_col]].rename(
    columns={sample_col: "sample", cluster_col: "cluster"}
)
labels["cluster"] = labels["cluster"].astype(str)

expr = vst.T.merge(labels, left_index=True, right_on="sample")

# ===============================
# Helper: detect columns
# ===============================
def detect_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise RuntimeError(f"None of {candidates} found in {df.columns.tolist()}")

# ===============================
# Load DEG
# ===============================
deg1 = pd.read_csv(DEG1)
deg2 = pd.read_csv(DEG2)

padj1 = detect_col(deg1, ["padj", "FDR", "adj.P.Val"])
lfc1  = detect_col(deg1, ["log2FoldChange", "log2FC", "logFC", "avg_log2FC"])

padj2 = detect_col(deg2, ["padj", "FDR", "adj.P.Val"])
lfc2  = detect_col(deg2, ["log2FoldChange", "log2FC", "logFC", "avg_log2FC"])

# ===============================
# Define branch-specific genes
# ===============================
up1 = set(deg1.loc[(deg1[padj1] < 0.05) & (deg1[lfc1] > 0), "gene"])
up2 = set(deg2.loc[(deg2[padj2] < 0.05) & (deg2[lfc2] > 0), "gene"])

branch1 = up1 - up2
branch2 = up2 - up1

print(f"[A2] Branch1 genes: {len(branch1)}")
print(f"[A2] Branch2 genes: {len(branch2)}")

pd.DataFrame({"gene": list(branch1)}).to_csv(
    os.path.join(OUT_DIR, "A2_branch1_genes.csv"), index=False
)
pd.DataFrame({"gene": list(branch2)}).to_csv(
    os.path.join(OUT_DIR, "A2_branch2_genes.csv"), index=False
)

# ===============================
# Compute cluster means
# ===============================
means = (
    expr.groupby("cluster")
        .mean(numeric_only=True)
        .loc[ORDER]
)

def zscore(v):
    return (v - v.mean()) / v.std()

b1 = means[list(branch1 & set(means.columns))].mean(axis=1)
b2 = means[list(branch2 & set(means.columns))].mean(axis=1)

b1z = zscore(b1)
b2z = zscore(b2)

# ===============================
# Plot branching trend
# ===============================
x = np.arange(len(ORDER))

plt.figure(figsize=(6,4))
plt.plot(x, b1z, marker="o", label="Cluster 1 branch")
plt.plot(x, b2z, marker="o", label="Cluster 2 branch")
plt.xticks(x, ORDER)
plt.axvline(1.5, linestyle="--", color="gray", alpha=0.6)
plt.ylabel("Mean Z-score")
plt.title("A2 | Divergent branch signatures after cluster 3")
plt.legend()
plt.tight_layout()

plt.savefig(
    os.path.join(OUT_DIR, "A2_branching_signature_trend.png"),
    dpi=300
)
plt.close()

print("[A2] Branching analysis completed.")
print("Results in:", OUT_DIR)
