#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import os
import sys

# ===============================
# CONFIG
# ===============================
if len(sys.argv) < 2:
    print("Usage: python B_TF_ORA_from_DEG.py <RUN_DIR>")
    sys.exit(1)

RUN_DIR = sys.argv[1]
DEG_DIR = os.path.join(RUN_DIR, "degs_deseq2")
OUT_DIR = os.path.join(RUN_DIR, "pseudotime", "B_TF_ORA")
os.makedirs(OUT_DIR, exist_ok=True)

PADJ_CUTOFF = 0.05
LFC_CUTOFF  = 0.0   # 上调基因

# ===============================
# 1. Load TRRUST TF → target gene sets
# ===============================
try:
    import gseapy as gp
except ImportError:
    print("Please install gseapy: pip install gseapy")
    sys.exit(1)

trrust = gp.get_library(
    name="TRRUST_Transcription_Factors_2019",
    organism="Human"
)

# TRRUST gene set名形如：STAT3 (targets)
tf_targets = {
    tf.split()[0]: set(genes)
    for tf, genes in trrust.items()
}

print(f"Loaded {len(tf_targets)} TF regulons from TRRUST")

# ===============================
# 2. Helper: ORA for one TF
# ===============================
def fisher_ora(targets, deg_genes, universe):
    targets = targets & universe
    deg_genes = deg_genes & universe

    if len(targets) < 5 or len(deg_genes) < 10:
        return np.nan, np.nan, np.nan

    a = len(targets & deg_genes)
    b = len(deg_genes - targets)
    c = len(targets - deg_genes)
    d = len(universe - (targets | deg_genes))

    if min(a, b, c, d) < 0:
        return np.nan, np.nan, np.nan

    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return a, odds, p

# ===============================
# 3. Run ORA per cluster
# ===============================
all_results = []

for k in ["0", "1", "2", "3"]:
    deg_path = os.path.join(DEG_DIR, f"DEG_cluster_{k}_vs_rest.csv")
    if not os.path.exists(deg_path):
        print(f"Skip cluster {k}: DEG file not found")
        continue

    deg = pd.read_csv(deg_path)

    # ---------- 自动识别 DEG 列名 ----------
    col_gene = [c for c in deg.columns if c.lower() in ["gene", "genes", "symbol"]]
    col_padj = [c for c in deg.columns if ("padj" in c.lower()) or ("fdr" in c.lower()) or ("adj" in c.lower())]
    col_lfc  = [c for c in deg.columns if ("log2foldchange" in c.lower()) or ("logfc" in c.lower()) or ("lfc" in c.lower())]

    if len(col_gene) == 0 or len(col_padj) == 0 or len(col_lfc) == 0:
        raise ValueError(
            f"Cannot detect required columns in {deg_path}. "
            f"Found columns: {list(deg.columns)}"
        )

    col_gene = col_gene[0]
    col_padj = col_padj[0]
    col_lfc  = col_lfc[0]

    deg = deg.rename(columns={
        col_gene: "gene",
        col_padj: "padj",
        col_lfc: "log2FoldChange"
    })

    deg = deg.dropna(subset=["gene", "padj", "log2FoldChange"])

    # 上调 DEG
    deg_up = set(
        deg.loc[
            (deg["padj"] < PADJ_CUTOFF) &
            (deg["log2FoldChange"] > LFC_CUTOFF),
            "gene"
        ]
    )

    universe = set(deg["gene"])

    records = []
    for tf, targets in tf_targets.items():
        overlap, odds, p = fisher_ora(targets, deg_up, universe)
        records.append({
            "cluster": k,
            "TF": tf,
            "overlap": overlap,
            "odds_ratio": odds,
            "pvalue": p
        })

    res = pd.DataFrame(records)
    res["padj"] = multipletests(res["pvalue"], method="fdr_bh")[1]
    res = res.sort_values("padj")

    out_path = os.path.join(OUT_DIR, f"TF_ORA_cluster_{k}.csv")
    res.to_csv(out_path, index=False)
    all_results.append(res)

    print(f"Cluster {k}: saved {out_path}")

# ===============================
# 4. Merge all clusters
# ===============================
if len(all_results) > 0:
    df_all = pd.concat(all_results, ignore_index=True)
    df_all.to_csv(
        os.path.join(OUT_DIR, "TF_ORA_all_clusters.csv"),
        index=False
    )

print("TF ORA finished.")
print("Results in:", OUT_DIR)


