#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import os, sys

if len(sys.argv) < 2:
    print("Usage: python B_final_core_TF_selection.py <RUN_DIR>")
    sys.exit(1)

RUN_DIR = sys.argv[1]
OUT_PATH = os.path.join(RUN_DIR, "pseudotime", "Core_TF_candidates.csv")

VST_PATH = os.path.join(RUN_DIR, "degs_deseq2", "vst_matrix.csv")
CLUSTER_PATH = os.path.join(RUN_DIR, "labels", "cluster_results_RNA_seed42.csv")

ORDER = ["0", "3", "1", "2"]
TIME = np.arange(len(ORDER))

# ===============================
# 1. Load expression + cluster
# ===============================
vst = pd.read_csv(VST_PATH, index_col=0)
clusters = pd.read_csv(CLUSTER_PATH)

sample_col = [c for c in clusters.columns if "sample" in c.lower() or "id" in c.lower()][0]
cluster_col = [c for c in clusters.columns if "cluster" in c.lower()][0]

clusters = clusters[[sample_col, cluster_col]].rename(
    columns={sample_col: "sample", cluster_col: "cluster"}
)
clusters["cluster"] = clusters["cluster"].astype(str)

common = vst.columns.intersection(clusters["sample"])
vst = vst[common]
clusters = clusters[clusters["sample"].isin(common)]

expr_long = vst.T.merge(clusters, left_index=True, right_on="sample")
cluster_means = expr_long.groupby("cluster").mean(numeric_only=True).loc[ORDER]

# ===============================
# 2. Load TRRUST regulons
# ===============================
import gseapy as gp

trrust = gp.get_library(
    name="TRRUST_Transcription_Factors_2019",
    organism="Human"
)

tf_targets = {
    tf.split()[0]: set(genes)
    for tf, genes in trrust.items()
}

print(f"Loaded {len(tf_targets)} TF regulons from TRRUST")

# ===============================
# 3. Expression-driven TF scoring
# ===============================
records = []

for tf, targets in tf_targets.items():
    targets = [g for g in targets if g in cluster_means.columns]
    if len(targets) < 2:
        continue

    # target module score
    target_score = cluster_means[targets].mean(axis=1)
    target_z = (target_score - target_score.mean()) / target_score.std()

    # TF expression trend
    if tf in cluster_means.columns:
        tf_expr = cluster_means[tf].values
        rho, p = spearmanr(tf_expr, TIME)
    else:
        rho, p = np.nan, np.nan

    for i, k in enumerate(ORDER):
        records.append({
            "cluster": k,
            "TF": tf,
            "target_module_score": float(target_z.loc[k]),
            "tf_expr": cluster_means.loc[k, tf] if tf in cluster_means.columns else np.nan,
            "tf_expr_spearman": rho,
            "tf_expr_pvalue": p
        })

df = pd.DataFrame(records)

# ===============================
# 4. Candidate selection (robust)
# ===============================
df["candidate"] = (
    (df["target_module_score"] > 0) |
    (df["tf_expr_spearman"].abs() > 0.3)
)

df = df[df["candidate"]].copy()

df = df.sort_values(
    ["cluster", "target_module_score", "tf_expr_spearman"],
    ascending=[True, False, False]
)

# 每个 cluster 保留 Top 5
df = df.groupby("cluster").head(5).reset_index(drop=True)

df.to_csv(OUT_PATH, index=False)

print("Core TF candidates written to:")
print(OUT_PATH)

