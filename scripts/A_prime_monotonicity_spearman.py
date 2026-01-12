#!/usr/bin/env python3
import os
import sys
import glob
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

# ============================================================
# 0. Parse arguments
# ============================================================
if len(sys.argv) < 2:
    print("Usage: python A_prime_monotonicity_spearman.py <RUN_DIR>")
    sys.exit(1)

RUN_DIR = sys.argv[1]

# ============================================================
# 1. Resolve paths (WRITE INTO pseudotime/A_prime)
# ============================================================
VST_PATH = os.path.join(RUN_DIR, "degs_deseq2", "vst_matrix.csv")
LABEL_DIR = os.path.join(RUN_DIR, "labels")
TOP20_GLOB = os.path.join(RUN_DIR, "heatmaps", "top_markers_TOP20*.csv")

label_files = glob.glob(os.path.join(LABEL_DIR, "cluster_results_*.csv"))
if len(label_files) == 0:
    raise FileNotFoundError("No cluster_results_*.csv found in labels/")

CLUSTER_PATH = label_files[0]
TOP20_PATHS = glob.glob(TOP20_GLOB)

# ---- OUTPUT DIR: pseudotime/A_prime ----
PSEUDO_DIR = os.path.join(RUN_DIR, "pseudotime")
OUT_DIR = os.path.join(PSEUDO_DIR, "A_prime")
os.makedirs(OUT_DIR, exist_ok=True)

OUT_PREFIX = os.path.join(OUT_DIR, "A_prime_monotonicity")

# ============================================================
# 2. Fixed biological configuration
# ============================================================
ORDER = ["0", "3", "1", "2"]
TIME = np.arange(len(ORDER))  # [0,1,2,3]

MODULES = {
    "Metabolism": ["HMGCS2", "CPT1A", "CD36", "PPARA", "ACOX1"],
    "Inflammation": ["LIF", "STAT3", "SOCS3", "IL6ST", "NFKBIA"],
    "ECM_Mechanical": ["MMP14", "COL1A1", "CTGF", "PIEZO1", "PIEZO2"]
}

# ============================================================
# 3. Load data
# ============================================================
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

# ============================================================
# 4. Cluster mean expression
# ============================================================
expr_long = vst.T.merge(clusters, left_index=True, right_on="sample")

cluster_means = (
    expr_long
    .groupby("cluster")
    .mean(numeric_only=True)
    .loc[ORDER]
)

# ============================================================
# 5. Gene-wise Spearman monotonicity
# ============================================================
records = []

for gene in cluster_means.columns:
    y = cluster_means[gene].values
    if np.any(np.isnan(y)):
        continue
    rho, p = spearmanr(y, TIME)
    records.append({
        "gene": gene,
        "spearman_rho": rho,
        "pvalue": p
    })

res = pd.DataFrame(records)
res["padj"] = multipletests(res["pvalue"], method="fdr_bh")[1]
res.to_csv(f"{OUT_PREFIX}_all_genes.csv", index=False)

sig_all = (res["padj"] < 0.05).mean()

# ============================================================
# 6. Top20 marker support (if available)
# ============================================================
sig_top20 = np.nan

if len(TOP20_PATHS) > 0:
    top20 = pd.read_csv(TOP20_PATHS[0])
    if "gene" in top20.columns:
        top20_genes = set(top20["gene"])
        res_top20 = res[res["gene"].isin(top20_genes)]
        sig_top20 = (res_top20["padj"] < 0.05).mean()
        res_top20.to_csv(f"{OUT_PREFIX}_top20.csv", index=False)

# ============================================================
# 7. Module-level monotonicity
# ============================================================
module_records = []

for module, genes in MODULES.items():
    genes = [g for g in genes if g in cluster_means.columns]
    if len(genes) < 3:
        continue
    module_score = cluster_means[genes].mean(axis=1).values
    rho, p = spearmanr(module_score, TIME)
    module_records.append({
        "module": module,
        "spearman_rho": rho,
        "pvalue": p
    })

mod = pd.DataFrame(module_records)
mod["padj"] = multipletests(mod["pvalue"], method="fdr_bh")[1]
mod.to_csv(f"{OUT_PREFIX}_modules.csv", index=False)

# ============================================================
# 8. Summary for manuscript
# ============================================================
summary_df = pd.DataFrame([{
    "all_genes_monotonic_fraction": sig_all,
    "top20_monotonic_fraction": sig_top20
}])
summary_df.to_csv(f"{OUT_PREFIX}_summary.csv", index=False)

# ============================================================
# 9. Console report
# ============================================================
print("\n=== A′ Monotonicity Summary (pseudotime) ===")
print(f"RUN_DIR: {RUN_DIR}")
print(f"All genes (FDR<0.05): {sig_all:.3f}")
if not np.isnan(sig_top20):
    print(f"Top20 markers (FDR<0.05): {sig_top20:.3f}")
else:
    print("Top20 markers: NA (file not found)")
print("\nModule-level results:")
print(mod)
print("\nResults written to:")
print(OUT_DIR)
