# src/pipelines/rna_spectral_pipeline.py
# -*- coding: utf-8 -*-

import os
import sys
import json
import time
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

# 当前文件: scripts/rna_spectral_pipeline.py
# 项目根目录: machine learning with medicine/
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

# ============================================================
# Project imports (methods library)
# ============================================================
from src.methods_utils import (
    BASE_DIR,
    RANDOM_STATE,
)

from src.multiomics_snf_v2 import (
    load_expression_only,
    select_top_variable_features,
    construct_similarity_matrix,
    snf_clustering,
)


def _make_run_id(cfg: dict) -> str:
    """
    为本次运行生成一个稳定的短 hash（只由关键参数决定）
    """
    key = json.dumps(cfg, sort_keys=True).encode("utf-8")
    return hashlib.sha1(key).hexdigest()[:10]


def run_rna_spectral_pipeline(
    top_n: int = 2000,
    knn_k: int = 20,
    mu: float = 0.5,
    candidate_K=(3, 4, 5),
    seed: int = RANDOM_STATE,
    # 输出根目录（推荐 results，而不是 data）
    results_root: str | None = None,
    # 是否保存 affinity 矩阵（可能很大，默认不存）
    save_affinity: bool = False,
):
    """
    RNA-only: RNA -> similarity graph S -> spectral clustering

    输出结构（干净、可复现）：
    results/clustering/RNA_only/
      run_<timestamp>__<hash>/
        labels/cluster_labels.csv
        metrics/summary.txt
        params/config.json
        (optional) affinity/S_rna.npy
    """
    print("===== RNA-only Spectral Clustering (Graph-based) =====")

    # ---------- 0) Resolve output root ----------
    if results_root is None:
        # BASE_DIR 通常指向项目根；若 BASE_DIR 不是项目根，可按需调整
        results_root = os.path.join(BASE_DIR, "results")

    cfg = {
        "method": "RNA_only_spectral",
        "top_n": int(top_n),
        "knn_k": int(knn_k),
        "mu": float(mu),
        "candidate_K": list(candidate_K),
        "seed": int(seed),
        "save_affinity": bool(save_affinity),
    }
    run_hash = _make_run_id(cfg)
    ts = time.strftime("%Y%m%d-%H%M%S")
    run_name = f"run_{ts}__{run_hash}"

    run_dir = Path(results_root) / "clustering" / "RNA_only" / run_name
    labels_dir = run_dir / "labels"
    metrics_dir = run_dir / "metrics"
    params_dir = run_dir / "params"
    affinity_dir = run_dir / "affinity"

    for d in [labels_dir, metrics_dir, params_dir]:
        d.mkdir(parents=True, exist_ok=True)
    if save_affinity:
        affinity_dir.mkdir(parents=True, exist_ok=True)

    # ---------- 1) Load RNA expression ----------
    X_rna = load_expression_only("RNA")

    # -------- 1.1) Remove fake samples like 'Unnamed:*' --------
    X_rna = X_rna.loc[
        ~X_rna.index.astype(str).str.match(r"^Unnamed")
    ]
    X_rna.index = X_rna.index.astype(str).str.strip()

    # ---------- 2) Feature selection ----------
    X_rna_var = select_top_variable_features(X_rna, top_n=top_n)

    # ---------- 3) Build similarity graph S ----------
    S_rna = construct_similarity_matrix(X_rna_var.values, K=knn_k, mu=mu)

    if save_affinity:
        np.save(str(affinity_dir / "S_rna.npy"), S_rna)

    # ---------- 4) Spectral clustering + choose K ----------
    labels, best_k, best_score = snf_clustering(
        S_rna, candidate_K=candidate_K, random_state=seed
    )

    # ---------- 5) Save labels ----------
    out_csv = labels_dir / "cluster_labels.csv"
    pd.DataFrame({
        "Sample_ID": X_rna_var.index.astype(str),
        "Cluster": labels
    }).to_csv(out_csv, index=False)

    # -------- 5.1) Save cluster sizes --------
    cluster_sizes = (
        pd.Series(labels, name="Cluster")
        .value_counts()
        .sort_index()
        .reset_index()
    )
    cluster_sizes.columns = ["Cluster", "N_samples"]

    cluster_sizes.to_csv(
        metrics_dir / "cluster_sizes.csv",
        index=False
    )


    # ---------- 6) Save params ----------
    cfg_out = cfg.copy()
    cfg_out.update({
        "n_samples": int(X_rna_var.shape[0]),
        "n_features_used": int(X_rna_var.shape[1]),
        "best_k": int(best_k) if best_k is not None else None,
        "best_silhouette_embedding": float(best_score),
        "run_name": run_name,
    })
    with open(params_dir / "config.json", "w", encoding="utf-8") as f:
        json.dump(cfg_out, f, ensure_ascii=False, indent=2)

    # ---------- 7) Save metrics summary (human-readable) ----------
    summary_txt = metrics_dir / "summary.txt"
    with open(summary_txt, "w", encoding="utf-8") as f:
        f.write("RNA-only Spectral Clustering (Graph-based)\n")
        f.write(f"run: {run_name}\n")
        f.write(f"samples: {X_rna_var.shape[0]}, features_used: {X_rna_var.shape[1]}\n")
        f.write(f"top_n={top_n}, knn_k={knn_k}, mu={mu}, candidate_K={tuple(candidate_K)}, seed={seed}\n")
        f.write(f"Best K = {best_k}, silhouette(embedding) = {best_score:.4f}\n")
        f.write(f"labels_csv: {out_csv}\n")
        if save_affinity:
            f.write(f"affinity_npy: {affinity_dir / 'S_rna.npy'}\n")

    # ---------- 8) Print concise outputs ----------
    print(f"[OUT] {run_dir}")
    print(f"[SAVE] labels -> {out_csv}")
    print(f"[SAVE] metrics -> {summary_txt}")
    print(f"[SAVE] params -> {params_dir / 'config.json'}")
    if save_affinity:
        print(f"[SAVE] affinity -> {affinity_dir / 'S_rna.npy'}")
    print(f"[DONE] Best K={best_k}, silhouette(embedding)={best_score:.4f}")

    return labels, best_k, best_score, str(run_dir)


if __name__ == "__main__":
    run_rna_spectral_pipeline()
