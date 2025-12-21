# scripts/rna_make_plots.py
# -*- coding: utf-8 -*-

import os
import sys
from pathlib import Path

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from src.multiomics_snf_v2 import load_expression_only
from src.rna_plotting import make_rna_cluster_plots


if __name__ == "__main__":
    # 你把这里改成你实际那次 run 的目录
    run_dir = "results/clustering/RNA_only/run_20251220-193800__3713d34e52"

    make_rna_cluster_plots(
        run_dir=run_dir,
        load_expression_only_fn=load_expression_only,
        top_m=50,
        fdr=0.05,
        lfc=1.0,
        max_genes_heatmap=200
    )
