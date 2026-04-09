#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

# make src importable
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.cps_scoring import (  # noqa: E402
    CPSConfig,
    align_expression_and_pheno,
    build_run_log,
    compute_signature_scores,
    get_present_genes,
    load_expression_matrix,
    load_pheno_table,
    summarize_by_group,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run CPS signature scoring.")
    p.add_argument(
        "--config",
        default=str(PROJECT_ROOT / "config" / "cps_signature.yaml"),
        help="Path to YAML config.",
    )
    return p.parse_args()


def load_config(path: str) -> CPSConfig:
    with open(path, "r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)
    return CPSConfig(**raw)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def make_boxplot(scores: pd.DataFrame, group_column: str, output_png: Path) -> None:
    groups = list(scores[group_column].dropna().astype(str).unique())
    groups = sorted(groups)

    data = [
        scores.loc[scores[group_column].astype(str) == g, "CPS"].dropna().values
        for g in groups
    ]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.boxplot(data, labels=groups, showfliers=True)
    ax.set_title("CPS by subgroup")
    ax.set_xlabel(group_column)
    ax.set_ylabel("CPS")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    fig.savefig(output_png, dpi=200)
    plt.close(fig)


def make_heatmap(scores: pd.DataFrame, group_column: str, output_png: Path) -> None:
    tmp = scores.copy()
    tmp[group_column] = tmp[group_column].astype(str)
    tmp = tmp.sort_values([group_column, "CPS"], ascending=[True, False]).reset_index(drop=True)

    mat = tmp[["ES_positive", "ES_control", "CPS", "CPS_zscore"]].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(12, max(4, 0.28 * len(tmp))))
    im = ax.imshow(mat, aspect="auto")

    ax.set_title("Sample-level CPS signature heatmap")
    ax.set_xlabel("Metrics")
    ax.set_ylabel("Samples")
    ax.set_xticks(range(4))
    ax.set_xticklabels(["ES_positive", "ES_control", "CPS", "CPS_zscore"], rotation=30, ha="right")
    ax.set_yticks(range(len(tmp)))
    ax.set_yticklabels(
        [f"{sid} | {grp}" for sid, grp in zip(tmp["sample_id"], tmp[group_column])],
        fontsize=7,
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("score")

    plt.tight_layout()
    fig.savefig(output_png, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    output_dir = PROJECT_ROOT / cfg.output_dir
    ensure_dir(output_dir)

    expr = load_expression_matrix(
        expression_file=str(PROJECT_ROOT / cfg.expression_file),
        gene_column_candidates=cfg.gene_column_candidates,
    )
    pheno = load_pheno_table(
        pheno_file=str(PROJECT_ROOT / cfg.pheno_file),
        sample_id_column=cfg.sample_id_column,
        group_column=cfg.group_column,
    )
    expr, pheno = align_expression_and_pheno(
        expr=expr,
        pheno=pheno,
        sample_id_column=cfg.sample_id_column,
    )

    positive_present = get_present_genes(expr, cfg.positive_genes)
    control_present = get_present_genes(expr, cfg.control_genes)

    if len(positive_present) < cfg.min_genes_required:
        raise ValueError(
            f"Too few positive genes present: {len(positive_present)} < {cfg.min_genes_required}"
        )
    if len(control_present) < min(cfg.min_genes_required, len(cfg.control_genes)):
        raise ValueError(
            f"Too few control genes present: {len(control_present)} < "
            f"{min(cfg.min_genes_required, len(cfg.control_genes))}"
        )

    scores = compute_signature_scores(
        expr=expr,
        positive_genes=positive_present,
        control_genes=control_present,
        method=cfg.scoring_method,
    )

    scores = scores.merge(
        pheno[[cfg.sample_id_column, cfg.group_column]],
        left_on="sample_id",
        right_on=cfg.sample_id_column,
        how="left",
    )

    if cfg.sample_id_column != "sample_id" and cfg.sample_id_column in scores.columns:
        scores = scores.drop(columns=[cfg.sample_id_column])

    scores["CPS_class"] = np.where(
        scores["CPS_zscore"] > cfg.zscore_threshold_high,
        "CPS-High",
        "CPS-LowOrMid",
    )

    summary = summarize_by_group(
        scores=scores,
        group_column=cfg.group_column,
        zscore_threshold_high=cfg.zscore_threshold_high,
    )

    sample_scores_file = output_dir / "sample_scores.csv"
    group_summary_file = output_dir / "group_summary.csv"
    run_log_file = output_dir / "run_log.txt"
    boxplot_file = output_dir / "cps_boxplot.png"
    heatmap_file = output_dir / "cps_heatmap.png"

    scores.to_csv(sample_scores_file, index=False)
    summary.to_csv(group_summary_file, index=False)

    run_log = build_run_log(
        cfg=cfg,
        expr=expr,
        positive_present=positive_present,
        control_present=control_present,
    )
    run_log_file.write_text(run_log, encoding="utf-8")

    make_boxplot(scores=scores, group_column=cfg.group_column, output_png=boxplot_file)
    make_heatmap(scores=scores, group_column=cfg.group_column, output_png=heatmap_file)

    print("Done.")
    print(f"Output dir: {output_dir}")
    print(f"Sample scores: {sample_scores_file}")
    print(f"Group summary: {group_summary_file}")
    print(f"Boxplot: {boxplot_file}")
    print(f"Heatmap: {heatmap_file}")


if __name__ == "__main__":
    main()