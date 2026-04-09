from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Dict, Tuple

import numpy as np
import pandas as pd


@dataclass
class CPSConfig:
    dataset: str
    expression_file: str
    pheno_file: str
    sample_id_column: str
    group_column: str
    gene_column_candidates: List[str]
    scoring_method: str
    positive_genes: List[str]
    control_genes: List[str]
    output_dir: str
    min_genes_required: int = 2
    zscore_threshold_high: float = 0.5


def _normalize_gene_name(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    s = s.replace('"', "").replace("'", "")
    return s.upper()


def _safe_zscore_1d(values: pd.Series) -> pd.Series:
    arr = values.astype(float)
    mu = arr.mean()
    sigma = arr.std(ddof=0)
    if sigma == 0 or np.isnan(sigma):
        return pd.Series(np.zeros(len(arr)), index=values.index)
    return (arr - mu) / sigma


def _safe_global_zscore(values: pd.Series) -> pd.Series:
    arr = values.astype(float)
    mu = arr.mean()
    sigma = arr.std(ddof=0)
    if sigma == 0 or np.isnan(sigma):
        return pd.Series(np.zeros(len(arr)), index=values.index)
    return (arr - mu) / sigma


def load_expression_matrix(
    expression_file: str,
    gene_column_candidates: Iterable[str],
) -> pd.DataFrame:
    """
    Returns expression dataframe:
      index   = gene symbols (upper-case)
      columns = sample IDs
      values  = numeric expression
    Supports a first gene column or an explicit named gene column.
    """
    path = Path(expression_file)
    if not path.exists():
        raise FileNotFoundError(f"Expression file not found: {expression_file}")

    df = pd.read_csv(path)

    gene_col = None
    for c in gene_column_candidates:
        if c in df.columns:
            gene_col = c
            break

    if gene_col is None:
        gene_col = df.columns[0]

    df = df.copy()
    df[gene_col] = df[gene_col].map(_normalize_gene_name)
    df = df[df[gene_col] != ""].copy()

    numeric_cols = [c for c in df.columns if c != gene_col]
    if not numeric_cols:
        raise ValueError("No sample columns found in expression matrix.")

    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(axis=0, how="all", subset=numeric_cols)
    df = df.groupby(gene_col, as_index=True)[numeric_cols].mean()
    return df


def load_pheno_table(pheno_file: str, sample_id_column: str, group_column: str) -> pd.DataFrame:
    path = Path(pheno_file)
    if not path.exists():
        raise FileNotFoundError(f"Pheno file not found: {pheno_file}")

    pheno = pd.read_csv(path)
    missing = [c for c in [sample_id_column, group_column] if c not in pheno.columns]
    if missing:
        raise ValueError(f"Missing required columns in pheno file: {missing}")

    pheno = pheno.copy()
    pheno[sample_id_column] = pheno[sample_id_column].astype(str)
    pheno[group_column] = pheno[group_column].astype(str)
    return pheno


def align_expression_and_pheno(
    expr: pd.DataFrame,
    pheno: pd.DataFrame,
    sample_id_column: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sample_ids = pheno[sample_id_column].astype(str).tolist()
    common = [s for s in sample_ids if s in expr.columns]

    if not common:
        raise ValueError(
            "No overlapping sample IDs between expression columns and pheno sample_id_column."
        )

    expr2 = expr.loc[:, common].copy()
    pheno2 = pheno[pheno[sample_id_column].isin(common)].copy()
    pheno2 = pheno2.set_index(sample_id_column).loc[common].reset_index()
    return expr2, pheno2


def get_present_genes(expr: pd.DataFrame, genes: Iterable[str]) -> List[str]:
    norm = [_normalize_gene_name(g) for g in genes]
    return [g for g in norm if g in expr.index]


def score_mean_z(
    expr: pd.DataFrame,
    genes: List[str],
) -> pd.Series:
    """
    For selected genes:
    1) z-score each gene across samples
    2) average across genes for each sample
    """
    if len(genes) == 0:
        return pd.Series(np.nan, index=expr.columns)

    sub = expr.loc[genes].copy()
    z = sub.apply(_safe_zscore_1d, axis=1)
    return z.mean(axis=0)


def score_ssgsea_like(
    expr: pd.DataFrame,
    genes: List[str],
) -> pd.Series:
    """
    A simple ssGSEA-like rank-based score:
    - rank genes within each sample
    - normalize ranks to [0,1]
    - average normalized ranks of selected genes
    This is not identical to GSVA::gsva(method='ssgsea'),
    but gives a stable within-sample rank-based enrichment signal.
    """
    if len(genes) == 0:
        return pd.Series(np.nan, index=expr.columns)

    genes = [g for g in genes if g in expr.index]
    if not genes:
        return pd.Series(np.nan, index=expr.columns)

    ranks = expr.rank(axis=0, method="average", ascending=True)
    denom = max(len(expr.index) - 1, 1)
    ranks = (ranks - 1) / denom
    return ranks.loc[genes].mean(axis=0)


def compute_signature_scores(
    expr: pd.DataFrame,
    positive_genes: List[str],
    control_genes: List[str],
    method: str = "mean_z",
) -> pd.DataFrame:
    if method not in {"mean_z", "ssgsea_like"}:
        raise ValueError(f"Unsupported scoring method: {method}")

    if method == "mean_z":
        pos = score_mean_z(expr, positive_genes)
        ctrl = score_mean_z(expr, control_genes)
    else:
        pos = score_ssgsea_like(expr, positive_genes)
        ctrl = score_ssgsea_like(expr, control_genes)

    out = pd.DataFrame(
        {
            "sample_id": expr.columns.astype(str),
            "ES_positive": pos.values,
            "ES_control": ctrl.values,
        }
    )
    out["CPS"] = out["ES_positive"] - out["ES_control"]
    out["CPS_zscore"] = _safe_global_zscore(out["CPS"])
    return out


def summarize_by_group(
    scores: pd.DataFrame,
    group_column: str,
    zscore_threshold_high: float,
) -> pd.DataFrame:
    tmp = scores.copy()
    tmp["CPS_class"] = np.where(
        tmp["CPS_zscore"] > zscore_threshold_high,
        "CPS-High",
        "CPS-LowOrMid",
    )

    summary = (
        tmp.groupby(group_column, dropna=False)
        .agg(
            n=("sample_id", "count"),
            mean_ES_positive=("ES_positive", "mean"),
            mean_ES_control=("ES_control", "mean"),
            mean_CPS=("CPS", "mean"),
            median_CPS=("CPS", "median"),
            std_CPS=("CPS", "std"),
            mean_CPS_zscore=("CPS_zscore", "mean"),
            n_CPS_high=("CPS_class", lambda s: int((s == "CPS-High").sum())),
        )
        .reset_index()
    )
    summary["frac_CPS_high"] = summary["n_CPS_high"] / summary["n"].replace(0, np.nan)
    return summary


def build_run_log(
    cfg: CPSConfig,
    expr: pd.DataFrame,
    positive_present: List[str],
    control_present: List[str],
) -> str:
    lines = [
        f"dataset: {cfg.dataset}",
        f"expression_file: {cfg.expression_file}",
        f"pheno_file: {cfg.pheno_file}",
        f"scoring_method: {cfg.scoring_method}",
        f"n_genes_total: {expr.shape[0]}",
        f"n_samples_total: {expr.shape[1]}",
        "",
        f"positive_genes_requested ({len(cfg.positive_genes)}): {', '.join(cfg.positive_genes)}",
        f"positive_genes_present   ({len(positive_present)}): {', '.join(positive_present)}",
        "",
        f"control_genes_requested  ({len(cfg.control_genes)}): {', '.join(cfg.control_genes)}",
        f"control_genes_present    ({len(control_present)}): {', '.join(control_present)}",
        "",
        f"min_genes_required: {cfg.min_genes_required}",
        f"zscore_threshold_high: {cfg.zscore_threshold_high}",
    ]
    return "\n".join(lines)