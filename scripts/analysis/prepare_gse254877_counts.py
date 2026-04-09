#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

project_root = Path(".").resolve()

infile = project_root / "data/raw/RNA_seq/GSE254877_raw_counts_expression.csv"
outfile = project_root / "data/raw/RNA_seq/GSE254877_counts_matrix_for_cps.csv"

# 读原始文件，不让 pandas 自动猜 header
raw = pd.read_csv(infile, header=None)

# 第0行：bam 文件名
# 第1行：真正列名（FeatureID, type, AB_3, AB_4, ...）
header_row = raw.iloc[1].tolist()

# 数据从第2行开始
df = raw.iloc[2:].copy()
df.columns = header_row

# 保留基因名列
if "FeatureID" not in df.columns:
    raise ValueError("Expected column 'FeatureID' not found.")

# 去掉 type 列，只保留 FeatureID + 样本列
keep_cols = ["FeatureID"] + [c for c in df.columns if c not in ["FeatureID", "type"]]
df = df[keep_cols].copy()

# 数值化样本列
sample_cols = [c for c in df.columns if c != "FeatureID"]
for c in sample_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# 去掉空基因名
df["FeatureID"] = df["FeatureID"].astype(str).str.strip()
df = df[df["FeatureID"] != ""].copy()

# 同名基因聚合
df = df.groupby("FeatureID", as_index=False)[sample_cols].mean()

df.to_csv(outfile, index=False)
print(f"written: {outfile}")
print(df.head())
print(f"shape: {df.shape}")