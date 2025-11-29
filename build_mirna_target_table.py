import os
import gzip
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

in_file = os.path.join(
    BASE_DIR,
    "data", "reference", "gene_attribute_edges.txt.gz"
)
out_file = os.path.join(
    BASE_DIR,
    "data", "reference", "mirna_target_human.csv"
)

print("[INFO] Loading:", in_file)

# 1. 读取 gzip 文本
with gzip.open(in_file, "rt") as f:
    df = pd.read_csv(f, sep="\t")

print("[INFO] Raw table shape:", df.shape)
print(df.head())

# 2. 删除第一行“说明行”（source = 'GeneSym' 的那一行）
df = df[df["source"] != "GeneSym"]

# 3. 把列映射成我们需要的格式：
#    source  -> TargetGene
#    target  -> miRNA
df = df.rename(columns={
    "source": "TargetGene",
    "target": "miRNA"
})

# 4. 只保留这两列，并去重
df2 = df[["miRNA", "TargetGene"]].dropna().drop_duplicates()

print("[INFO] Processed table shape:", df2.shape)
print(df2.head())

# 5. 保存成 CSV，给 pipeline 用
os.makedirs(os.path.dirname(out_file), exist_ok=True)
df2.to_csv(out_file, index=False)

print("[SAVE] miRNA–target mapping saved to:")
print("      ", out_file)
