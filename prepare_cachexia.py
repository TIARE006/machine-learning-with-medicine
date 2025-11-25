import os
import re
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

# ========= 1. 路径 =========
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.join(ROOT_DIR, "data")

# 原始数据
RNA_EXPR_FILE  = os.path.join(BASE_PATH, "RNA-seq",       "GSE254877_raw_counts_expression.csv")
SM_EXPR_FILE   = os.path.join(BASE_PATH, "small RNA-seq", "GSE254878_smallRNAs_raw_counts_expression.csv")
RNA_PHENO_FILE = os.path.join(BASE_PATH, "RNA-seq",       "GSE254877_pheno.csv")
SM_PHENO_FILE  = os.path.join(BASE_PATH, "small RNA-seq", "GSE254878_pheno.csv")

# 预处理结果输出目录
PREPARED_DIR = os.path.join(BASE_PATH, "prepared data")
os.makedirs(PREPARED_DIR, exist_ok=True)

TOPK_RNA   = 2000   # 选变异度最大的基因
TOPK_SMRNA = 300    # 选变异度最大的 miRNA

print("Script dir:", ROOT_DIR)
print("Data dir:", BASE_PATH)
print("RNA expr exists:", os.path.exists(RNA_EXPR_FILE))
print("smRNA expr exists:", os.path.exists(SM_EXPR_FILE))

# ========= 2. 读取 =========
# low_memory=False + to_numeric 解决 dtype warning
rna_expr  = pd.read_csv(RNA_EXPR_FILE, index_col=0, low_memory=False)
rna_expr  = rna_expr.apply(pd.to_numeric, errors="coerce")

sm_expr   = pd.read_csv(SM_EXPR_FILE,  index_col=0, low_memory=False)
sm_expr   = sm_expr.apply(pd.to_numeric, errors="coerce")

rna_pheno = pd.read_csv(RNA_PHENO_FILE, index_col=0)
sm_pheno  = pd.read_csv(SM_PHENO_FILE,  index_col=0)

print("RNA first 5 columns:", list(rna_expr.columns[:5]))
print("smRNA first 5 columns:", list(sm_expr.columns[:5]))
print("RNA expr shape:", rna_expr.shape)
print("smRNA expr shape:", sm_expr.shape)

# ========= 3. 对齐样本（用列名里的 4 位数字ID） =========

def extract_sample_id(name: str):
    """
    从列名中提取“最后一个”4位数字作为样本ID，比如 2066, 3007, 3016
    """
    nums = re.findall(r"\d{4}", str(name))
    return nums[-1] if nums else None


# 删掉明显不是样本的列
rna_expr = rna_expr.drop(columns=["Unnamed: 1"], errors="ignore")
sm_expr  = sm_expr.drop(columns=["type"],        errors="ignore")

# 构建：样本ID -> 原列名
rna_id_map = {}
for c in rna_expr.columns:
    sid = extract_sample_id(c)
    if sid is not None and sid not in rna_id_map:
        rna_id_map[sid] = c

sm_id_map = {}
for c in sm_expr.columns:
    sid = extract_sample_id(c)
    if sid is not None and sid not in sm_id_map:
        sm_id_map[sid] = c

# 找共同样本 ID
common_ids = sorted(set(rna_id_map.keys()) & set(sm_id_map.keys()))
print("Common sample IDs (first 10):", common_ids[:10])
print("Total common samples:", len(common_ids))

if len(common_ids) == 0:
    raise RuntimeError("没有找到共同样本ID，检查列名或 extract_sample_id 规则。")

# 用共同 ID 对齐表达矩阵（先用原列名取，再重命名）
rna_expr = rna_expr[[rna_id_map[i] for i in common_ids]]
sm_expr  = sm_expr[[sm_id_map[i]  for i in common_ids]]

rna_expr.columns = common_ids
sm_expr.columns  = common_ids

# ========= 3.1 pheno 也按相同 ID 对齐（尽量） =========
def align_pheno_by_ids(pheno_df, common_ids):
    # 把索引里的样本ID提出来做一列，然后设成 index
    pheno = pheno_df.copy()
    pheno["sample_id"] = pheno.index.map(extract_sample_id)
    pheno = pheno.dropna(subset=["sample_id"]).set_index("sample_id")
    # 只保留有表达数据的样本
    keep_ids = [i for i in common_ids if i in pheno.index]
    return pheno.loc[keep_ids], keep_ids

rna_pheno_aligned, ids_with_rna_pheno = align_pheno_by_ids(rna_pheno, common_ids)
sm_pheno_aligned, ids_with_sm_pheno  = align_pheno_by_ids(sm_pheno,  common_ids)

# 只有同时在 RNA+smRNA+pheno 中都出现的样本才保留
final_ids = [i for i in common_ids
             if (i in ids_with_rna_pheno) and (i in ids_with_sm_pheno)]

print("Samples with both omics + pheno:", len(final_ids))

if len(final_ids) == 0:
    # 如果 pheno 对不齐，就只按表达数据的 common_ids 来，pheno 原样保存
    print("Warning: pheno 不能完全对齐，后面只根据表达矩阵的 common_ids 保存。")
    final_ids = common_ids
    rna_pheno_subset = rna_pheno
    sm_pheno_subset  = sm_pheno
else:
    # 按 final_ids 重新对子集
    rna_expr = rna_expr[final_ids]
    sm_expr  = sm_expr[final_ids]
    rna_pheno_subset = rna_pheno_aligned.loc[final_ids]
    sm_pheno_subset  = sm_pheno_aligned.loc[final_ids]

# ========= 4. 选变异度最大的特征 =========
def select_top_var(df, topk):
    var = df.var(axis=1)
    idx = var.sort_values(ascending=False).head(topk).index
    return df.loc[idx]

rna_expr = select_top_var(rna_expr, TOPK_RNA)
sm_expr  = select_top_var(sm_expr,  TOPK_SMRNA)

print("Filtered RNA:", rna_expr.shape)
print("Filtered smRNA:", sm_expr.shape)

# ========= 5. log1p + z-score 标准化 =========
def log_zscore(df):
    X = np.log1p(df.values.astype(float))  # log1p
    scaler = StandardScaler()
    X = scaler.fit_transform(X.T).T        # 按样本标准化 → 先转置再转回来
    return X

rna_mat = log_zscore(rna_expr).T   # (samples, features)
sm_mat  = log_zscore(sm_expr).T

samples = np.array(final_ids)

print("Final RNA matrix:", rna_mat.shape)
print("Final smRNA matrix:", sm_mat.shape)
print("Final samples:", samples.shape)

# ========= 6. 保存到 prepared data =========
np.save(os.path.join(PREPARED_DIR, "cachexia_rna.npy"),   rna_mat)
np.save(os.path.join(PREPARED_DIR, "cachexia_smrna.npy"), sm_mat)
np.save(os.path.join(PREPARED_DIR, "cachexia_samples.npy"), samples)

rna_expr.index.to_series().to_csv(
    os.path.join(PREPARED_DIR, "cachexia_rna_genes.txt"),
    index=False
)
sm_expr.index.to_series().to_csv(
    os.path.join(PREPARED_DIR, "cachexia_smrna_ids.txt"),
    index=False
)

rna_pheno_subset.to_csv(
    os.path.join(PREPARED_DIR, "cachexia_rna_pheno_subset.csv")
)
sm_pheno_subset.to_csv(
    os.path.join(PREPARED_DIR, "cachexia_smrna_pheno_subset.csv")
)

print("Done preprocessing.")
print("RNA pheno columns:", rna_pheno_subset.columns.tolist())
print("smRNA pheno columns:", sm_pheno_subset.columns.tolist())

