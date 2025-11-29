import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN

# =========================
# 选择器
# =========================
# 可选: "RNA" 或 "smallRNA"
DATA_TYPE = "RNA"

# 可选: "dbscan" 或 "hdbscan"
METHOD = "dbscan"


# =========================
# 1. 根据开关选择数据路径
# =========================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

if DATA_TYPE == "smallRNA":
    file_path = os.path.join(
        BASE_DIR,
        "data",
        "small RNA-seq",
        "GSE254878_smallRNAs_raw_counts_expression.csv"
    )
elif DATA_TYPE == "RNA":
    file_path = os.path.join(
        BASE_DIR,
        "data",
        "RNA-seq",
        "GSE254877_raw_counts_expression.csv"
    )
else:
    raise ValueError("DATA_TYPE 必须是 'smallRNA' 或 'RNA'")

print(f"Current data type: {DATA_TYPE}")
print(f"Using file: {file_path}")
print(f"Clustering method: {METHOD}")


# =========================
# 2. 读取并清洗数据
# =========================
df = pd.read_csv(file_path, low_memory=False)

# 删除第一行异常描述行
df = df.drop(index=0)

# 自动使用第一列作为基因ID
gene_col = df.columns[0]
print("Using gene ID column:", gene_col)
df = df.set_index(gene_col)

# smallRNA 才有 type 列
if "type" in df.columns:
    df = df.drop(columns=["type"])

# 强制转为数值
df = df.apply(pd.to_numeric, errors="coerce")

# 只删除整行全空
df = df.dropna(how="all")

# 剩余 NaN 用 0 填补
df = df.fillna(0)

# 转置：样本为行
X = df.T
print("Expression matrix shape:", X.shape)


# =========================
# 3. 标准化
# =========================
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)


# =========================
# 4. 聚类：DBSCAN / HDBSCAN
# =========================
labels = None

if METHOD.lower() == "dbscan":
    # ===== 自动扫一遍 eps，选择能分出 cluster 的参数 =====
    eps_list = [0.5, 1, 2, 3, 5, 8, 10]   # 可以根据需要再加
    best_eps = None
    best_n_clusters = -1
    best_labels = None

    for eps in eps_list:
        db = DBSCAN(eps=eps, min_samples=3)   # min_samples 放宽一点
        tmp_labels = db.fit_predict(X_scaled)

        # 计算该 eps 下的簇数（排除噪声 -1）
        uniq = set(tmp_labels)
        n_clusters_tmp = len([l for l in uniq if l != -1])
        n_noise_tmp = np.sum(tmp_labels == -1)

        print(f"eps={eps}: clusters={n_clusters_tmp}, noise={n_noise_tmp}")

        # 选出簇数最多、且至少有 1 个簇的 eps
        if n_clusters_tmp > best_n_clusters and n_clusters_tmp > 0:
            best_n_clusters = n_clusters_tmp
            best_eps = eps
            best_labels = tmp_labels

    if best_labels is None:
        # 如果所有 eps 都没分出簇，就退回最后一个 eps 的结果
        print("所有 eps 都没有形成有效簇，使用最后一次结果（可能全是噪声）")
        db = DBSCAN(eps=eps_list[-1], min_samples=3)
        best_labels = db.fit_predict(X_scaled)
        best_eps = eps_list[-1]
        uniq = set(best_labels)
        best_n_clusters = len([l for l in uniq if l != -1])

    labels = np.array(best_labels)
    print(f"最终选择 eps = {best_eps}, clusters = {best_n_clusters}")

elif METHOD.lower() == "hdbscan":
    try:
        import hdbscan
    except ImportError:
        raise ImportError(
            "需要先安装 hdbscan: pip install hdbscan "
            "（建议在 subtypegan_py 环境里安装）"
        )
    clusterer = hdbscan.HDBSCAN(min_cluster_size=5)
    labels = clusterer.fit_predict(X_scaled)

else:
    raise ValueError("METHOD 必须是 'dbscan' 或 'hdbscan'")

labels = np.array(labels)
unique_labels = set(labels)
n_clusters = len([l for l in unique_labels if l != -1])
n_noise = np.sum(labels == -1)

print(f"{METHOD.upper()} found {n_clusters} clusters")
print(f"Number of noise samples: {n_noise}")



# =========================
# 5. PCA 可视化
# =========================
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

plt.figure()
mask_core = labels != -1
mask_noise = labels == -1

# 聚类点
plt.scatter(
    X_pca[mask_core, 0],
    X_pca[mask_core, 1],
    c=labels[mask_core],
    cmap="viridis",
    label="clusters"
)

# 噪声点（如果有）
if n_noise > 0:
    plt.scatter(
        X_pca[mask_noise, 0],
        X_pca[mask_noise, 1],
        c="lightgray",
        marker="x",
        label="noise"
    )

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title(f"PCA Visualization with {METHOD.upper()} ({DATA_TYPE})")
plt.legend()
plt.show()


# =========================
# 6. 保存结果到 data 下（按类型分类 + 不同命名）
# =========================
if DATA_TYPE == "smallRNA":
    output_dir = os.path.join(BASE_DIR, "data", "small RNA-seq")
else:
    output_dir = os.path.join(BASE_DIR, "data", "RNA-seq")

os.makedirs(output_dir, exist_ok=True)

method_tag = METHOD.lower()  # dbscan / hdbscan

result = pd.DataFrame({
    "Sample_ID": X.index,
    "Cluster": labels    # 注意: -1 表示噪声
})

out_name = f"cluster_{method_tag}_{DATA_TYPE}.csv"
result_path = os.path.join(output_dir, out_name)
result.to_csv(result_path, index=False)

print(f"Cluster results saved to: {result_path}")
