import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


# =========================
# 聚类数据类型开关 
# 可选： "smallRNA" 或 "RNA"
# =========================
DATA_TYPE = "RNA"   # "smallRNA" 或 "RNA"

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)


# =========================
# 1. 路径设置：raw 输入 & clustering 输出
# =========================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

if DATA_TYPE == "smallRNA":
    expr_file = os.path.join(
        BASE_DIR,
        "data",
        "small RNA-seq",
        "raw",
        "GSE254878_smallRNAs_raw_counts_expression.csv"
    )
    out_base = os.path.join(BASE_DIR, "data", "small RNA-seq")
elif DATA_TYPE == "RNA":
    expr_file = os.path.join(
        BASE_DIR,
        "data",
        "RNA-seq",
        "raw",
        "GSE254877_raw_counts_expression.csv"
    )
    out_base = os.path.join(BASE_DIR, "data", "RNA-seq")
else:
    raise ValueError("DATA_TYPE 必须是 'smallRNA' 或 'RNA'")

# 聚类结果和图统一放在 clustering 子目录
clustering_dir = os.path.join(out_base, "clustering")
os.makedirs(clustering_dir, exist_ok=True)

print(f"Current clustering mode: {DATA_TYPE}")
print(f"Using file: {expr_file}")


# =========================
# 2. 读取并清洗数据
# =========================
df = pd.read_csv(expr_file, low_memory=False)

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

# 只删除“整行全空”的情况
df = df.dropna(how="all")

# 用 0 填补剩余 NaN（表达为 0 合理）
df = df.fillna(0)

# 转置：样本作为行
X = df.T

print("Expression matrix shape:", X.shape)


# =========================
# 3. 标准化
# =========================
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)


# =========================
# 4. 共识聚类
# =========================
def consensus_clustering(X, k_min=2, k_max=7, n_iter=50, subsample_rate=0.8):
    """
    简单版共识聚类：
    - 对每个 K 进行多次子采样 + KMeans 聚类
    - 统计样本两两“被分到同一类”的频率
    - 计算每个 K 的平均共识得分（越高越稳定）
    """
    n_samples = X.shape[0]
    consensus_scores = {}

    for k in range(k_min, k_max + 1):
        consensus_matrix = np.zeros((n_samples, n_samples), dtype=float)

        for _ in range(n_iter):
            # 随机抽一部分样本做子聚类
            idx = np.random.choice(
                n_samples,
                int(subsample_rate * n_samples),
                replace=False
            )
            X_sub = X[idx]

            labels = KMeans(
                n_clusters=k,
                random_state=None
            ).fit_predict(X_sub)

            # 在子样本内部更新共识矩阵
            for i in range(len(idx)):
                for j in range(len(idx)):
                    if labels[i] == labels[j]:
                        consensus_matrix[idx[i], idx[j]] += 1

        # 归一化到 [0,1]
        consensus_matrix /= n_iter
        # 计算整个矩阵的平均共识作为该 K 的稳定性指标
        consensus_scores[k] = np.mean(consensus_matrix)

    return consensus_scores


print("Running Consensus Clustering to choose K ...")
consensus_scores = consensus_clustering(
    X_scaled,
    k_min=2,
    k_max=7,
    n_iter=50,
    subsample_rate=0.8
)
print("Consensus scores (K -> score):", consensus_scores)

# 选共识得分最高的 K
best_k = max(consensus_scores, key=consensus_scores.get)
print("Best number of clusters selected by consensus clustering:", best_k)

# 画一下共识得分随 K 的变化，并保存到 clustering 目录
plt.figure()
plt.plot(list(consensus_scores.keys()),
         list(consensus_scores.values()),
         marker='o')
plt.xlabel("Number of clusters (K)")
plt.ylabel("Mean consensus score")
plt.title(f"Consensus clustering stability ({DATA_TYPE})")
plt.show()


# =========================
# 5. 用最佳 K 做最终聚类   
# =========================
kmeans = KMeans(
    n_clusters=best_k,
    random_state=RANDOM_STATE,
    n_init=50
)
cluster_labels = kmeans.fit_predict(X_scaled)


# =========================
# 6. PCA 可视化
# =========================
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

plt.figure()
plt.scatter(X_pca[:, 0], X_pca[:, 1], c=cluster_labels)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title(f"PCA Visualization ({DATA_TYPE}, k={best_k})")
plt.show()


# =========================
# 7. 保存聚类结果（到 clustering 目录）
# =========================##
result = pd.DataFrame({
    "Sample_ID": X.index,
    "Cluster": cluster_labels
})

result_path = os.path.join(
    clustering_dir,
    f"cluster_results_{DATA_TYPE}_seed{RANDOM_STATE}.csv"
)
result.to_csv(result_path, index=False)
print(f"Cluster results saved to: {result_path}")
