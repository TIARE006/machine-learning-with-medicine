import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


# =========================
# 🔽 聚类数据类型开关 🔽
# 可选： "smallRNA" 或 "RNA"
# =========================
DATA_TYPE = "smallRNA"   # smallRNA 或 RNA


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

print(f"✅ Current clustering mode: {DATA_TYPE}")
print(f"✅ Using file: {file_path}")


# =========================
# 2. 读取并清洗数据（完全修复版）
# =========================
df = pd.read_csv(file_path, low_memory=False)

# 删除第一行异常描述行
df = df.drop(index=0)

# 自动使用第一列作为基因ID
gene_col = df.columns[0]
print("✅ Using gene ID column:", gene_col)
df = df.set_index(gene_col)

# smallRNA 才有 type 列
if "type" in df.columns:
    df = df.drop(columns=["type"])

# 强制转为数值
df = df.apply(pd.to_numeric, errors="coerce")

# ✅ 只删除“整行全空”的情况
df = df.dropna(how="all")

# ✅ 用 0 填补剩余 NaN（RNA表达中0是合理的）
df = df.fillna(0)

# 转置：样本为行
X = df.T

print("✅ Expression matrix shape:", X.shape)


# =========================
# 3. 标准化
# =========================
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)


# =========================
# 4. 寻找最佳聚类数
# =========================
range_n_clusters = range(2, 8)
silhouette_scores = []

for k in range_n_clusters:
    kmeans = KMeans(n_clusters=k, random_state=42)
    labels = kmeans.fit_predict(X_scaled)
    silhouette_scores.append(silhouette_score(X_scaled, labels))

plt.figure()
plt.plot(range_n_clusters, silhouette_scores, marker='o')
plt.xlabel("Number of clusters (k)")
plt.ylabel("Silhouette Score")
plt.title(f"Silhouette Analysis ({DATA_TYPE})")
plt.show()

best_k = range_n_clusters[np.argmax(silhouette_scores)]
print("✅ Best number of clusters:", best_k)


# =========================
# 5. 执行最终聚类
# =========================
kmeans = KMeans(n_clusters=best_k, random_state=42)
cluster_labels = kmeans.fit_predict(X_scaled)


# =========================
# 6. PCA 可视化
# =========================
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

plt.figure()
plt.scatter(X_pca[:, 0], X_pca[:, 1])
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title(f"PCA Visualization ({DATA_TYPE}, k={best_k})")
plt.show()


# =========================
# 7. 保存结果（分类存放到 data 目录）
# =========================

# 根据数据类型选择输出文件夹
if DATA_TYPE == "smallRNA":
    output_dir = os.path.join(BASE_DIR, "data", "small RNA-seq")
else:
    output_dir = os.path.join(BASE_DIR, "data", "RNA-seq")

# 如果文件夹不存在则创建
os.makedirs(output_dir, exist_ok=True)

result = pd.DataFrame({
    "Sample_ID": X.index,
    "Cluster": cluster_labels
})

result_path = os.path.join(
    output_dir,
    f"cluster_results_{DATA_TYPE}.csv"
)

result.to_csv(result_path, index=False)

print(f"✅ Cluster results saved to: {result_path}")



