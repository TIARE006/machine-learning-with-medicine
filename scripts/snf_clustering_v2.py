#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
snf_clustering_v2.py

根目录入口脚本：
- 调用 src.multiomics_snf_v2 中封装好的 SNF + EarlyFusion 全流程
- 自动输出：
    - SNF 与 EarlyFusion 聚类结果（CSV）
    - UMAP / PCA 可视化图
    - silhouette 对比（TXT）

真正的代码逻辑全部放在 src/multiomics_snf_v2.py 中。
"""

from src.multiomics_snf_v2 import run_snf_v2_pipeline


def main():
    print("===== Running SNF Clustering V2.0 Pipeline =====")
    run_snf_v2_pipeline()
    print("===== SNF Clustering V2.0 Finished =====")


if __name__ == "__main__":
    main()

