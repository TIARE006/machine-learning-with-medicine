# RNA & smallRNA Analysis Pipeline: Clustering, Differential Expression, and Pathway Enrichment


# Contents
- [Project Overview](#project-overview)
- [Installations](#installations)
- [Full Pipeline Overview](#Full-Pipeline-Overview)
- [Data Folder Overview](#data-folder-overview)
- [Quick Start](#quick-start)
- [Datasets](#datasets)
- [Tasks \& Baselines](#tasks--baselines)
- [References](#references)


# Project Overview

This project implements a complete analysis pipeline for two GEO datasets:

- **small RNA-seq (GSE254878)**
- **RNA-seq (GSE254877)**

The pipeline performs four major steps:

1. **Unsupervised clustering** based on normalized expression data  
2. **Differential expression analysis (DEG)** using strict and relaxed thresholds  
3. **Target gene mapping** for small RNA (miRNA → mRNA) through curated reference tables  
4. **GO Biological Process pathway enrichment**, for both directly differential genes and miRNA-derived target genes

All intermediate and final outputs — clustering labels, DEG tables, gene lists, and pathway enrichment results — are saved under a structured `data/` directory.


# Installation

## 1. Clone the repository
```bash
git clone https://github.com/TIARE006/machine-learning-with-medicine
cd machine-learning-with-medicine
```

## Step 2. Set up the environment:
```bash
# Set up the environment
conda create -n mlomics python=3.14.0
conda activate 
```

## Step 3. Install requirements:
```bash
pip install -r requirements.txt
```

## Step 4. Download datasets:
```bash
# Option 1: automatically download GEO datasets
python scripts/download_data.py

# Option 2: manually download from GEO
# GSE254877: RNA-seq
# GSE254878: small RNA-seq
```


# Full Pipeline Overview

The script `run_full_pipeline.py` performs the complete multiomics workflow:

---

## 1. Clustering (SNF)
- Reads multi-omics matrices
- Performs Similarity Network Fusion
- Outputs cluster assignments

## 2. Differential Expression Analysis (DEG)
- Computes DEGs using multiple thresholds:
  - Full list
  - Strict threshold
  - Relaxed threshold
  - Top N genes
- Outputs DEG tables for each comparison

## 3. miRNA Target Mapping (smallRNA-seq only)
- Maps differential miRNAs to target mRNAs
- Generates up-/down-regulated target lists

## 4. GO Biological Process Enrichment
- GO-BP enrichment for:
  - RNA-seq DEGs
  - miRNA-derived target genes
- Saves enriched pathway tables

## 5. Final structured output
All results are saved under `data/<DATA_TYPE>/integrated_results/`







# Data Folder Overview

This folder contains all data used by the clustering and DEG + pathway analysis pipelines for:

- **small RNA-seq**: GSE254878  
- **RNA-seq**: GSE254877  

The Python scripts (e.g. `cluster_analysis1.py`, `de_analysis_pipeline.py`) assume the following structure **relative to the project root**:

```text
data/
├── small RNA-seq/
│   ├── raw/
│   │   └── GSE254878_smallRNAs_raw_counts_expression.csv
│   ├── clustering/
│   │   └── cluster_results_smallRNA_seed42.csv
│   ├── deg/
│   │   ├── DEG_full_smallRNA_seed42.csv
│   │   ├── DEG_sig_strict_smallRNA_FDR0.05_log2FC1.0_seed42.csv
│   │   ├── DEG_sig_relaxed_smallRNA_FDR0.1_log2FC0.5_seed42.csv
│   │   ├── DEG_top200_smallRNA_seed42.csv
│   │   ├── DEG_up_genes_smallRNA_seed42.txt
│   │   └── DEG_down_genes_smallRNA_seed42.txt
│   ├── pathway/
│   │   ├── Pathway_up_smallRNA_seed42.csv
│   │   ├── Pathway_down_smallRNA_seed42.csv
│   │   ├── Pathway_targets_up_smallRNA_seed42.csv
│   │   └── Pathway_targets_down_smallRNA_seed42.csv
│   └── targets/
│       ├── Targets_up_from_smallRNA_seed42.txt
│       └── Targets_down_from_smallRNA_seed42.txt
│
├── RNA-seq/
│   ├── raw/
│   │   └── GSE254877_raw_counts_expression.csv
│   ├── clustering/
│   │   └── cluster_results_RNA_seed42.csv
│   ├── deg/
│   │   ├── DEG_full_RNA_seed42.csv
│   │   ├── DEG_sig_strict_RNA_FDR0.05_log2FC1.0_seed42.csv
│   │   ├── DEG_sig_relaxed_RNA_FDR0.1_log2FC0.5_seed42.csv
│   │   ├── DEG_top200_RNA_seed42.csv
│   │   ├── DEG_up_genes_RNA_seed42.txt
│   │   └── DEG_down_genes_RNA_seed42.txt
│   └── pathway/
│       ├── Pathway_up_RNA_seed42.csv
│       └── Pathway_down_RNA_seed42.csv
│
└── reference/
    ├── gene_attribute_edges.txt.gz
    └── mirna_target_human.csv

```

## Project Structure

```
MACHINE-LEARNING-WITH-...

├── config/
│   └── settings.yaml                 # 项目参数配置

├── data/
│   ├── integrated_results/           # 整合分析结果
│   ├── lncRNA-seq/                   # lncRNA 数据
│   ├── reference/                    # 参考数据库（miRNA/基因映射等）
│   ├── RNA-seq/                      # RNA-seq 数据
│   └── small RNA-seq/                # small RNA 数据
│
│   └── DOE_Cachexia Molecular Typing.md  # 说明文档（命名中含空格）

├── R/
│   └── run_deseq2_by_snf.R           # R 语言差异分析脚本

├── scripts/
│   ├── cluster_analysis_singleomics.py   # 单组学聚类分析
│   └── run_full_pipeline.py              # 全流程分析脚本

├── src/
│   ├── __pycache__/                      # Python 缓存文件
│   ├── __init__.py                       # 包初始化
│   ├── methods_utils.py                  # 工具函数
│   └── multiomics_snf_v2.py              # SNF 多组学算法实现

└── README.md                              # 项目文档
```



## Folder-by-Folder Description

## 1. `small RNA-seq/`
用于处理 **GSE254878 small RNA-seq** 数据的所有文件。

###  `raw/` — 原始表达矩阵

**GSE254878_smallRNAs_raw_counts_expression.csv**  
- Small RNA 原始 count 矩阵  
- 行：miRNA / snoRNA / tRNA 等  
- 列：样本  
- 第一行为描述行（自动移除）  
- 第一列为 feature ID  

###  `clustering/` — 聚类结果

**cluster_results_smallRNA_seed42.csv**  
包含以下列：  
- `Sample_ID`  
- `Cluster`（0/1/...）  

用于后续 DEG 分析。

###  `deg/` — 差异表达（DEG）结果  
由 `de_analysis_pipeline.py`（DATA_TYPE="smallRNA"）自动生成。

#### Main DEG tables
- **DEG_full_smallRNA_seed42.csv**  
  全部 small RNAs 的差异表达结果（gene, log2FC, p_value, FDR）。

- **DEG_sig_strict_smallRNA_FDR0.05_log2FC1.0_seed42.csv**  
  严格阈值：  
  - FDR < 0.05  
  - |log2FC| > 1.0  

- **DEG_sig_relaxed_smallRNA_FDR0.1_log2FC0.5_seed42.csv**  
  宽松阈值（用于富集和 target mapping）：  
  - FDR < 0.10  
  - |log2FC| > 0.5  

- **DEG_top200_smallRNA_seed42.csv**  
  p-value 最小的前 200 个 feature。

#### Up/down lists
- **DEG_up_genes_smallRNA_seed42.txt** — 上调 smallRNAs  
- **DEG_down_genes_smallRNA_seed42.txt** — 下调 smallRNAs  

###  `pathway/` — SmallRNA 与其 Target 的 GO-BP 富集

- **Pathway_up_smallRNA_seed42.csv** — 上调 smallRNAs 的 GO-BP 富集  
- **Pathway_down_smallRNA_seed42.csv** — 下调 smallRNAs 的 GO-BP 富集  
- **Pathway_targets_up_smallRNA_seed42.csv** — smallRNA → target mRNA → GO-BP 富集  
- **Pathway_targets_down_smallRNA_seed42.csv** — 下调 smallRNAs target mRNA 富集  

###  `targets/` — miRNA → mRNA Target 映射

- **Targets_up_from_smallRNA_seed42.txt** — 上调 smallRNAs 的 mRNA targets  
- **Targets_down_from_smallRNA_seed42.txt** — 下调 smallRNAs 的 mRNA targets  

## 2. `RNA-seq/`
用于处理 **GSE254877 RNA-seq** 数据。

###  `raw/` — 原始表达矩阵

**GSE254877_raw_counts_expression.csv**  
- 原始 mRNA count 矩阵  
- 行：基因  
- 列：样本  
- 第一行描述行移除  
- 第一列为 gene ID  

###  `clustering/` — RNA-seq 聚类结果

**cluster_results_RNA_seed42.csv**  
包含：  
- `Sample_ID`  
- `Cluster`  

###  `deg/` — mRNA 差异表达结果  
由 `de_analysis_pipeline.py` 生成。

- **DEG_full_RNA_seed42.csv**  
- **DEG_sig_strict_RNA_FDR0.05_log2FC1.0_seed42.csv**  
- **DEG_sig_relaxed_RNA_FDR0.1_log2FC0.5_seed42.csv**  
- **DEG_top200_RNA_seed42.csv**  
- **DEG_up_genes_RNA_seed42.txt** / **DEG_down_genes_RNA_seed42.txt**

###  `pathway/` — mRNA 的 GO-BP 富集

- **Pathway_up_RNA_seed42.csv**  
- **Pathway_down_RNA_seed42.csv**

## 3. `reference/` — SmallRNA Target 映射所需数据库

### **gene_attribute_edges.txt.gz**
来自 miRTarBase / Enrichr 的原始边表：  
包含 `source`, `target`, `weight` 等字段。

用于构建 miRNA → mRNA 映射。

### **mirna_target_human.csv**
由脚本构建后得到的人类 miRNA → TargetGene 映射表：

| miRNA | TargetGene |
|-------|-------------|
| hsa-miR-154-5p | DICER1 |
| ... | ... |

供 smallRNA target mapping 与 pathway enrichment 使用。

## Quick Usage Summary

在 `de_analysis_pipeline.py` 设置：

```python
DATA_TYPE = "smallRNA"
# 或
DATA_TYPE = "RNA"



