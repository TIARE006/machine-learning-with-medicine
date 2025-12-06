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



## Folder-by-folder description

1. small RNA-seq/
raw/

GSE254878_smallRNAs_raw_counts_expression.csv
Raw small RNA counts from GEO (GSE254878).

Rows: features (miRNA / snoRNA / tRNA, etc.)

Columns: samples (patients).

The first row is a description row and is removed in the code.

The first column is used as the feature ID.

clustering/

cluster_results_smallRNA_seed42.csv
Final clustering result for small RNA-seq.

Columns:

Sample_ID: sample name (matches the columns of the raw expression file)

Cluster: cluster label (0 / 1, etc.), used later for DEG analysis.

deg/ (Differential Expression Results for small RNA)

All these files are generated by de_analysis_pipeline.py with DATA_TYPE = "smallRNA".

DEG_full_smallRNA_seed42.csv
Full differential expression table for all small RNA features. Columns:

gene (feature ID)

log2FC (Cluster 0 vs Cluster 1)

p_value

FDR (BH-corrected)

DEG_sig_strict_smallRNA_FDR0.05_log2FC1.0_seed42.csv
“Strict” significant small RNAs:

FDR < 0.05

|log2FC| > 1.0
→ This is the main DEG list for reporting.

DEG_sig_relaxed_smallRNA_FDR0.1_log2FC0.5_seed42.csv
“Relaxed” DEGs, used mostly for enrichment and target mapping:

FDR < 0.10

|log2FC| > 0.5

DEG_top200_smallRNA_seed42.csv
Top 200 features ranked by raw p-value (exploratory use).

DEG_up_genes_smallRNA_seed42.txt
One ID per line: small RNAs up-regulated in cluster 0 vs cluster 1 (relaxed cutoff).
Used as input smallRNA list for miRNA→mRNA mapping.

DEG_down_genes_smallRNA_seed42.txt
Down-regulated small RNAs (relaxed cutoff).

pathway/ (GO enrichment for smallRNA and their targets)

Pathway_up_smallRNA_seed42.csv
GO Biological Process enrichment for up-regulated small RNAs
(directly using small RNA IDs; mainly for reference).

Pathway_down_smallRNA_seed42.csv
GO-BP enrichment for down-regulated small RNAs.

Pathway_targets_up_smallRNA_seed42.csv
GO-BP enrichment for target mRNAs of up-regulated smallRNAs
(i.e., after miRNA→mRNA mapping, enrichment is run on the mRNA gene list).

Pathway_targets_down_smallRNA_seed42.csv
GO-BP enrichment for target mRNAs of down-regulated smallRNAs (may be empty).

targets/ (miRNA → mRNA target lists)

Targets_up_from_smallRNA_seed42.txt
Unique mRNA target genes of up-regulated smallRNAs (one gene symbol per line).

Targets_down_from_smallRNA_seed42.txt
Unique mRNA target genes of down-regulated smallRNAs.

These files are used as the gene lists for Pathway_targets_*.csv.

2. RNA-seq/
raw/

GSE254877_raw_counts_expression.csv
Raw mRNA counts from GEO (GSE254877).

Rows: genes

Columns: samples

First row is a description row (removed in the code).

First column is used as gene ID.

clustering/

cluster_results_RNA_seed42.csv
Clustering result for RNA-seq, same format as for small RNA:

Sample_ID

Cluster

deg/ (Differential Expression Results for mRNA)

Generated by de_analysis_pipeline.py with DATA_TYPE = "RNA".

DEG_full_RNA_seed42.csv
DEG table for all genes.

DEG_sig_strict_RNA_FDR0.05_log2FC1.0_seed42.csv
Strict-significance mRNA DEGs (FDR < 0.05 and |log2FC| > 1.0).

DEG_sig_relaxed_RNA_FDR0.1_log2FC0.5_seed42.csv
Relaxed mRNA DEGs (FDR < 0.10 and |log2FC| > 0.5).

DEG_top200_RNA_seed42.csv
Top 200 genes by p-value.

DEG_up_genes_RNA_seed42.txt / DEG_down_genes_RNA_seed42.txt
Up-/down-regulated genes (relaxed cutoff). Used as input for enrichment.

pathway/ (GO enrichment for mRNA)

Pathway_up_RNA_seed42.csv
GO-BP enrichment result for up-regulated mRNA genes.

Pathway_down_RNA_seed42.csv
GO-BP enrichment result for down-regulated mRNA genes.

3. reference/

Files used to build the miRNA→target gene mapping.

gene_attribute_edges.txt.gz
Raw edge list downloaded from the miRTarBase / Enrichr resource.
Contains columns like:
source, source_desc, source_id, target, target_desc, target_id, weight.

mirna_target_human.csv
Preprocessed human miRNA–target table.
Built from gene_attribute_edges.txt.gz by build_mirna_target_table.py.

Columns:

miRNA – miRNA ID (e.g. hsa-miR-154-5p)

TargetGene – human gene symbol (e.g. DICER1)

This file is used by de_analysis_pipeline.py when DATA_TYPE = "smallRNA" to map
differential miRNAs to their target mRNAs before running GO-BP enrichment.

Quick usage note

Set DATA_TYPE = "smallRNA" or "RNA" in de_analysis_pipeline.py.

The script will:

Read the corresponding raw/ expression file and clustering/ labels.

Compute DEGs (full / strict / relaxed / Top N).

For RNA-seq: run GO-BP enrichment directly on up/down DEG gene lists.

For small RNA-seq:

map differential miRNAs to target mRNAs using reference/mirna_target_human.csv;

run GO-BP enrichment on the target mRNA gene lists.

---


