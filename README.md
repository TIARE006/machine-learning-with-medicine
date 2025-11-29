Cancer Cachexia RNA-seq & smallRNA Clustering and Differential Expression Pipeline

This project implements a complete analysis workflow for RNA-seq and small RNA-seq (miRNA) data related to cancer cachexia. The pipeline includes preprocessing, clustering, differential gene expression (DEG), and pathway enrichment. All outputs are automatically organized into clean folder structures for easy interpretation.

Directory Structure

RNA-seq/
│
├── raw/                ← 原始数据（不要动）
│     ├── *_expression.csv
│     ├── *_features.csv
│     ├── *_pheno.csv
│     └── *_raw_counts_expression.csv
│
├── clustering/         ← 聚类结果（输入给后续 DEG）
│     ├── cluster_results_RNA.csv
│     ├── cluster_results_RNA_seed42.csv
│     ├── cluster_dbscan_RNA.csv
│     ├── cluster_hdbscan_RNA.csv
│
├── deg/                ← 差异分析结果（重点！！！）
│     ├── DEG_full_RNA_seed42.csv
│     ├── DEG_sig_RNA_*.csv
│     ├── DEG_up_genes_RNA_seed42.txt
│     ├── DEG_down_genes_RNA_seed42.txt
│
├── pathway/            ← 富集分析（论文图/表来自这里）
│     ├── Pathway_up_RNA_seed42.csv
│     └── Pathway_down_RNA_seed42.csv
│
└── other/              ← 不重要，系统自动归类


Explanation of each folder:

raw/
Contains original GEO files: raw counts, processed expression matrices, phenotype metadata, and annotation files. These are inputs and should never be modified manually.

clustering/
Contains clustering assignments and PCA visualization results. Example:
cluster_results_RNA_seed42.csv
cluster_dbscan_RNA.csv
cluster_hdbscan_RNA.csv

deg/
Contains differential expression results including:

Full DEG table

Significant DEGs

Up-regulated gene list

Down-regulated gene list

Example files:
DEG_full_RNA_seed42.csv
DEG_sig_RNA_FDR0.05_log2FC1.0_seed42.csv
DEG_up_genes_RNA_seed42.txt
DEG_down_genes_RNA_seed42.txt

pathway/
Contains pathway enrichment results from GO Biological Process (GO-BP) or KEGG.
Files include enriched pathways for both up-regulated and down-regulated gene sets.
Example:
Pathway_up_RNA_seed42.csv
Pathway_down_RNA_seed42.csv

other/
Automatically sorted files that do not match classification rules. Usually not used in analysis.

Analysis Pipeline Overview

Step 1 — Preprocessing

Load raw counts expression matrix

Remove non-numeric lines

Set first column as gene ID

Transpose so rows = samples, columns = genes

Apply Z-score normalization

Step 2 — Clustering

Consensus clustering is used to determine the stable number of clusters

Final K-means implemented with n_init=50 for higher stability

PCA visualization is provided

Output saved to clustering/ directory

Step 3 — Differential Expression (DEG)
Script: de_analysis_pipeline.py

Loads clustering labels

Performs Welch’s t-test

Computes log2FC

Adjusts p-values using Benjamini–Hochberg FDR

Saves:

Full DEG table

Significant DEGs

Up-regulated list

Down-regulated list
Outputs saved under deg/ directory.

Step 4 — Pathway Enrichment
Performed using gseapy (Enrichr) with GO Biological Process (GO-BP).
Uses up/down gene lists generated from DEG.
Outputs saved under pathway/ directory.

How to Run

Run clustering:
python cluster_analysis1.py

Output:
data/<dataset>/clustering/cluster_results_<DATA_TYPE>_seed42.csv

Run DEG + Pathway analysis:
python de_analysis_pipeline.py

Output:
data/<dataset>/deg/
data/<dataset>/pathway/

Key Output Files to Use for Interpretation

Clustering (for downstream grouping):
cluster_results_<DATA_TYPE>_seed42.csv

Differential Expression:
DEG_sig_<DATA_TYPE>FDR0.05_log2FC1.0_seed42.csv
DEG_up_genes<DATA_TYPE>seed42.txt
DEG_down_genes<DATA_TYPE>_seed42.txt

Pathway Interpretation:
Pathway_up_<DATA_TYPE>seed42.csv
Pathway_down<DATA_TYPE>_seed42.csv

These are the files needed for biological interpretation, scientific figures, and writing the Results section of a paper.

Recommended Next Steps

Possible visualizations and extensions of the pipeline include:

Volcano plot for DEGs

Heatmap of top DEGs

Dotplot or barplot for GO/KEGG enrichment

PCA visualization grouped by clusters

miRNA–mRNA integration analysis

I can provide ready-to-run scripts for all these if needed.

Citation Notes

If used in publications, please acknowledge:

GEO datasets GSE254877 & GSE254878

gseapy

scikit-learn

End of README.