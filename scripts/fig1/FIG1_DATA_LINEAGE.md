# Fig. 1 Data Lineage

The final target is:

`results/fig1_final/intNMF_K4_PA_PQ_final_main_figure/Fig1_intNMF_K4_PA_PQ_final_main_figure.pdf`

Every plotted value must be connected to a retained producer script and a retained root input.

| Panels | Direct plotted data | Producer | Root source |
|---|---|---|---|
| A, D, E | `results/fig1_inputs/generated/GSE254877_14PCD_ssGSEA_zscores.csv` | `01_expression_inputs/10_build_fig1_expression_inputs.R` | GSE254877 raw RNA counts and versioned `pcd_gene_sets.csv` |
| B top | `results/intNMF_reproducible_two_block/CPI_summary.csv` | `02_intnmf/20_train_reproducible_two_block_intnmf.R` | Raw RNA, lncRNA, small-RNA, and locked feature specification |
| B bottom, C | `results/intNMF_reproducible_two_block/final_consensus_matrix.csv` | `02_intnmf/20_train_reproducible_two_block_intnmf.R` | Same raw-to-fit training chain |
| D-H grouping | `results/intNMF_reproducible_two_block/two_block_final_K4_labels.csv` | `04_annotate_reproducible_K4.R` | Reproducible K4 fit and locked pyroptosis annotation rule |
| F, G | `results/fig1_final/all84_pyro_module_Euclidean_K2_main/all84_pyro_module_scores_by_sample.csv` | `03_pyroptosis/32_all84_pyro_module_main_figure.R` | GSE254877 RNA counts, final K4 labels, and locked module genes |
| F, G | `results/fig1_final/all84_pyro_module_Euclidean_K2_main/all84_K2_validation_scores_by_sample.csv` | `03_pyroptosis/32_all84_pyro_module_main_figure.R` | GSE254877 RNA counts, final K4 labels, and held-out genes |
| G, H annotations | `results/fig1_final/all84_pyro_module_Euclidean_K2_main/all84_sample_metadata_used.csv` | `03_pyroptosis/32_all84_pyro_module_main_figure.R` | GEO phenotype, expression-inferred sex, and final K4 labels |
| H expression | `results/fig1_inputs/generated/GSE254877_pyroptosis_marker_gene_zmatrix.csv` | `01_expression_inputs/10_build_fig1_expression_inputs.R` | Raw RNA counts and fixed marker list |

`12_audit_fig1_lineage.R` verifies that all direct inputs, producer scripts, root files, and generated-source hashes match their manifests. It writes machine-readable audits under `results/fig1_inputs/`.

The 14-PCD gene sets are project-versioned analysis definitions. Their source references are stored per gene in `pcd_gene_sets.csv`; changing that file changes the analysis and must regenerate panels A, D, E, and G.

## Reproducibility Status

The active model is trained end to end from retained raw count files. The 4500-feature model specification is versioned in `locked_intnmf_feature_spec.csv`. K=4 is selected by fixed-seed CPI analysis. Twenty fixed K4 seed runs have pairwise ARI mean and minimum equal to 1. The medoid-seed fit, PA/PQ annotation, downstream scores, and final figure are reproducible from retained code and inputs.

An independent seed-11 rerun is recorded under `results/intNMF_reproducible_two_block/reproduction_validation/`: adjusted Rand index is 1 and the consensus-matrix maximum absolute difference is 0 relative to the active final fit.
