# Fig. 1 Workflow

The active target is:

`results/fig1_final/intNMF_K4_PA_PQ_final_main_figure/Fig1_intNMF_K4_PA_PQ_final_main_figure.pdf`

Run all scripts from the project root. The active Fig. 1 workflow is R-only and uses the project-local `Rlibs/` library.

## Active Scripts

### Package setup

- `00_setup/00_setup_R_packages.R`: installs the local R dependencies under `Rlibs/`.

### Source-input generation

- `01_expression_inputs/10_build_fig1_expression_inputs.R`: generates 14-PCD ssGSEA scores and the pyroptosis marker matrix directly from raw RNA counts.
- `01_expression_inputs/11_build_intnmf_source_inputs.R`: generates sample metadata and the technical-adjusted two-block intNMF inputs directly from raw RNA, small-RNA, and GEO phenotype files.
- `01_expression_inputs/pcd_gene_sets.csv`: versioned 14-PCD gene definitions and source references.
- `01_expression_inputs/locked_intnmf_feature_spec.csv`: locked 4500-feature intNMF model specification.
- `12_audit_fig1_lineage.R`: validates panel A-H data lineage and source hashes.

See `FIG1_DATA_LINEAGE.md` for the panel-level provenance table.

### intNMF clustering

- `02_intnmf/20_train_reproducible_two_block_intnmf.R`: runs the active fixed-seed CPI selection and 20-seed final K4 fit from raw-derived inputs.
- `04_annotate_reproducible_K4.R`: applies the locked PA/PQ annotation rule to the reproducible K4 solution.
- `02_intnmf/43_validate_raw_to_fit_reproduction.R`: records the independent seed-11 raw-to-fit reproduction check.

The final plotting workflow reads:

- `results/intNMF_reproducible_two_block/CPI_summary.csv`
- `results/intNMF_reproducible_two_block/two_block_final_K4_labels.csv`
- `results/intNMF_reproducible_two_block/final_consensus_matrix.csv`

### Pyroptosis inputs

- `03_pyroptosis/32_all84_pyro_module_main_figure.R`: generates all-84 pyroptosis module scores, validation scores, and metadata from RNA counts and final K4 labels.

### Final plotting

- `05_plotting/46_make_final_Fig1_intNMF_K4_PA_PQ.R`: generates panels A-H and assembles the current final Fig. 1 PDF/PNG.

## Reproduce Current Fig. 1

```powershell
Rscript scripts/fig1/00_setup/00_setup_R_packages.R
Rscript scripts/fig1/01_expression_inputs/11_build_intnmf_source_inputs.R
Rscript scripts/fig1/02_intnmf/20_train_reproducible_two_block_intnmf.R
Rscript scripts/fig1/04_annotate_reproducible_K4.R
Rscript scripts/fig1/01_expression_inputs/10_build_fig1_expression_inputs.R
Rscript scripts/fig1/12_audit_fig1_lineage.R
Rscript scripts/fig1/03_pyroptosis/32_all84_pyro_module_main_figure.R
Rscript scripts/fig1/05_plotting/46_make_final_Fig1_intNMF_K4_PA_PQ.R
```

## Audit

```powershell
Rscript scripts/fig1/12_audit_fig1_lineage.R
```

Expected output:

```text
Fig. 1 lineage audit passed for panels A-H.
```

## Retired Material

Older Fig. 2-Fig. 7, exploratory, legacy Fig. 1 layouts, CPS scoring, historical K2/K3 attempts, and external-cohort figure drafts are not part of the current active dependency chain.
