# Fig. 1 intNMF K4 PA/PQ Reproducible Analysis

This project is currently scoped to the Fig. 1 analysis only.

Active final output:

`results/fig1_final/intNMF_K4_PA_PQ_final_main_figure/Fig1_intNMF_K4_PA_PQ_final_main_figure.pdf`

## Reproducibility Status

The active Fig. 1 workflow is designed to be traceable from retained raw inputs to the final plotted values.

- Root inputs, generated inputs, model outputs, and final panel inputs are recorded with SHA-256 hashes.
- The active intNMF solution uses a locked two-block feature specification.
- K=4 is selected by fixed-seed CPI screening.
- Twenty fixed K4 seed runs are stable with pairwise ARI mean and minimum equal to 1.
- An independent seed-11 raw-to-fit rerun reproduces the active result with ARI 1 and consensus-matrix maximum absolute difference 0.

Machine-readable audit outputs:

- `results/fig1_inputs/Fig1_panel_data_lineage.csv`
- `results/fig1_inputs/Fig1_source_hash_audit.csv`
- `results/intNMF_reproducible_two_block/reproducibility_summary.csv`
- `results/intNMF_reproducible_two_block/reproduction_validation/raw_to_fit_reproduction_metrics.csv`

## Setup

Run all commands from the project root.

```powershell
Rscript scripts/fig1/00_setup/00_setup_R_packages.R
```

The setup script installs required R packages into the project-local `Rlibs/` directory. The project `.Rprofile` also prepends `Rlibs/` to `.libPaths()` when R starts in this directory.

## Reproduce Fig. 1

```powershell
Rscript scripts/fig1/01_expression_inputs/11_build_intnmf_source_inputs.R
Rscript scripts/fig1/02_intnmf/20_train_reproducible_two_block_intnmf.R
Rscript scripts/fig1/04_annotate_reproducible_K4.R
Rscript scripts/fig1/01_expression_inputs/10_build_fig1_expression_inputs.R
Rscript scripts/fig1/12_audit_fig1_lineage.R
Rscript scripts/fig1/03_pyroptosis/32_all84_pyro_module_main_figure.R
Rscript scripts/fig1/05_plotting/46_make_final_Fig1_intNMF_K4_PA_PQ.R
```

The lineage audit can be rerun at any time:

```powershell
Rscript scripts/fig1/12_audit_fig1_lineage.R
```

Expected message:

```text
Fig. 1 lineage audit passed for panels A-H.
```

## Active Directories

- `scripts/fig1/`: active Fig. 1 workflow scripts.
- `scripts/fig1/01_expression_inputs/`: raw-derived Fig. 1 input builders and locked gene/feature definitions.
- `scripts/fig1/02_intnmf/`: reproducible intNMF training and validation scripts.
- `scripts/fig1/03_pyroptosis/`: pyroptosis module scoring for Fig. 1.
- `scripts/fig1/05_plotting/`: final Fig. 1 plotting script.
- `results/fig1_inputs/`: generated Fig. 1 inputs and lineage audits.
- `results/intNMF_reproducible_two_block/`: active reproducible K4 model outputs.
- `results/fig1_final/intNMF_K4_PA_PQ_final_main_figure/`: final Fig. 1 panels and assembled figure.

## Scope Note

Older Fig. 2-Fig. 7, exploratory, and legacy figure materials are not part of the active reproducible analysis scope.
