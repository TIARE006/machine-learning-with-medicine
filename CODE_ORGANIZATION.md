# Code Organization

The project is currently organized around the active Fig. 1 reproducibility workflow.

## Active Code

- `scripts/fig1/00_setup/`: project-local R package setup.
- `scripts/fig1/01_expression_inputs/`: raw-derived input generation, locked feature specification, and fixed gene-set definitions.
- `scripts/fig1/02_intnmf/`: reproducible two-block intNMF training, K selection, stability checks, and validation.
- `scripts/fig1/03_pyroptosis/`: pyroptosis module and validation-score generation.
- `scripts/fig1/04_annotate_reproducible_K4.R`: locked PA/PQ annotation.
- `scripts/fig1/05_plotting/46_make_final_Fig1_intNMF_K4_PA_PQ.R`: final Fig. 1 panel and assembly script.
- `scripts/fig1/12_audit_fig1_lineage.R`: machine-checkable Fig. 1 lineage and hash audit.

## Active Outputs

- `results/fig1_inputs/`: generated Fig. 1 inputs plus source-hash and panel-lineage audits.
- `results/intNMF_reproducible_two_block/`: active model fit, labels, consensus matrix, stability summaries, package versions, and reproduction validation.
- `results/fig1_final/all84_pyro_module_Euclidean_K2_main/`: Fig. 1 pyroptosis module inputs used by final plotting.
- `results/fig1_final/intNMF_K4_PA_PQ_final_main_figure/`: current final Fig. 1 deliverables.

## Rules

1. Run scripts from the project root unless a script explicitly says otherwise.
2. Keep generated outputs under `results/`, not under code directories.
3. Keep root inputs, locked definitions, generated inputs, and final plotted data connected through manifest hashes.
4. Rerun `scripts/fig1/12_audit_fig1_lineage.R` after changing any Fig. 1 input, producer script, or generated source.
5. Fig. 2-Fig. 7 and exploratory figure materials are outside the current active scope.
