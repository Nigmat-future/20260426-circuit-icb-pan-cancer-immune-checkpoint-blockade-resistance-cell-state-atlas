# Supplementary Materials Manifest

## Supplementary Tables

- Table S1: `results/tables/melanoma_only/melanoma_state_signatures.yml`
  - Complete gene signatures for all six CIRCUIT states, including 20 up-regulation genes and 10 down-regulation genes per state.
- Table S2: `results/tables/melanoma_only/validation/melanoma_univariable_ci.tsv`
  - Per-cohort univariable AUROC and AUPRC with 95% bootstrap confidence intervals for all CIRCUIT states and benchmark signatures.
- Table S3: `results/tables/melanoma_only/validation/melanoma_cv_folds.tsv`
  - Fold-level metrics from repeated stratified 5-fold cross-validation across 50 repeats (250 total folds).

## Supplementary Figures

- Figure S1 source data: `results/tables/melanoma_only/melanoma_k_scan.tsv`
  - NMF reconstruction error and stability scan across component counts k = 2 to 8.
- Figure S2 source data: `results/tables/melanoma_only/melanoma_lineage_info.tsv`
  - Lineage overlap fractions for all six CIRCUIT states against curated lineage marker sets.

## Notes

- The main manuscript focuses on melanoma discovery and validation.
- Cross-cancer exploratory outputs are retained in `results/tables/melanoma_only/validation/generalization_univariable_ci.tsv` but are not part of the main manuscript narrative.
