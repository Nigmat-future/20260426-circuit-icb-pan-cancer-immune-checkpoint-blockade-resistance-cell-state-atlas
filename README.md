# CIRCUIT-ICB

**C**ell-state **I**dentification and **R**esistance **C**lassification **U**sing **I**ntegrated **T**ranscriptomic Inference for **I**mmune **C**heckpoint **B**lockade

[![Python](https://img.shields.io/badge/Python-3.11-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Workflow: Snakemake](https://img.shields.io/badge/Workflow-Snakemake-brightgreen.svg)](https://snakemake.readthedocs.io/)

This repository contains the code and result tables for the manuscript

> Rahim N. *CIRCUIT-ICB: NMF-based discovery of malignant cell states predicting immune checkpoint blockade resistance in melanoma.* Manuscript under review at *Bioinformatics*.

The manuscript source is [`manuscript/manuscript.md`](manuscript/manuscript.md). Every number in the main text and tables can be traced to a file under [`results/tables/melanoma_only/`](results/tables/melanoma_only/) (see the mapping below).

## What the manuscript reports

- **Discovery.** Non-negative matrix factorization (NMF) of 43 malignant-cell pseudobulk profiles from three melanoma scRNA-seq studies (Jerby-Arnon 2018, Pozniak 2022, Tirosh 2016), after excluding housekeeping genes and filtering markers by a specificity ratio of at least 2. Stability-guided selection gives k = 6 states.
- **Validation.** The six states are scored with a within-sample rank method in two pre-treatment anti-PD-1 melanoma cohorts, GSE78220 and GSE91061 (76 patients, 25 responders).
- **Results.** STATE_01 and STATE_05 are higher in non-responders, each with pooled AUROC 0.729 (permutation P = 0.049 and 0.047 against 1,000 size-matched random signatures). The best benchmark signature (IFNG / antigen presentation) reaches 0.616. A six-state elastic-net model reaches a repeated cross-validated AUROC of 0.694 (95% interval 0.600 to 0.761), against 0.560 for the benchmark signatures.

## Manuscript items and their source files

All paths are relative to `results/tables/melanoma_only/`.

| Manuscript item | Source file(s) |
|---|---|
| Table 1, Figure S1 (NMF k scan and stability) | `melanoma_k_scan.tsv` |
| Table 2 (top markers per state) | `melanoma_top15_markers.tsv`, `melanoma_nmf_markers.tsv` |
| Table 3 (pooled univariable AUROC, permutation P) | `validation/melanoma_univariable_ci.tsv`, `validation/melanoma_random_signature_pvalues.tsv` |
| Table 4 (repeated 5-fold CV, 50 repeats) | `validation/melanoma_cv_summary.tsv`, `validation/melanoma_cv_benchmark_summary.tsv`, `validation/melanoma_cv_combined_summary.tsv` |
| Section 3.4 calibration | `validation/melanoma_calibration.tsv` |
| Sections 3.3 and 3.5 direction stability | `validation/melanoma_direction_stability.tsv` |
| Random-effects meta-analysis (Methods 2.9) | `validation/melanoma_random_effects_meta.tsv` |
| Table S1 (full state signatures) | `melanoma_state_signatures.yml` |
| Table S2 (per-cohort AUROC / AUPRC with CI) | `validation/melanoma_univariable_ci.tsv` |
| Table S3 (fold-level CV metrics) | `validation/melanoma_cv_folds.tsv` |
| Figure S2 (lineage overlap) | `melanoma_lineage_info.tsv` |

The same mapping is listed in [`manuscript/supplementary_manifest.md`](manuscript/supplementary_manifest.md).

## Reproducing the manuscript

Tested with Python 3.11.

```bash
conda env create -f environment.yml
conda activate circuit-icb
python -m pip install -e ".[dev,manuscript]"
```

Input data (not tracked in git):

1. Malignant-cell pseudobulk CSV files from the public Zenodo compilation of ICB-treated scRNA-seq studies, placed under `data/raw/zenodo/malignant_pseudobulk/`.
2. Bulk RNA-seq and clinical tables from GEO, prepared by `scripts/07_prepare_geo_bulk.py` into `data/processed/bulk_icb/<cohort>/`. The cohort registry used by the melanoma analysis (`data/manifests/publish_v1_bulk_icb_cohorts.tsv`) is built by `scripts/21_build_publish_v1_cohorts.py`.

Then either run the workflow target

```bash
snakemake --cores 1 manuscript
```

or the three scripts directly:

```bash
python scripts/30_melanoma_only_pipeline.py   # discovery, scoring and validation -> results/tables/melanoma_only/
python scripts/40_draw_figures.py             # figures -> manuscript/figures/
python scripts/41_generate_docx.py            # manuscript/CIRCUIT_ICB_manuscript.docx
```

Unit tests: `pytest`.

## Analyses in this repository that are not part of the manuscript

These outputs are kept for transparency. They are not reported in the manuscript and should not be read as part of its claims.

- **Leave-one-cohort-out model (melanoma).** `scripts/30_melanoma_only_pipeline.py` also fits the elastic-net model on one melanoma cohort and tests it on the other, written to `validation/melanoma_leave_one_cohort.tsv`. With only two cohorts of different size and responder rate, this check is unstable: AUROC 0.60 when testing on GSE78220 and 0.37 when testing on GSE91061.
- **Melanoma states outside melanoma.** The same script scores the melanoma states in three non-melanoma cohorts (GSE126044, GSE135222, GSE176307), written to `validation/generalization_univariable_ci.tsv`. Pooled AUROCs are 0.52 to 0.59, and STATE_05 changes direction, so the melanoma states are not expected to transfer across cancer types.
- **Earlier pan-cancer analysis (`publish_v1`).** Before the melanoma-focused study, the framework was run on pan-cancer discovery data with a separate 7-state NMF and validated in six bulk cohorts (212 samples). Pooled AUROCs were 0.55 to 0.56, leave-one-cohort mean AUROC was 0.48, and none of the four primary states was significant against size-matched random signatures. This analysis is documented in [`docs/publish_v1_validation_report.md`](docs/publish_v1_validation_report.md) and [`docs/manuscript_ready_results_summary.md`](docs/manuscript_ready_results_summary.md). Its result tables are not committed. The repository name reflects this original pan-cancer scope.
- **Target-evidence modules** (DepMap, PRISM, LINCS, HPA, TCGA Xena, CPTAC adapter; scripts 09 to 27) come from the pan-cancer scope and are not used in the manuscript.
- **Project planning notes** from the original pan-cancer design are in [`docs/planning/`](docs/planning/).

## Repository layout

```text
config/              state, benchmark and pathway signature files; project config
docs/                method notes; pan-cancer (publish_v1) reports; planning/ notes
examples/            toy inputs for the smoke test
manuscript/          manuscript source, supplementary manifest, cover letter summary
results/tables/      melanoma_only/ holds every table behind the manuscript
scripts/             numbered pipeline steps (30 = manuscript analysis, 40-41 = figures and docx)
src/circuit_icb/     core package: NMF states, scoring, validation
tests/               pytest unit tests
workflow/Snakefile   Snakemake workflow (target `manuscript` reproduces the paper)
```

`data/` is created locally when the pipeline runs and is not tracked.

## Scope and limitations

The validation in the manuscript is retrospective and uses 76 patients from two cohorts. The states are candidates for testing in larger melanoma cohorts. They are not a clinical decision tool.
