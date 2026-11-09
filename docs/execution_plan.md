# Execution Plan

> **Scope note.** This document belongs to the earlier pan-cancer `publish_v1` analysis and is not part of the melanoma manuscript (`manuscript/manuscript.md`). See the README for the files behind the manuscript.

## 72-Hour Milestone

1. Freeze public source registry and Zenodo metadata manifest.
2. Generate normalized discovery, file, and treatment manifests.
3. Run example state-score and target-priority smoke tests.
4. Curate first 3 to 5 external bulk ICB cohorts and response-label rules.
5. Define the primary malignant-state discovery method and one sensitivity method.

Exit criterion:

- `pytest` passes
- example workflow outputs exist
- `data/manifests/icb_sc_discovery_datasets.tsv` is generated from the public workbook

## 2-Week Milestone

1. Download selected discovery RDS or pseudobulk files.
2. Build patient/sample-level malignant state activity table.
3. Run preliminary consensus NMF or curated program scoring.
4. Curate bulk ICB cohorts and compute state scores.
5. Produce Figure 1 and draft Figure 2/3 tables.

Exit criterion:

- at least two states reproduce across studies
- at least three external bulk cohorts score successfully
- response models are evaluated with locked endpoints

## 6-8 Week Milestone

1. Complete leave-one-study-out and leave-one-cancer-type-out validation.
2. Benchmark against immune and tumor-intrinsic baselines.
3. Add TCGA pan-cancer characterization.
4. Add CPTAC, DepMap, PRISM/GDSC, LINCS, GTEx/HPA evidence.
5. Generate target and compound priority tables.
6. Prepare manuscript draft, supplement, and reproducibility release.

Exit criterion:

- all main figures have source tables
- all target-priority claims have evidence rows
- workflow can be re-run from manifests on a clean environment
