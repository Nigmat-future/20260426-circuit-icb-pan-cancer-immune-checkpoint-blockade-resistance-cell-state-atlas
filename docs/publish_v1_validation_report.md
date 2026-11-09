# Publish V1 Validation Report

> **Scope note.** This document belongs to the earlier pan-cancer `publish_v1` analysis and is not part of the melanoma manuscript (`manuscript/manuscript.md`). See the README for the files behind the manuscript.

Generated on 2026-04-26.

## Cohort Gate

- Included cohorts: 6
- Binary response/DCB samples: 212
- Minimum gate: met (`>=4` cohorts and `>=150` binary samples)

Included cohorts are GSE78220, GSE91061, GSE115821, GSE176307, GSE135222, and GSE126044. IMvigor210 and Gide2019 are retained in the curation log with explicit non-inclusion reasons.

## Primary CIRCUIT States

Primary states for main-text validation:

- `CIRCUIT_STATE_01_score`: IFNG / antigen-presentation response-context state
- `CIRCUIT_STATE_02_score`: core hypoxia-stress resistance state
- `CIRCUIT_STATE_04_score`: core EMT/TGF-beta exclusion state
- `CIRCUIT_STATE_06_score`: core hypoxia-stress resistance state

Pooled AUROC-best-direction values are modest:

- `CIRCUIT_STATE_04_score`: 0.562
- `CIRCUIT_STATE_06_score`: 0.556
- `CIRCUIT_STATE_02_score`: 0.553
- `CIRCUIT_STATE_01_score`: 0.552

The 1000x size-matched random signature tests are not significant for the four primary states. Leave-one-cohort logistic validation has mean AUROC 0.483.

## Manuscript Boundary

The publish_v1 external validation is adequate for a public-data resource manuscript, but not for a clinical classifier claim. Main text should say that CIRCUIT-ICB maps reproducible malignant-cell states and nominates candidate target/drug hypotheses. Predictive performance should be presented as retrospective association evidence and benchmark context.

## Key Tables

- `results/tables/publish_v1/cohort_summary.tsv`
- `results/tables/publish_v1/validation/publish_v1_univariable_ci.tsv`
- `results/tables/publish_v1/validation/publish_v1_random_effects_meta.tsv`
- `results/tables/publish_v1/validation/publish_v1_leave_one_cohort.tsv`
- `results/tables/publish_v1/validation/publish_v1_leave_one_cancer_type.tsv`
- `results/tables/publish_v1/validation/publish_v1_random_signature_pvalues.tsv`
