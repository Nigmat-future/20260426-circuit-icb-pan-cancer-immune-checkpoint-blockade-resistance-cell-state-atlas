# Current Outputs

> **Scope note.** This document belongs to the earlier pan-cancer `publish_v1` analysis and is not part of the melanoma manuscript (`manuscript/manuscript.md`). See the README for the files behind the manuscript.

Generated on 2026-04-26.

## Metadata And Manifests

- `data/raw/zenodo/v4_ICB_datasets.xlsx`: public Zenodo workbook
- `data/raw/zenodo/zenodo_14511579_files.tsv`: Zenodo file index with sizes, checksums, and download URLs
- `data/raw/zenodo/zenodo_10407126_files.tsv`: malignant-only Zenodo v1 file index
- `data/manifests/icb_sc_discovery_datasets.tsv`: normalized discovery dataset manifest
- `data/manifests/icb_sc_discovery_files.tsv`: normalized file-location manifest
- `data/manifests/icb_sc_treatment_information.tsv`: normalized treatment metadata
- `data/manifests/zenodo_pseudobulk_profile.tsv`: profile of extracted pseudobulk matrices
- `data/manifests/bulk_icb_candidate_cohorts.csv`: locked and candidate external bulk cohorts
- `data/manifests/target_evidence_sources.csv`: target evidence layer registry

## Smoke-Test Results

- `results/tables/example_state_scores.tsv`: rank-based state scores for the toy bulk matrix
- `results/tables/example_target_priority.tsv`: weighted target priority score for the toy evidence table
- `results/tables/zenodo_pseudobulk_state_scores.tsv`: state scores for 935 rows from 11 public pseudobulk matrices
- `results/tables/release_index.tsv`: compact index of generated result tables

## Main Analysis Outputs

- `results/tables/nmf_states/nmf_state_activity.tsv`: 7 malignant-only NMF state activities for 174 pseudobulk rows
- `results/tables/nmf_states/nmf_state_markers.tsv`: top 50 genes for each NMF state
- `results/tables/nmf_states/nmf_state_signatures.yml`: bulk-compatible signature file generated from NMF markers
- `results/tables/bulk_icb/GSE78220.state_scores.tsv`: NMF state scores in GSE78220
- `results/tables/bulk_icb/GSE91061.state_scores.tsv`: NMF state scores in pretreatment GSE91061
- `results/tables/bulk_icb/locked_external_bulk_state_scores.tsv`: combined locked external bulk validation table
- `results/tables/validation/locked_external_bulk_nmf.univariable.tsv`: per-state external AUROC/AUPRC table
- `results/tables/validation/locked_external_bulk_nmf.leave_one_cohort.tsv`: leave-one-cohort external validation table
- `results/tables/target_evidence/circuit_icb_target_evidence.tsv`: current target evidence matrix
- `results/tables/target_evidence/circuit_icb_target_priority.tsv`: current target priority table
- `results/tables/target_evidence/hpa_safety.tsv`: HPA normal tissue expression safety evidence
- `results/tables/target_evidence/hpa_atlas_tcga_druggability.tsv`: HPA proteinatlas TCGA prognostic and druggability evidence
- `results/tables/target_evidence/depmap_dependency.tsv`: DepMap Public 26Q1 CRISPR dependency evidence
- `results/tables/target_evidence/prism_drug_sensitivity.tsv`: PRISM secondary-screen target-compound sensitivity evidence
- `results/tables/target_evidence/lincs_state_reversal.tsv`: L1000CDS2 state-reversal perturbation table
- `results/tables/target_evidence/lincs_gene_reversal.tsv`: LINCS-to-gene reversal bridge through PRISM target annotations
- `results/tables/target_evidence/circuit_icb_target_evidence.merged.tsv`: merged multi-omics evidence matrix
- `results/tables/target_evidence/circuit_icb_target_priority.annotated.tsv`: target priority table with state relevance labels
- `results/tables/target_evidence/core_resistance_target_priority.tsv`: filtered target table for hypoxia-stress and EMT/TGF-beta exclusion states

## Publish V1 Outputs

- `data/manifests/publish_v1_bulk_icb_cohorts.tsv`: locked publish_v1 external bulk ICB registry
- `results/tables/publish_v1/cohort_curation_log.tsv`: included/excluded cohort audit trail, including IMvigor210 and Gide2019 failure reasons
- `results/tables/publish_v1/cohort_summary.tsv`: 6 included cohorts and 212 binary response/DCB samples
- `results/tables/publish_v1/publish_v1_all_scores.tsv`: CIRCUIT state and benchmark signature scores for 216 bulk ICB samples
- `results/tables/publish_v1/validation/publish_v1_univariable_ci.tsv`: per-cohort and pooled AUROC/AUPRC with bootstrap CI
- `results/tables/publish_v1/validation/publish_v1_random_effects_meta.tsv`: random-effects meta-analysis across cohorts
- `results/tables/publish_v1/validation/publish_v1_leave_one_cohort.tsv`: leave-one-cohort logistic validation
- `results/tables/publish_v1/validation/publish_v1_leave_one_cancer_type.tsv`: leave-one-cancer-type validation
- `results/tables/publish_v1/validation/publish_v1_random_signature_pvalues.tsv`: 1000x size-matched random signature null tests for primary CIRCUIT states
- `results/tables/publish_v1/target_evidence/publish_v1_core_resistance_target_priority.tsv`: main-text target table after common-essential, safety, and lineage-context down-weighting
- `results/tables/publish_v1/tcga/tcga_state_scores.tsv`: TCGA Xena z-score state landscape for 9163 primary tumors
- `results/tables/publish_v1/tcga/tcga_state_survival_os.tsv`: TCGA OS Cox associations
- `results/tables/publish_v1/tcga/tcga_state_survival_pfi.tsv`: TCGA PFI Cox associations
- `results/tables/publish_v1/tcga/tcga_state_immune_proxy_correlations.tsv`: TCGA state correlations with IFNG/T-cell/cytolytic proxy signatures
- `results/tables/publish_v1/cptac/cptac_protein_evidence_status.tsv`: CPTAC/PDC adapter and runtime status; HPA is not substituted for CPTAC
- `results/tables/publish_v1/manuscript_ready_summary.tsv`: compact manuscript claim-boundary summary

## Verification

- Unit tests: 32 passed with `pytest -q`
- Compile check: passed with `python -m compileall -q src scripts`
- Full Zenodo pseudobulk scoring completed locally
- Malignant-only NMF state extraction completed locally
- GSE78220, GSE91061, GSE115821, GSE176307, GSE135222, and GSE126044 publish_v1 scoring and validation completed locally
- TCGA Xena API-based state landscape completed locally
- Snakemake was not executed because the local environment did not have the `snakemake` command installed; `environment.yml` includes `snakemake-minimal` for a reproducible environment

## Interpretation Caution

The publish_v1 external-validation gate is met by cohort/sample count, but predictive performance remains modest. The manuscript should frame CIRCUIT-ICB as a reproducible malignant-cell state resource and target-prioritization framework, not a clinical-use classifier. The target matrix has public-data evidence from HPA/TCGA summary fields, TCGA Xena state landscape, DepMap, PRISM, LINCS, and HPA safety; CPTAC/PDC protein abundance remains adapter-ready but not populated.
