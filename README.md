# CIRCUIT-ICB

CIRCUIT-ICB is a reproducible research framework for discovering malignant-cell immune checkpoint blockade resistance states and linking them to clinical response, protein evidence, functional dependency, drug sensitivity, and perturbational reversal.

## Working Hypothesis

ICB failure is associated with recurrent malignant-cell-intrinsic states that can be learned from patient scRNA-seq, transferred to bulk RNA-seq, and connected to druggable vulnerabilities through public multi-omics perturbation resources.

## Repository Layout

```text
config/                 Project configuration
data/manifests/         Curated source and cohort manifests
data/raw/               Downloaded public metadata or raw files
data/processed/         Derived analysis-ready tables
docs/                   Method notes and execution documentation
notebooks/              Exploratory and figure notebooks
results/                Analysis outputs
scripts/                CLI wrappers for pipeline steps
src/circuit_icb/        Python package
tests/                  Unit tests for scoring and prioritization logic
workflow/               Snakemake workflow
```

## Current Implementation Status

- study design pack: created
- Zenodo ICB scRNA-seq metadata workbook: downloaded
- manifest schema: included
- rank-based state scoring code: included
- target prioritization code: included
- workflow skeleton: included
- Zenodo pseudobulk package: downloaded and scored with example state signatures
- malignant-only NMF states: extracted from public malignant pseudobulk matrices
- external bulk ICB validation: GSE78220 and GSE91061 prepared, scored, and evaluated
- publish_v1 external validation: 6 included public bulk ICB cohorts, 212 binary response/DCB samples, benchmark signatures, bootstrap CI, random-effects meta-analysis, leave-one-cohort and leave-one-cancer validation, calibration, and 1000x random-signature null tests
- malignant state annotation: NMF states labeled into IFNG/antigen presentation, hypoxia-stress, EMT/TGF-beta exclusion, and lineage-context programs
- target evidence matrix: state specificity, external ICB clinical association, HPA/TCGA prognosis, DepMap dependency, PRISM drug sensitivity, LINCS reversal bridge, HPA safety, and HPA-derived druggability populated
- TCGA Xena state landscape: 9163 primary tumors scored with z-score state signatures, plus OS/PFI Cox and immune-proxy correlation tables
- CPTAC/PDC protein evidence: adapter-ready with explicit status table; not populated unless a transparent tumor protein abundance table is supplied
- full biological results: current validation is a first-pass public-data benchmark, not a final clinical biomarker claim

## Quick Start

```bash
python -m pip install -e ".[dev]"
python scripts/00_build_icb_manifest.py
python scripts/01_score_bulk_states.py --expression examples/bulk_expression_tpm.tsv --signatures config/state_signatures.example.yml --out results/tables/example_state_scores.tsv
python scripts/02_prioritize_targets.py --evidence examples/target_evidence.tsv --weights config/target_weights.yml --out results/tables/example_target_priority.tsv
python scripts/04_score_zenodo_pseudobulk.py
python scripts/05_extract_nmf_states.py
python scripts/07_prepare_geo_bulk.py --cohort GSE78220
python scripts/07_prepare_geo_bulk.py --cohort GSE91061
python scripts/07_prepare_geo_bulk.py --cohort GSE115821
python scripts/07_prepare_geo_bulk.py --cohort GSE176307
python scripts/07_prepare_geo_bulk.py --cohort GSE135222
python scripts/07_prepare_geo_bulk.py --cohort GSE126044
python scripts/08_combine_locked_external_bulk.py
python scripts/09_build_target_evidence_matrix.py
python scripts/12_build_hpa_safety.py
python scripts/14_build_depmap_dependency.py
python scripts/16_build_prism_sensitivity.py
python scripts/17_query_l1000cds2_reversal.py
python scripts/18_build_lincs_gene_reversal.py
python scripts/19_build_hpa_atlas_evidence.py
python scripts/10_merge_external_evidence.py --tcga results/tables/target_evidence/hpa_atlas_tcga_druggability.tsv --depmap results/tables/target_evidence/depmap_dependency.tsv --prism results/tables/target_evidence/prism_drug_sensitivity.tsv --lincs results/tables/target_evidence/lincs_gene_reversal.tsv --safety results/tables/target_evidence/hpa_safety.tsv
python scripts/15_annotate_target_priority.py
python scripts/21_build_publish_v1_cohorts.py
python scripts/22_score_publish_v1_cohorts.py
python scripts/23_validate_publish_v1.py
python scripts/24_build_publish_v1_target_priority.py
python scripts/26_build_tcga_xena_state_landscape.py
python scripts/27_write_cptac_status.py
python scripts/25_build_publish_v1_manuscript_summary.py
pytest
```

## Primary Data Anchor

The discovery atlas is anchored to Zenodo record `14511579`, version 4 of "Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint blockade-treated patients".

The current smoke test also downloads and scores the public pseudobulk package from the same record when present locally.

## Claim Boundary

CIRCUIT-ICB is a retrospective computational framework. Candidate targets and compounds are hypotheses for follow-up, not clinical recommendations.

The publish_v1 result set meets the minimum external-validation sample gate, but the leave-one-cohort classifier remains modest. Manuscript language should emphasize a malignant-cell state resource and target-prioritization framework, not a clinical-use predictor.
