# External Validation And Target Evidence

> **Scope note.** This document belongs to the earlier pan-cancer `publish_v1` analysis and is not part of the melanoma manuscript (`manuscript/manuscript.md`). See the README for the files behind the manuscript.

## Locked External Bulk ICB Cohorts

Publish_v1 locked validation cohorts:

- GSE78220: melanoma anti-PD-1, processed FPKM with response labels from GEO metadata
- GSE91061: melanoma anti-PD-1, pretreatment-only primary validation, Entrez IDs mapped through NCBI Homo_sapiens.gene_info
- GSE115821 / Liu2019: melanoma ICB, featureCounts converted to log2(CPM+1), pretreatment samples collapsed to patient level
- GSE176307: urothelial carcinoma ICB, Salmon gene TPM converted to log2(TPM+1)
- GSE135222: NSCLC anti-PD-1/anti-PD-L1, Ensembl IDs mapped through HGNC and DCB/NDB derived from PFS >= 180 days
- GSE126044: NSCLC anti-PD-1, counts converted to log2(CPM+1)

Audited but not included:

- IMvigor210: official R package downloaded and inspected, but `cds.RData` is a Bioconductor S4 CountDataSet unsupported by pyreadr/rdata here and no Rscript runtime is available.
- Gide2019 / PRJEB23709: processed expression table not identified; raw FASTQ/SRA processing is outside the locked assumptions.

Current response validation is deliberately simple:

- score NMF state signatures in bulk expression
- harmonize response labels
- evaluate univariable AUROC/AUPRC per state
- run bootstrap CI, random-effects meta-analysis, leave-one-cohort, leave-one-cancer-type, calibration, benchmark scoring, and 1000x random-signature null tests

## Current External Result

The publish_v1 registry contains 6 included cohorts and 212 binary response/DCB samples:

- best primary pooled CIRCUIT state: `CIRCUIT_STATE_04_score`, AUROC-best-direction 0.562
- 1000x size-matched random signature p-values for primary states are not significant
- mean leave-one-cohort logistic AUROC is 0.483
- conclusion: use state-association/resource language, not clinical classifier language

## Target Evidence Matrix

The current matrix includes:

- `state_specificity`: from malignant-only NMF state marker weights
- `clinical_association`: gene-level association in GSE78220 and GSE91061
- `tcga_prognosis`: adverse prognostic fraction from HPA proteinatlas TCGA prognostic summaries
- `dependency`: DepMap Public 26Q1 CRISPRGeneDependency, with common-essential genes capped conservatively
- `drug_sensitivity`: PRISM secondary screen target-compound median AUC converted to a rank score
- `perturbation_reversal`: L1000CDS2 core-state reversal perturbations bridged to genes through PRISM target annotations
- `safety`: HPA RNA tissue consensus normal-expression burden
- `druggability`: HPA protein-class heuristic, including FDA target, enzyme, transporter, membrane, and protein-evidence classes

Reserved column still present for:

- `protein_evidence`: CPTAC/PDC

Use `scripts/20_build_cptac_protein_evidence.py` to adapt a supplied CPTAC/PDC tumor protein abundance table into the `protein_evidence` column. We deliberately do not use HPA protein existence as a substitute for CPTAC tumor proteomics.

Publish_v1 also includes:

- `results/tables/publish_v1/tcga/tcga_state_scores.tsv`: TCGA Xena z-score state scores for 9163 primary tumors
- `results/tables/publish_v1/tcga/tcga_state_survival_os.tsv`: OS Cox associations
- `results/tables/publish_v1/tcga/tcga_state_survival_pfi.tsv`: PFI Cox associations
- `results/tables/publish_v1/tcga/tcga_state_immune_proxy_correlations.tsv`: IFNG/T-cell/cytolytic immune-proxy correlations
- `results/tables/publish_v1/cptac/cptac_protein_evidence_status.tsv`: CPTAC/PDC status and runtime blocker record

The current merged target tables are:

- `results/tables/target_evidence/circuit_icb_target_evidence.merged.tsv`
- `results/tables/target_evidence/circuit_icb_target_priority.merged.tsv`
- `results/tables/target_evidence/circuit_icb_target_priority.annotated.tsv`
- `results/tables/target_evidence/core_resistance_target_priority.tsv`
- `results/tables/publish_v1/target_evidence/publish_v1_target_priority_all.tsv`
- `results/tables/publish_v1/target_evidence/publish_v1_core_resistance_target_priority.tsv`

## Interpretation

The current target priority table is a partially populated but runnable public-data evidence matrix. Scores should be interpreted as prioritization hypotheses. `BREAST_LUMINAL`, `LIVER_SECRETORY_COMPLEMENT`, and `EPITHELIAL_KERATIN` are explicitly marked as lineage-context states, so the primary resistance table focuses on `HYPOXIA_STRESS` and `EMT_TGFB_EXCLUSION` candidates. Common-essential genes, low HPA safety scores, and lineage-context states are explicitly down-weighted in publish_v1 target tables.
