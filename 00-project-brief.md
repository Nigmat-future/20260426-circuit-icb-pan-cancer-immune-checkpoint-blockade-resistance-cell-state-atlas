# Project Brief

## Working Title

CIRCUIT-ICB: A cell-state-resolved and perturbation-aware pan-cancer atlas of immune checkpoint blockade resistance

## Study Thesis

Immune checkpoint blockade resistance is partly encoded in malignant tumor cell states that recur across cancer types. CIRCUIT-ICB will discover these states from ICB-treated single-cell RNA-seq cohorts, transfer them into bulk clinical cohorts, and prioritize druggable vulnerabilities using proteomic, dependency, drug-screen, normal-tissue, and perturbational evidence.

## Core Contribution

CIRCUIT-ICB is designed as a resource and reproducible framework, not only a predictive signature. The expected contribution is a linked chain:

1. malignant-cell resistance states in ICB-treated patients
2. patient-level state activity scores transferable to bulk RNA-seq
3. cross-cohort and cross-cancer response benchmarks
4. target and compound prioritization with multi-omics evidence
5. an open workflow, release tables, and manuscript-ready figures

## Minimum Publishable Unit

The smallest defensible paper contains:

- a curated discovery atlas from the public ICB-treated scRNA-seq compilation
- 5 to 8 interpretable malignant-cell programs or states
- patient-level association with response, timepoint, and cancer type
- bulk-compatible state scores tested in at least 3 independent ICB response cohorts
- TCGA pan-cancer characterization for prognosis and immune context
- a transparent target-priority table combining expression, dependency, drug, perturbation, safety, and protein evidence

## High-Value Extension

The high-score version adds:

- leave-one-study-out and leave-one-cancer-type-out validation
- calibration, decision-curve analysis, and benchmark against published immune signatures
- CPTAC protein and phosphoprotein confirmation
- PRISM/GDSC sensitivity association in high-state cell lines
- LINCS/CLUE perturbation reversal analysis
- optional spatial transcriptomics validation of T-cell exclusion or suppressive niches

## Success Conditions

- All analyses are reproducible from manifests and scripts.
- All patient-response tests use patient-level units, not cell-level pseudo-replication.
- State scores improve or complement standard immune response markers in external cohorts.
- Target prioritization is multi-evidence and does not rely only on differential expression.
- Claims are framed as retrospective biomarker and hypothesis-generation evidence.

## Top Risks

- Response labels differ across studies and can create endpoint noise.
- Single-cell discovery cohorts are heterogeneous in cancer type, platform, treatment timepoint, and treatment regimen.
- Bulk transfer can capture tumor purity or immune infiltration rather than malignant-intrinsic programs.
- Drug perturbation reversal does not prove clinical synergy.
- Large RDS and proteomic files may require staged download and local storage planning.
