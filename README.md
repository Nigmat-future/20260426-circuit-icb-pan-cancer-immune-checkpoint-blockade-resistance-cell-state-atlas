<div align="center">
  
# 🧬 CIRCUIT-ICB

**C**omputational **I**nference of **R**esistance **C**ell-states **U**sing **I**ntegrated **T**ranscriptomics in **I**mmune **C**heckpoint **B**lockade

[![Python](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/Workflow-Snakemake-brightgreen.svg)](https://snakemake.readthedocs.io/)
[![Zenodo](https://img.shields.io/badge/Data_Anchor-Zenodo_14511579-blue.svg)](https://zenodo.org/record/14511579)

*A reproducible, scalable research framework for discovering malignant-cell immune checkpoint blockade (ICB) resistance states and linking them to clinical response, multi-omic evidence, and druggable vulnerabilities across pan-cancer datasets.*

</div>

---

## 🎯 Working Hypothesis

ICB therapy failure is strongly associated with recurrent **malignant-cell-intrinsic states**. By learning these states from patient-derived single-cell RNA-seq (scRNA-seq), transferring them to bulk RNA-seq cohorts, and mapping them to expansive molecular networks, **CIRCUIT-ICB** connects transcriptional resistance signatures directly to **druggable vulnerabilities** via comprehensive evidence synthesis.

## 🚀 Key Features & Implementation Status

- **Robust State Inference:** Extracts malignant-only NMF states from public pseudobulk matrices and annotates them into programs (e.g., IFNG/antigen presentation, hypoxia-stress, EMT/TGF-beta exclusion, lineage-context).
- **Extensive External Validation:** Prepares, scores, and evaluates rigorous external bulk ICB validation across 6 cohorts (GSE78220, GSE91061, GSE115821, GSE176307, GSE135222, GSE126044) comprising 212 binary response/DCB samples. Includes bootstrap CI, random-effects meta-analysis, and leave-one-cohort cross-validation.
- **Multimodal Target Prioritization:** Synthesizes a comprehensive evidence matrix integrating:
  - State specificity & External ICB clinical association
  - HPA/TCGA prognosis profiles
  - DepMap genetic dependencies
  - PRISM drug sensitivity panels
  - LINCS L1000 reversal bridge profiling
  - HPA safety and druggability mappings
- **Large-Scale Multi-Omic Mapping:** Scores 9,163 primary tumors from the TCGA Xena landscape using z-score state signatures, coupled with OS/PFI Cox regressions and immune-proxy correlations.
- **Proteomic Adapter Readiness:** Built-in hooks for CPTAC/PDC protein evidence integration.

## 📂 Repository Architecture

Structured for scalability, reproducibility, and clarity:

```text
├── config/              # Project configuration & state signature YAMLs
├── data/
│   ├── manifests/       # Curated source and cohort manifests
│   ├── raw/             # Downloaded public metadata and raw files
│   └── processed/       # Derived analysis-ready matrices and tables
├── docs/                # Methodological notes and execution documentation
├── notebooks/           # Jupyter notebooks for exploratory data analysis (EDA) and figure generation
├── results/             # Final analysis outputs, evaluation metrics, and plots
├── scripts/             # Python CLI wrappers for individual pipeline steps
├── src/circuit_icb/     # Core Python package (scoring, prioritization, inference logic)
├── tests/               # Pytest-based unit testing for scoring and prioritization logic
└── workflow/            # Snakemake workflow orchestrator
```

## 💻 Quick Start

Set up your environment and execute the pipeline via the unified CLI scripts:

```bash
# 1. Install the core package in development mode
python -m pip install -e ".[dev]"

# 2. Build the primary data manifests
python scripts/00_build_icb_manifest.py

# 3. State Scoring & Target Prioritization (Example Run)
python scripts/01_score_bulk_states.py \
    --expression examples/bulk_expression_tpm.tsv \
    --signatures config/state_signatures.example.yml \
    --out results/tables/example_state_scores.tsv

python scripts/02_prioritize_targets.py \
    --evidence examples/target_evidence.tsv \
    --weights config/target_weights.yml \
    --out results/tables/example_target_priority.tsv

# 4. End-to-End Discovery Pipeline
python scripts/04_score_zenodo_pseudobulk.py
python scripts/05_extract_nmf_states.py

# ... cohort preparation, matrix building, evidence mapping, and evaluation ...
# (Refer to standard operating procedures in docs/ for full pipeline execution, including scripts 07-27)

# 5. Run test suite to ensure computational integrity
pytest
```

## 📊 Primary Data Anchor

The discovery atlas is deeply anchored to **Zenodo record `14511579`** (Version 4): *"Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint blockade-treated patients"*. 

> **Note:** The current smoke test automatically downloads and scores the public pseudobulk package from this Zenodo record when present locally.

## ⚠️ Claim Boundary & Disclaimer

**CIRCUIT-ICB is a retrospective computational research framework.** Candidate targets, states, and compounds generated by this pipeline represent robust, data-driven **hypotheses for downstream experimental follow-up**, and do not constitute clinical recommendations.

*While the `publish_v1` result set meets minimum external-validation sample gates, the leave-one-cohort classifier performance remains modest. Publications utilizing this framework should emphasize its value as a malignant-cell state resource and hypothesis-generation engine rather than a definitive clinical biomarker claim.*

---
<div align="center">
  <i>Built for the future of precision immuno-oncology.</i>
</div>
