# Cohort And Phenotype Plan

## Discovery Cohort

Primary discovery uses the public Zenodo resource "Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint blockade-treated patients" version 4, DOI `10.5281/zenodo.14511579`.

The version 4 metadata reports 8 scRNA-seq datasets, 9 cancer types, 223 patients, 90,270 cancer cells, and 265,671 other cells. The project stores the public metadata workbook at `data/raw/zenodo/v4_ICB_datasets.xlsx`.

## Discovery Inclusion Criteria

- human tumor scRNA-seq data from patients exposed to immune checkpoint blockade
- malignant cells available or recoverable from curated annotations
- patient/sample metadata contains at least cancer type, sample identifier, treatment timing, and response-related labels
- response label can be mapped to a common controlled vocabulary
- sample provenance is traceable to the original source

## Discovery Exclusion Criteria

- no patient or sample identifier
- no malignant-cell annotation and no feasible malignant-cell recovery
- response label cannot be mapped or audited
- duplicate samples without resolvable provenance
- cell-level metadata lacking cancer type or study label

## Bulk Validation Cohorts

Bulk ICB cohorts are split into:

- response cohorts with CR/PR/SD/PD or durable clinical benefit labels
- survival cohorts with PFS or OS after ICB
- cancer-specific cohorts for melanoma, urothelial cancer, kidney cancer, lung cancer, breast cancer, liver cancer, and other available tumor types
- negative-control non-ICB or TCGA cohorts for specificity analysis

Seed examples to curate include GSE78220, GSE91061, PRJEB23709/ERP105482, GSE115821, and IMvigor210-derived resources where legally available.

## Endpoint Hierarchy

Primary endpoint:

- objective response: CR/PR versus SD/PD

Secondary endpoints:

- durable clinical benefit, using the original study definition when available
- PFS and OS
- response at pretreatment samples only
- acquired resistance when pre/post samples exist
- immune-hot versus immune-cold proxy in TCGA

## Harmonized Response Vocabulary

- `responder`: CR, PR, durable clinical benefit, original responder label
- `non_responder`: PD, primary resistance, original non-responder label
- `stable_or_mixed`: SD, mixed response, indeterminate benefit
- `unknown`: unavailable, not assessable, or not harmonizable

Primary binary models use `responder` versus `non_responder`. Sensitivity analyses evaluate whether SD should be included with non-response or excluded.

## Confounders

- cancer type and histology
- treatment class, including anti-PD-1, anti-PD-L1, anti-CTLA-4, and combinations
- treatment timing: pretreatment, on-treatment, post-treatment, resistant lesion
- prior therapy and treatment line
- tumor site: primary versus metastasis
- sequencing platform and library protocol
- tumor purity and malignant-cell fraction
- TMB proxy, PD-L1 proxy, and interferon/T-cell inflamed signatures

## Leakage Controls

- no marker selection using external validation labels
- no patient split after cell-level processing that mixes the same patient across train and test
- leave-one-study-out validation for single-cell state robustness
- leave-one-cohort-out validation for bulk ICB prediction
- leave-one-cancer-type-out validation for pan-cancer claims
