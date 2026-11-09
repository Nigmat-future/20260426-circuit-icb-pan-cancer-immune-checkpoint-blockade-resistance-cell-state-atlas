# Validation Plan

## Internal Single-Cell Validation

- leave-one-study-out state discovery and marker stability
- leave-one-cancer-type-out state recovery
- patient-level bootstrap stability for state-response association
- compare malignant-only states with all-cell states to identify contamination
- verify that cell-cycle-only programs do not explain all response signal

## External Bulk ICB Validation

For every cohort:

- harmonize genes and response labels
- compute state scores without refitting marker sets
- evaluate response association with AUROC and AUPRC
- fit simple adjusted models when covariates are available
- report cohort-specific and pooled random-effects estimates

## Pan-Cancer TCGA Validation

TCGA is not an ICB response cohort. It is used for:

- prognosis association
- immune subtype context
- tumor purity and immune infiltration correlation
- TMB-like proxy, CNV burden, and antigen-presentation context
- cancer-type heterogeneity and specificity

TCGA results are interpreted as biological characterization, not direct ICB prediction.

## Orthogonal Target Validation

For candidate targets:

- CPTAC protein abundance supports transcript-level state association
- DepMap dependency is stronger in relevant lineage or high-state cell lines
- PRISM or GDSC compounds show stronger sensitivity in high-state models
- LINCS perturbagens reverse the state signature direction
- GTEx/HPA normal tissue expression flags safety concerns

## Robustness Checks

- random gene-set score distribution
- housekeeping gene-set score distribution
- alternative marker sizes per state
- alternative rank-score definitions
- exclusion of one large cancer type or one large study
- pretreatment-only analysis
- response-label sensitivity: exclude SD, count SD as non-response, or keep three-class models

## Failure Analysis

If validation fails, record:

- which state failed
- whether failure is cohort-specific or cancer-type-specific
- whether gene coverage was insufficient
- whether the state captured immune infiltration, tumor purity, or platform effects
- whether a state should be reclassified as cancer-specific rather than pan-cancer
