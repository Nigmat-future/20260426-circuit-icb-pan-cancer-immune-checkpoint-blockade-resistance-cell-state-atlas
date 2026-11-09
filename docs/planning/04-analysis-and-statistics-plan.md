# Analysis And Statistics Plan

## Aim 1: Malignant-Cell State Atlas

Inputs:

- curated ICB-treated scRNA-seq objects or pseudobulk data
- patient/sample-level metadata
- malignant-cell annotation

Primary analysis:

- discover malignant programs with consensus NMF or metacell-level matrix factorization
- label programs with pathway and TF enrichment
- aggregate cell-state activity to patient/sample level
- test association with response using patient-level models

Statistical model:

- continuous state activity: linear mixed model or Wilcoxon at patient level
- binary response: logistic regression adjusted for cancer type and study when sample size permits
- timepoint analysis: paired or mixed model for pre/post samples

## Aim 2: Bulk Transfer And ICB Prediction

Inputs:

- state marker sets from Aim 1
- bulk ICB expression matrices
- response and survival metadata

Primary score:

- rank-based up-minus-down state activity
- report missing-gene fraction and score confidence

Predictive model:

- univariable state score models
- low-parameter elastic net logistic regression for response
- Cox model for PFS/OS when follow-up is available

Validation:

- leave-one-cohort-out
- leave-one-cancer-type-out
- bootstrap confidence intervals
- calibration and decision-curve analysis

## Aim 3: Target And Drug Prioritization

Target priority score:

`state_specificity + clinical_association + protein_evidence + dependency + drug_sensitivity + perturbation_reversal + safety + druggability`

Each evidence component is stored separately and normalized to 0 to 1 before weighting. Missing evidence is not imputed as positive. Sensitivity analyses compare equal weighting, conservative weighting, and evidence-completeness weighting.

## Aim 4: Optional Spatial Validation

Questions:

- are resistant malignant states spatially enriched in T-cell excluded niches?
- are they adjacent to CAF, myeloid-suppressive, hypoxic, or ECM-rich regions?
- are top target genes expressed in invasive fronts or immune-excluded regions?

Spatial analyses are orthogonal validation only unless the dataset includes ICB response labels.

## Baseline Comparators

- T-cell inflamed gene expression profile
- cytolytic score, such as GZMA/PRF1
- IFN-gamma response score
- antigen-presentation score
- cell-cycle score
- EMT/TGF-beta score
- TIDE or TIDE-like published scores where available
- tumor purity and immune infiltration proxies
- random same-size gene sets
- housekeeping gene sets

## Leakage Prevention

- derive state markers only from discovery data
- do not tune marker thresholds on external validation response labels
- do not split cells from the same patient across training and validation
- nested resampling for any hyperparameter selection
- keep TCGA pan-cancer survival analysis separate from ICB response model fitting

## Minimum Decisive Experiment

The project remains worth pursuing if at least two malignant-cell states:

- are reproducible across discovery studies,
- have patient-level response association,
- transfer to independent bulk ICB cohorts,
- remain informative after adjusting for cancer type, tumor purity proxy, and T-cell inflamed score,
- nominate biologically coherent targets with at least two independent non-expression evidence layers.
