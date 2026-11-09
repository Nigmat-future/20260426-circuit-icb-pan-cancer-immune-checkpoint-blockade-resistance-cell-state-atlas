# Data Modality And QC

## Data Modalities

- ICB-treated scRNA-seq for malignant-cell state discovery
- public tumor scRNA-seq resources for cell-type specificity checks
- bulk RNA-seq and microarray ICB cohorts for response validation
- TCGA/GDC for pan-cancer survival, immune context, purity, and genomic background
- CPTAC/PDC for protein and phosphoprotein support
- HPA and GTEx for normal-tissue safety filtering
- DepMap CRISPR and RNAi for dependency evidence
- PRISM and GDSC for drug-sensitivity evidence
- LINCS/CLUE L1000 for perturbational reversal
- optional spatial transcriptomics for niche-level validation

## scRNA-seq QC Plan

Initial implementation uses curated Seurat objects and pseudobulk files from the Zenodo resource to avoid reprocessing all raw count matrices. A later strict-reproducibility branch can rebuild objects from original sources.

QC checks:

- verify sample, patient, cancer type, treatment timing, and response columns
- check malignant-cell annotation source and consistency
- filter genes and cells only with documented thresholds
- quantify per-study cell counts and patient counts
- aggregate all inferential tests to patient or sample level
- compare malignant-state markers against immune/stromal markers to detect contamination

## Malignant-Cell Identification

Preferred order:

1. use curated malignant-cell annotations from the source resource
2. confirm with epithelial or lineage markers and tumor-specific markers
3. when raw data are available, run CopyKAT or inferCNV as sensitivity analysis
4. exclude ambiguous cells from primary state discovery and use them only in sensitivity checks

## State Discovery

Candidate methods:

- non-negative matrix factorization on malignant-cell pseudobulk or metacells
- consensus NMF within study and across merged feature space
- topic modeling or ICA as sensitivity analysis
- mixed-effect modeling for responder versus non-responder associations

Program labels are assigned by Hallmark, Reactome, GO, TF activity, and known immune-resistance literature. Labels remain descriptive, not causal.

## Bulk Transfer QC

- map all expression matrices to HGNC gene symbols
- drop duplicated genes by median expression or maximum variance rule
- use rank-based scoring to reduce platform and normalization sensitivity
- evaluate tumor purity and immune infiltration correlation
- include housekeeping and random signatures as negative controls
- report missing signature-gene fraction per sample

## Target Evidence QC

For each target, store a directional evidence table:

- state specificity in malignant cells
- clinical association in ICB cohorts
- TCGA immune/prognostic association
- CPTAC protein support
- DepMap dependency
- drug sensitivity in high-state cell lines
- LINCS perturbation reversal
- GTEx/HPA normal-tissue safety
- druggability annotation

No single evidence layer is treated as sufficient.
