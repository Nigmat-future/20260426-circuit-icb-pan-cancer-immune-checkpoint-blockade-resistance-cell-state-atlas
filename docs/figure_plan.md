# Figure Plan

## Figure 1: Study Design

Input:

- source registry
- discovery and validation manifests

Output:

- schematic of data flow from ICB scRNA-seq to target prioritization
- cohort counts and modality counts

## Figure 2: Malignant-Cell States

Input:

- malignant-cell expression or metacell matrix
- patient metadata
- discovered programs and markers

Output:

- UMAP or low-dimensional map
- state marker heatmap
- pathway enrichment
- response and timepoint association at patient level

## Figure 3: Bulk ICB Validation

Input:

- state signatures
- bulk ICB expression matrices
- harmonized response labels

Output:

- responder versus non-responder state score distributions
- AUROC/AUPRC by cohort
- leave-one-cohort-out performance
- benchmark comparison

## Figure 4: TCGA/CPTAC Characterization

Input:

- TCGA expression and clinical data
- CPTAC protein/phosphoprotein tables

Output:

- pan-cancer state activity landscape
- survival forest plot
- immune context associations
- protein support heatmap

## Figure 5: Target And Drug Priority

Input:

- target evidence matrix
- compound reversal and sensitivity matrix

Output:

- top target table
- multi-evidence radar or tile plot
- dependency and sensitivity panels
- compound-state reversal network

## Figure 6: Spatial And Resource

Input:

- spatial transcriptomics datasets
- state score maps
- top target expression maps

Output:

- spatial state enrichment
- T-cell exclusion or CAF/myeloid niche proximity
- screenshot or schema of resource release
