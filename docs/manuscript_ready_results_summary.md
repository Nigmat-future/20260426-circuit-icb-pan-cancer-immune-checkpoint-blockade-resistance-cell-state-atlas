# Manuscript-Ready Results Summary

> **Scope note.** This document belongs to the earlier pan-cancer `publish_v1` analysis and is not part of the melanoma manuscript (`manuscript/manuscript.md`). See the README for the files behind the manuscript. Here "manuscript" means a planned pan-cancer resource paper, not the submitted melanoma manuscript.

## Recommended Framing

CIRCUIT-ICB should be framed as a reproducible malignant-cell ICB resistance-state atlas and multi-omics target-prioritization resource. The strongest claim is not classifier performance; it is the end-to-end resource linking patient-derived malignant states, external ICB cohorts, TCGA state context, DepMap/PRISM/LINCS evidence, and transparent candidate-target ranking.

## Results To Emphasize

- Six public bulk ICB cohorts were harmonized under a single clinical schema, with 212 binary response/DCB samples.
- Benchmark signatures and random signature controls were added, giving reviewers a clear baseline.
- TCGA Xena scoring covers 9163 primary tumors and provides cancer-stratified state landscape, OS/PFI Cox tables, and immune-proxy correlations.
- The publish_v1 target table contains 33 core resistance-state candidates after explicit down-weighting for common essentiality, broad normal expression risk, and lineage-context programs.
- CPTAC/PDC is handled transparently: an adapter exists, but no tumor protein matrix is filled in this run, and HPA protein existence is not used as a substitute.

## Results To Phrase Conservatively

- The best primary pooled CIRCUIT state AUROC-best-direction is 0.562.
- Leave-one-cohort logistic validation remains modest, with mean AUROC 0.483.
- Random signature null tests do not support a strong classifier claim.

Recommended wording: "CIRCUIT-ICB provides a retrospective, cell-state-resolved resource and candidate mechanism/target prioritization framework." Avoid "clinically actionable predictor" or "ready-to-use diagnostic model."

## Figure Mapping

- Figure 1: data flow and cohort curation, including explicit excluded resources
- Figure 2: NMF state atlas and core versus lineage/context state labeling
- Figure 3: publish_v1 external validation, benchmark comparison, and random signature controls
- Figure 4: TCGA state landscape, survival associations, and immune-proxy correlations
- Figure 5: target/drug evidence tile plot with DepMap, PRISM, LINCS, HPA safety, and CPTAC status
