# Cover Letter Summary

Dear Editor,

We submit the manuscript "CIRCUIT-ICB: NMF-based discovery of malignant cell states predicting immune checkpoint blockade resistance in melanoma" for consideration in Bioinformatics.

This study presents CIRCUIT-ICB, a reproducible computational framework for discovering malignant-cell-intrinsic transcriptional states from scRNA-seq pseudobulk data and transferring them into bulk RNA-seq immunotherapy response prediction. The framework addresses a practical failure mode of unsupervised decomposition in this setting: signatures dominated by housekeeping and non-specific genes that do not outperform random gene sets.

Using melanoma as a focused use case, we combined extensive housekeeping-gene exclusion, specificity-based marker filtering, stability-guided NMF component selection, lineage-state detection, and rank-based bulk scoring. Across 43 malignant-cell pseudobulk samples from three melanoma scRNA-seq datasets, CIRCUIT-ICB identified six stable transcriptional states. In two independent melanoma bulk RNA-seq ICB cohorts comprising 76 evaluable patients, two resistance-associated states reached pooled AUROC 0.729 and outperformed matched random signatures by permutation testing (P < 0.05). A multi-state elastic-net model achieved cross-validated AUROC 0.694, exceeding benchmark transcriptomic signatures evaluated in the same framework.

We believe the manuscript will interest Bioinformatics readers because it contributes both a practical analytical workflow and a biologically interpretable melanoma-specific application. The work is method-centered, fully reproducible, and directly addresses the translational challenge of linking single-cell transcriptional programs to clinically accessible bulk RNA-seq predictors.

This manuscript is original, is not under consideration elsewhere, and all authors have approved its submission.
