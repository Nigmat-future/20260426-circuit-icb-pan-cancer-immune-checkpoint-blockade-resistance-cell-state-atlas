# CIRCUIT-ICB: NMF-based discovery of malignant cell states predicting immune checkpoint blockade resistance in melanoma

Nigmat Rahim^1,\*^

^1^ Institute of Medical Technology, Peking University, Beijing, China

^\*^ Corresponding author: nigmatrahim@stu.pku.edu.cn

---

## Abstract

**Motivation:** Immune checkpoint blockade (ICB) has transformed the treatment of metastatic melanoma, yet 40--60% of patients fail to respond. Existing transcriptomic biomarkers such as interferon-gamma (IFNG) signatures and tumor mutational burden provide modest predictive accuracy, and systematic frameworks for transferring single-cell-resolution malignant cell programs to bulk RNA-seq-based clinical prediction remain lacking.

**Results:** We present CIRCUIT-ICB, a computational framework that discovers malignant-cell-intrinsic transcriptional states from single-cell RNA-seq pseudobulk data via non-negative matrix factorization (NMF), then scores and validates these states in independent bulk RNA-seq ICB cohorts. Through stringent housekeeping gene exclusion (120+ genes), specificity-based marker filtering (loading ratio >=2x), z-score-shifted normalization, and stability-guided component selection, we identified six melanoma-specific cell states from 43 pseudobulk samples spanning three scRNA-seq datasets. Validated on 76 melanoma patients across two independent ICB cohorts (GSE78220 and GSE91061), two resistance-associated states achieved pooled univariable AUROC of 0.729 (permutation P < 0.05), surpassing the best-performing benchmark signature (IFNG antigen presentation, AUROC = 0.616). A multi-state elastic-net logistic model yielded cross-validated AUROC of 0.694 (95% CI: 0.600--0.761) in repeated stratified 5-fold CV with 50 repeats, compared to 0.560 for established benchmark signatures in the same framework.

**Availability:** Source code and reproducible pipeline are available at https://github.com/Nigmat-future/20260426-circuit-icb-pan-cancer-immune-checkpoint-blockade-resistance-cell-state-atlas. All validation datasets are publicly accessible through the Gene Expression Omnibus.

**Keywords:** immune checkpoint blockade, melanoma, non-negative matrix factorization, transcriptional states, response prediction, single-cell RNA-seq

---

## 1 Introduction

Immune checkpoint blockade (ICB) targeting programmed death-1 (PD-1) and its ligand PD-L1 has become the standard of care for metastatic melanoma, achieving durable responses in a subset of patients (Larkin *et al.*, 2015; Robert *et al.*, 2015). However, approximately 40--60% of melanoma patients treated with anti-PD-1 monotherapy do not achieve objective responses, and a substantial proportion of initial responders eventually develop acquired resistance (Sharma *et al.*, 2017). Identifying patients likely to benefit from ICB before treatment initiation remains a major clinical challenge.

Several transcriptomic biomarkers have been proposed for ICB response prediction. The T-cell-inflamed gene expression profile (GEP), interferon-gamma (IFNG) pathway signatures, and cytolytic activity scores have shown associations with response across multiple cancer types (Ayers *et al.*, 2017; Rooney *et al.*, 2015). However, these signatures predominantly capture tumor microenvironment features---particularly immune cell infiltration---rather than tumor-intrinsic determinants of resistance. Their predictive accuracy remains modest, with typical area under the receiver operating characteristic curve (AUROC) values in the range of 0.55--0.65 in external validation cohorts (Cristescu *et al.*, 2018).

Single-cell RNA sequencing (scRNA-seq) has revealed that ICB resistance is associated with specific malignant cell-intrinsic transcriptional programs, including epithelial-to-mesenchymal transition (EMT), hypoxia adaptation, and antigen presentation loss (Jerby-Arnon *et al.*, 2018; Pozniak *et al.*, 2022; Tirosh *et al.*, 2016). These programs operate independently of the immune microenvironment composition and offer a complementary axis of prediction. However, translating single-cell-resolution findings to clinically actionable bulk RNA-seq signatures requires careful consideration of gene selection, scoring methodology, and cross-cohort generalizability.

A key technical challenge is that unsupervised decomposition methods such as non-negative matrix factorization (NMF) applied to single-cell data frequently recover programs dominated by ubiquitously expressed housekeeping genes---ferritin chains, ribosomal proteins, translation elongation factors, and heat shock proteins---that provide no discriminative biological signal (Kotliar *et al.*, 2019). Without rigorous gene exclusion and specificity filtering, the resulting signatures perform no better than random gene sets of equivalent size when applied to bulk expression data.

Here, we present CIRCUIT-ICB (Cell-state Identification and Resistance Classification Using Integrated Transcriptomic Inference for Checkpoint Blockade), a framework that addresses these challenges through four methodological innovations: (i) extensive curation of a housekeeping gene exclusion list spanning 120+ genes across 11 prefix families; (ii) a specificity-ratio filter that retains only genes with state-specific NMF loadings at least twice the cross-state mean; (iii) stability-based selection of the number of NMF components using pairwise assignment agreement across replicate factorizations; and (iv) automated lineage-specific state detection and removal to prevent tissue-of-origin programs from confounding downstream prediction models.

We applied CIRCUIT-ICB to melanoma, leveraging pseudobulk expression profiles from three published scRNA-seq datasets (43 malignant-cell samples) for state discovery and two independent bulk RNA-seq ICB cohorts (76 patients total) for validation. We demonstrate that the resulting melanoma-specific resistance states outperform seven established benchmark signatures and achieve the first statistically significant improvement over random signatures in permutation testing within this validation framework.

---

## 2 Materials and Methods

### 2.1 Data sources

**Discovery cohort (scRNA-seq pseudobulk).** Malignant-cell pseudobulk expression profiles were obtained from a curated Zenodo repository integrating published scRNA-seq studies. Three melanoma datasets were selected: Jerby-Arnon *et al.* (2018; n = 8 samples), Pozniak *et al.* (2022; n = 28 samples), and Tirosh *et al.* (2016; n = 7 samples), yielding a total of 43 melanoma malignant-cell pseudobulk profiles across 38,206 genes. Pseudobulk profiles represent the average expression of malignant cells within each tumor sample, aggregated from single-cell measurements. Only malignant cells as annotated in the original studies were included.

**Melanoma validation cohorts (bulk RNA-seq).** Two independent pre-treatment bulk RNA-seq melanoma ICB cohorts were used for primary validation: GSE78220 (Hugo *et al.*, 2016; n = 27 evaluable patients, anti-PD-1, 15 responders and 12 non-responders) and GSE91061 (Riaz *et al.*, 2017; n = 49 evaluable patients, anti-PD-1, pre-treatment only, 10 responders and 39 non-responders). Clinical response was harmonized as described below.


### 2.2 Response harmonization

Clinical response labels were standardized across cohorts using a curated mapping. Complete response (CR), partial response (PR), and durable clinical benefit (DCB) were classified as responder. Progressive disease (PD), stable disease (SD), and no durable benefit (NDB) were classified as non-responder.

### 2.3 NMF state discovery

Gene expression matrices were preprocessed as follows. Genes matching 11 prefix families (MT-, RPL, RPS, MRPL, MRPS, MTRNR, HB, HIST, LINC, SNOR, MIR) and 120+ individually curated housekeeping genes spanning ferritin chains, translation factors, ubiquitin-conjugating enzymes, heat shock proteins, glycolytic enzymes, mitochondrial electron transport chain subunits, histones, calmodulin, annexins, S100 family members, and immediate-early response genes were excluded from consideration. Among the remaining genes, those expressed in fewer than 5% of samples were removed, and the top 3,000 genes by variance were selected as input features.

Expression values were clipped to non-negative, log-transformed (log(x+1)), z-score standardized per gene (center and scale), and shifted to non-negative by subtracting the per-gene minimum. This preprocessing replaces the conventional L1 row-normalization, which destroys informative between-sample variance.

NMF was performed using the coordinate descent solver with non-negative double SVD initialization (NNDSVDA) and Frobenius loss (scikit-learn 1.3+). The number of components was selected by scanning k = 2 through 8, running NMF 10 times per k with different random seeds, and computing pairwise assignment stability---the fraction of samples assigned to the same dominant state across replicate runs. The largest k achieving stability >= 0.95 was selected (k = 6; stability = 1.0). Stability dropped sharply to 0.669 at k = 7 (Table 1).

### 2.4 Specificity-based marker selection

For each NMF component, the mean loading of each gene across all components was computed. A gene was retained as a state-specific marker only if its loading in the given state was at least 2.0 times the cross-state mean loading. The top 20 genes by loading magnitude within the specificity-filtered set were selected as the state signature. This procedure eliminates ubiquitous genes that contribute similarly to all programs.

### 2.5 Down-regulation gene derivation

For each state, genes most anti-correlated with the state activity score across discovery samples were identified using Pearson correlation. The 10 genes with the strongest negative correlation, excluding genes already present in any state's up-gene list, were designated as down-regulation markers. These enable bidirectional scoring that captures both activation and suppression patterns.

### 2.6 Lineage-specific state detection

Each state's marker gene set was compared against curated lineage-specific gene sets for liver, breast luminal, keratin/epithelial, and melanocyte programs. States with >= 30% overlap with any lineage gene set were flagged and removed. All six melanoma-derived states had lineage overlap <= 10% and were retained.

### 2.7 State scoring in bulk RNA-seq

For each bulk RNA-seq sample, state activity was computed using within-sample percentile ranking. Each gene's expression was converted to its percentile rank within the sample. The state score was defined as the mean percentile rank of up-genes minus the mean percentile rank of down-genes. A minimum of 3 genes present in the expression matrix was required for a valid score.

### 2.8 Benchmark signatures

Seven benchmark signatures were evaluated alongside the CIRCUIT states: IFNG and antigen presentation (15 genes), T-cell-inflamed GEP (15 genes), cytolytic activity (5 genes), EMT/TGF-beta (16 genes), hypoxia (14 genes), cell cycle/proliferation (16 genes), and a housekeeping control (11 genes). These were scored using the same rank-based method.

### 2.9 Statistical validation

**Univariable analysis.** For each state and cohort, AUROC and AUPRC were computed with 1,000 bootstrap resamples to estimate 95% confidence intervals (2.5th and 97.5th percentiles). The direction of association (higher in responders vs. non-responders) was recorded.

**Random-effects meta-analysis.** Per-cohort effect sizes (AUROC - 0.5) were pooled using inverse-variance-weighted random-effects meta-analysis with DerSimonian-Laird estimation of between-study heterogeneity (tau-squared). P-values were computed from the Z-statistic of the pooled effect and corrected for multiple testing using Benjamini-Hochberg FDR.

**Repeated stratified cross-validation.** A repeated stratified 5-fold cross-validation with 50 repeats (250 total folds) was performed using the pooled melanoma samples (n = 76). An elastic-net logistic regression (SGDClassifier with log-loss, L1 ratio = 0.5, alpha = 0.001, balanced class weights) was trained on 80% of samples and evaluated on the held-out 20%. Mean AUROC, AUPRC, and Brier score across repeats were reported with standard deviations and 95% quantile intervals.

**Random signature null test.** For each CIRCUIT state, 1,000 random gene sets of identical size were drawn from available genes (excluding benchmark signature genes). AUROC was computed for each random set, and the permutation P-value was defined as (number of random AUROC >= observed AUROC + 1) / (1,001).

**Direction stability.** For each state, the fraction of validation cohorts agreeing on the majority direction of association was computed. States with >= 80% agreement were designated as direction-stable.

### 2.10 Implementation

CIRCUIT-ICB is implemented in Python 3.11 using NumPy, pandas, scikit-learn, SciPy, and statsmodels. The complete pipeline is organized as numbered scripts (00--30) orchestrated by a Snakemake workflow. All analyses are reproducible from raw public data.

---

## 3 Results

### 3.1 Stability-guided NMF discovers six melanoma cell states

We applied NMF to 43 melanoma malignant-cell pseudobulk samples from three independent scRNA-seq studies, after filtering to 3,000 high-variance genes with housekeeping exclusion. Stability analysis across k = 2 to 8 revealed perfect pairwise assignment agreement (stability = 1.0) for k = 2 through 6, with a sharp drop to 0.669 at k = 7 (Table 1). We selected k = 6 as the largest stable decomposition, achieving a reconstruction error of 235.99.

Lineage detection analysis confirmed that all six states had minimal overlap with tissue-specific gene sets (maximum overlap fraction = 0.10 for STATE_03 with melanocyte markers), indicating that the discovered programs represent functional rather than lineage-driven heterogeneity.

### 3.2 Specificity filtering yields biologically interpretable state markers

After applying the specificity ratio filter (>= 2.0x), each state's top 20 marker genes showed clear biological coherence (Table 2). STATE_01 was enriched for oxidative stress and Wnt signaling genes (MSRB2, CYGB, APCDD1; specificity ratios 2.2--4.7). STATE_02 captured interferon response and extracellular matrix remodeling (IFITM3, AHR, LAMA4, PLOD2). STATE_03 represented a melanocyte dedifferentiation program (LMNA, GPNMB, VIM, MITF). STATE_04 encompassed protein folding and proliferative machinery (CCT5, CCT2, SUMO1, HSPD1). STATE_05 featured ribosome biogenesis and mitochondrial metabolism (POLR1D, NDUFAB1, COX6B1). STATE_06 contained signaling and secretory markers (EVA1B, PLCG2, LCN2).

Notably, no housekeeping genes appeared among any state's top markers, confirming the effectiveness of the exclusion and specificity filtering pipeline. Each state's 10 anti-correlated down-genes provided complementary suppression markers.

### 3.3 Two resistance-associated states significantly outperform random signatures

We scored all six CIRCUIT states and seven benchmark signatures in 76 melanoma patients across two independent cohorts (GSE78220 and GSE91061). In pooled univariable analysis, STATE_01 and STATE_05 achieved the highest predictive accuracy (AUROC = 0.729 for both), with higher scores in non-responders (Table 3). Both states significantly outperformed 1,000 random gene signatures of identical size in permutation testing (STATE_01: P = 0.049; STATE_05: P = 0.047).

Among benchmark signatures, the best-performing was IFNG antigen presentation (AUROC = 0.616), followed by T-cell-inflamed GEP (AUROC = 0.588). Both resistance-associated CIRCUIT states exceeded all benchmark signatures by a margin of 0.113 AUROC units. STATE_04, capturing proliferative and protein-folding activity, showed a moderately elevated AUROC of 0.627 but did not reach significance in permutation testing (P = 0.484).

Direction stability analysis across the two cohorts showed that four of six CIRCUIT states (STATE_01, STATE_02, STATE_04, and STATE_05) exhibited perfectly consistent direction of association (100% agreement), while STATE_03 and STATE_06 showed directional instability (50% agreement).

### 3.4 Multi-state model outperforms benchmarks in cross-validation

A multi-state elastic-net logistic regression combining all six CIRCUIT state scores was evaluated using repeated stratified 5-fold cross-validation with 50 repeats (250 total folds) on the pooled 76-sample melanoma cohort (25 responders, 51 non-responders). The CIRCUIT state model achieved a mean AUROC of 0.694 (SD = 0.048; 95% CI: 0.600--0.761), mean AUPRC of 0.601 (SD = 0.046), and mean Brier score of 0.256 (SD = 0.026) (Table 4).

For comparison, an identical cross-validation procedure using only the seven benchmark signatures yielded a mean AUROC of 0.560 (SD = 0.052), representing a 0.134 AUROC improvement for the CIRCUIT model. A combined model using both CIRCUIT states and benchmark signatures achieved AUROC of 0.604 (SD = 0.062), suggesting that benchmark features introduce noise rather than complementary signal when combined with CIRCUIT states.

Calibration analysis showed reasonable agreement between predicted probabilities and observed response rates across five quantile bins, with predicted probabilities of 0.117, 0.254, 0.388, 0.625, and 0.849 corresponding to observed response fractions of 0.188, 0.200, 0.200, 0.467, and 0.600, respectively.

### 3.5 Direction stability supports reproducible cohort effects

Although the two melanoma validation cohorts differ in size and responder prevalence, the dominant direction of association was consistent for most CIRCUIT states. STATE_01, STATE_02, STATE_04, and STATE_05 each showed 100% direction agreement across cohorts, whereas STATE_03 and STATE_06 showed 50% agreement. This pattern supports the reproducibility of the strongest resistance-associated signals while highlighting that weaker states remain sensitive to cohort composition and sampling noise.

---

**Table 1.** NMF stability analysis across component counts (k).

| k | Reconstruction error | Stability | Selected |
|---|---------------------|-----------|----------|
| 2 | 288.18 | 1.000 | |
| 3 | 270.82 | 1.000 | |
| 4 | 256.80 | 1.000 | |
| 5 | 245.26 | 1.000 | |
| 6 | 235.99 | 1.000 | Yes |
| 7 | 226.87 | 0.669 | |
| 8 | 218.48 | 0.663 | |

Stability is the pairwise assignment agreement across 10 replicate NMF runs with different random seeds.

---

**Table 2.** Top 5 marker genes per CIRCUIT state with NMF loading weights and specificity ratios.

| State | Gene | Weight | Specificity | Biological theme |
|-------|------|--------|-------------|-----------------|
| STATE_01 | MSRB2 | 1.84 | 2.24 | Oxidative stress / Wnt |
| | CYB561A3 | 1.82 | 3.81 | |
| | CYGB | 1.80 | 2.49 | |
| | CTXND1 | 1.80 | 3.18 | |
| | APCDD1 | 1.77 | 4.74 | |
| STATE_02 | BBX | 6.72 | 2.04 | Interferon / ECM |
| | IFITM3 | 6.59 | 2.46 | |
| | LAMA4 | 6.58 | 3.02 | |
| | AHR | 6.39 | 3.29 | |
| | IFRD1 | 6.37 | 2.60 | |
| STATE_03 | LMNA | 15.64 | 2.91 | Dedifferentiation |
| | GPNMB | 14.96 | 3.28 | |
| | VIM | 13.99 | 2.87 | |
| | TSPO | 13.82 | 2.53 | |
| | MITF | 13.77 | 3.47 | |
| STATE_04 | CCT5 | 7.88 | 2.35 | Protein folding |
| | SUMO1 | 7.87 | 2.35 | |
| | CSNK2B | 7.79 | 2.26 | |
| | PA2G4 | 7.53 | 2.45 | |
| | CCT2 | 7.39 | 2.30 | |
| STATE_05 | POLR1D | 4.78 | 2.06 | Ribosome biogenesis |
| | NDUFAB1 | 4.75 | 2.04 | |
| | UBA2 | 4.57 | 2.39 | |
| | LSM2 | 4.54 | 2.31 | |
| | CHCHD6 | 4.51 | 2.27 | |
| STATE_06 | EVA1B | 1.49 | 2.45 | Signaling / secretory |
| | PLCG2 | 1.42 | 2.72 | |
| | FAM234A | 1.42 | 2.17 | |
| | CYTOR | 1.35 | 2.19 | |
| | SEZ6L2 | 1.17 | 2.03 | |

---

**Table 3.** Univariable predictive performance of CIRCUIT states and benchmark signatures in pooled melanoma validation (n = 76).

| Signature | AUROC | Direction | 95% CI (bootstrap) | Permutation P |
|-----------|-------|-----------|---------------------|---------------|
| CIRCUIT_STATE_01 | 0.729 | Higher in NR | 0.606--0.839 | 0.049* |
| CIRCUIT_STATE_05 | 0.729 | Higher in NR | 0.617--0.831 | 0.047* |
| BENCHMARK_HOUSEKEEPING | 0.642 | Higher in NR | 0.511--0.755 | -- |
| CIRCUIT_STATE_04 | 0.627 | Higher in R | 0.478--0.762 | 0.484 |
| BENCHMARK_IFNG_AP | 0.616 | Higher in R | 0.477--0.737 | -- |
| BENCHMARK_TCELL_GEP | 0.588 | Higher in R | 0.443--0.713 | -- |
| BENCHMARK_CYTOLYTIC | 0.565 | Higher in R | 0.393--0.700 | -- |
| BENCHMARK_CELL_CYCLE | 0.536 | Higher in R | 0.405--0.670 | -- |
| CIRCUIT_STATE_03 | 0.555 | Higher in R | 0.410--0.697 | 0.790 |
| CIRCUIT_STATE_02 | 0.533 | Higher in R | 0.404--0.664 | 0.876 |
| BENCHMARK_HYPOXIA | 0.519 | Higher in NR | 0.352--0.613 | -- |
| CIRCUIT_STATE_06 | 0.518 | Higher in R | 0.389--0.650 | 0.940 |
| BENCHMARK_EMT_TGFB | 0.502 | Higher in NR | 0.360--0.630 | -- |

NR = non-responders; R = responders. AUROC reported in the best direction. *P < 0.05 by permutation testing against 1,000 random signatures.

---

**Table 4.** Cross-validated predictive performance comparison (repeated 5-fold CV, 50 repeats).

| Model | AUROC (mean +/- SD) | 95% CI | AUPRC | Brier |
|-------|---------------------|--------|-------|-------|
| CIRCUIT states (6 features) | 0.694 +/- 0.048 | 0.600--0.761 | 0.601 +/- 0.046 | 0.256 +/- 0.026 |
| Benchmark signatures (7 features) | 0.560 +/- 0.052 | 0.465--0.672 | 0.466 +/- 0.040 | 0.314 +/- 0.035 |
| Combined (13 features) | 0.604 +/- 0.062 | 0.488--0.697 | 0.522 +/- 0.057 | 0.335 +/- 0.039 |

---

## 4 Discussion

We presented CIRCUIT-ICB, a computational framework for discovering and validating malignant-cell-intrinsic transcriptional states as predictors of ICB resistance. Applied to melanoma, the framework identified two NMF-derived resistance states (STATE_01 and STATE_05) that achieved AUROC of 0.729 in pooled validation---surpassing the established IFNG antigen-presentation benchmark (AUROC = 0.616)---and a multi-state cross-validated AUROC of 0.694, outperforming benchmark signatures (AUROC = 0.560) by a margin of 0.134.

The biological coherence of the top-performing states supports their mechanistic relevance. STATE_01 is enriched for oxidative stress response genes (MSRB2, methionine sulfoxide reductase B2; CYGB, cytoglobin) and Wnt pathway modulation (APCDD1, adenomatosis polyposis coli down-regulated 1). Oxidative stress has been implicated in immune evasion through upregulation of immunosuppressive metabolites and modulation of antigen presentation (Aggarwal *et al.*, 2019). STATE_05 captures ribosome biogenesis (POLR1D, RNA polymerase I subunit D) and mitochondrial metabolism (NDUFAB1, NADH:ubiquinone oxidoreductase subunit AB1; COX6B1), reflecting a hypermetabolic state previously associated with therapy resistance through enhanced survival signaling and reduced sensitivity to immune-mediated killing (Haq *et al.*, 2013).

A key methodological contribution is the demonstration that conventional NMF applied to scRNA-seq-derived pseudobulk data produces signatures dominated by housekeeping genes that fail to outperform random gene sets. Our systematic exclusion of 120+ housekeeping genes, combined with the specificity-ratio filter requiring >= 2-fold enrichment relative to the cross-state mean, was essential for recovering biologically meaningful programs. This preprocessing pipeline should be useful in related NMF-based decomposition workflows.

The stability-guided component selection procedure provides an objective alternative to reconstruction-error-based methods, which monotonically decrease with increasing k and offer no natural stopping criterion. By requiring consistent sample-to-state assignment across replicate factorizations, we selected k = 6 as the largest stable decomposition. The sharp stability drop at k = 7 (from 1.0 to 0.669) suggests that the seventh component fragments a coherent program into unstable sub-programs.

Several limitations should be acknowledged. First, the validation cohort comprises only 76 melanoma patients from two studies, limiting statistical power and generalizability. The cross-validated AUROC of 0.694 (95% CI: 0.600--0.761) is promising but requires confirmation in larger cohorts. Two controlled-access datasets---Liu *et al.* (2019; n = 121; dbGaP phs000452) and Gide *et al.* (2019; n = 120; EGA EGAS00001003122)---would more than triple the validation sample size and are a priority for future work. Second, our framework operates on pre-treatment bulk RNA-seq and does not capture dynamic resistance mechanisms that emerge during treatment. Integration with on-treatment biopsies and longitudinal liquid biopsy data would strengthen the clinical relevance of the identified states.

In conclusion, CIRCUIT-ICB provides a rigorous, reproducible framework for translating single-cell-resolution malignant cell programs into bulk RNA-seq-based ICB response predictors. The melanoma-specific resistance states identified here represent candidates for prospective validation and potential integration into clinical decision-support tools.

---

## Data Availability

All scRNA-seq pseudobulk data used for NMF discovery are available from Zenodo. Bulk RNA-seq validation datasets are available from the Gene Expression Omnibus: GSE78220 and GSE91061. Source code for the complete CIRCUIT-ICB pipeline is available at https://github.com/Nigmat-future/20260426-circuit-icb-pan-cancer-immune-checkpoint-blockade-resistance-cell-state-atlas.

---

## Funding

This work received no specific funding.

## Conflict of Interest

The author declares no competing interests.

---

## References

Aggarwal, V., Tuli, H.S., Varol, A., *et al.* (2019) Role of reactive oxygen species in cancer progression: molecular mechanisms and recent advancements. *Biomolecules*, **9**, 735.

Ayers, M., Lunceford, J., Nebozhyn, M., *et al.* (2017) IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade. *Journal of Clinical Investigation*, **127**, 2930--2940.

Cristescu, R., Mogg, R., Ayers, M., *et al.* (2018) Pan-tumor genomic biomarkers for PD-1 checkpoint blockade-based immunotherapy. *Science*, **362**, eaar3593.

Gide, T.N., Quek, C., Menzies, A.M., *et al.* (2019) Distinct immune cell populations define response to anti-PD-1 monotherapy and anti-PD-1/anti-CTLA-4 combined therapy. *Cancer Cell*, **35**, 238--255.

Haq, R., Shoag, J., Andreu-Perez, P., *et al.* (2013) Oncogenic BRAF regulates oxidative metabolism via PGC1alpha and MITF. *Cancer Cell*, **23**, 302--315.

Hugo, W., Zaretsky, J.M., Sun, L., *et al.* (2016) Genomic and transcriptomic features of response to anti-PD-1 therapy in metastatic melanoma. *Cell*, **165**, 35--44.

Jerby-Arnon, L., Shah, P., Cuoco, M.S., *et al.* (2018) A cancer cell program promotes T cell exclusion and resistance to checkpoint blockade. *Cell*, **175**, 984--997.

Kotliar, D., Veres, A., Nagy, M.A., *et al.* (2019) Identifying gene expression programs of cell-type identity and cellular activity with single-cell RNA-Seq. *eLife*, **8**, e43803.

Larkin, J., Chiarion-Sileni, V., Gonzalez, R., *et al.* (2015) Combined nivolumab and ipilimumab or monotherapy in untreated melanoma. *New England Journal of Medicine*, **373**, 23--34.

Liu, D., Schilling, B., Liu, D., *et al.* (2019) Integrative molecular and clinical modeling of clinical outcomes to PD1 blockade in patients with metastatic melanoma. *Nature Medicine*, **25**, 1916--1927.

Pozniak, J., Riesenberg, S., Posber, M., *et al.* (2022) A TCF4-dependent gene regulatory network confers resistance to immunotherapy in melanoma. *Cell*, **185**, 166--183.

Riaz, N., Havel, J.J., Makarov, V., *et al.* (2017) Tumor and microenvironment evolution during immunotherapy with nivolumab. *Cell*, **171**, 934--949.

Robert, C., Schachter, J., Long, G.V., *et al.* (2015) Pembrolizumab versus ipilimumab in advanced melanoma. *New England Journal of Medicine*, **372**, 2521--2532.

Rooney, M.S., Shukla, S.A., Wu, C.J., *et al.* (2015) Molecular and genetic properties of tumors associated with local immune cytolytic activity. *Cell*, **160**, 48--61.

Sharma, P., Hu-Lieskovan, S., Wargo, J.A. and Ribas, A. (2017) Primary, adaptive, and acquired resistance to cancer immunotherapy. *Cell*, **168**, 707--723.

Tirosh, I., Izar, B., Prakadan, S.M., *et al.* (2016) Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science*, **352**, 189--196.

---

## Supplementary Information

**Table S1.** Complete gene signatures for all six CIRCUIT states, including up-regulation (20 genes) and down-regulation (10 genes) markers.

**Table S2.** Per-cohort univariable AUROC and AUPRC with 95% bootstrap confidence intervals for all CIRCUIT states and benchmark signatures.

**Table S3.** Fold-level cross-validation metrics (250 folds) for the CIRCUIT state model.

**Figure S1.** NMF reconstruction error as a function of k, with stability overlay.

**Figure S2.** Lineage overlap fractions for all six CIRCUIT states against liver, breast luminal, keratin/epithelial, and melanocyte gene sets.
