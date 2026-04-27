"""Generate a complete Word document (.docx) for the CIRCUIT-ICB manuscript."""
from __future__ import annotations

import sys
from pathlib import Path

from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, Cm, RGBColor
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "manuscript" / "figures"
MANSRC = ROOT / "manuscript" / "manuscript.md"
RES = ROOT / "results" / "tables" / "melanoma_only"
VAL = RES / "validation"

doc = Document()

# ── Page setup ──
for section in doc.sections:
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)
    section.left_margin = Cm(2.54)
    section.right_margin = Cm(2.54)

# ── Style helpers ──
style_normal = doc.styles["Normal"]
style_normal.font.name = "Times New Roman"
style_normal.font.size = Pt(10)
style_normal.paragraph_format.space_after = Pt(4)
style_normal.paragraph_format.line_spacing = 1.15
# Set East Asian font
rpr = style_normal.element.get_or_add_rPr()
rFonts = rpr.find(qn("w:rFonts"))
if rFonts is None:
    from docx.oxml import OxmlElement
    rFonts = OxmlElement("w:rFonts")
    rpr.insert(0, rFonts)
rFonts.set(qn("w:eastAsia"), "Times New Roman")


def add_heading(text: str, level: int = 1):
    h = doc.add_heading(text, level=level)
    for run in h.runs:
        run.font.name = "Times New Roman"
        run.font.color.rgb = RGBColor(0, 0, 0)


def add_para(text: str, bold: bool = False, italic: bool = False, size: int = 10,
             alignment=WD_ALIGN_PARAGRAPH.JUSTIFY):
    p = doc.add_paragraph()
    p.alignment = alignment
    run = p.add_run(text)
    run.font.name = "Times New Roman"
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    return p


def add_rich_para(segments: list[tuple[str, bool, bool]], alignment=WD_ALIGN_PARAGRAPH.JUSTIFY):
    """Add paragraph with mixed formatting. segments = [(text, bold, italic), ...]"""
    p = doc.add_paragraph()
    p.alignment = alignment
    for text, bold, italic in segments:
        run = p.add_run(text)
        run.font.name = "Times New Roman"
        run.font.size = Pt(10)
        run.bold = bold
        run.italic = italic
    return p


def add_figure(filename: str, caption: str, width: float = 5.8):
    path = FIG / filename
    if path.exists():
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run()
        run.add_picture(str(path), width=Inches(width))
    cap = doc.add_paragraph()
    cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = cap.add_run(caption)
    run.font.name = "Times New Roman"
    run.font.size = Pt(9)
    run.italic = True


def format_auc_ci(low: float, high: float) -> str:
    return f"{low:.3f}–{high:.3f}"


def format_mean_sd(mean: float, sd: float) -> str:
    return f"{mean:.3f} ± {sd:.3f}"


def signature_display_name(name: str) -> str:
    mapping = {
        "BENCHMARK_IFNG_ANTIGEN_PRESENTATION_score": "BENCH_IFNG_AP",
        "BENCHMARK_TCELL_INFLAMED_GEP_score": "BENCH_TCELL_GEP",
        "BENCHMARK_CYTOLYTIC_SCORE_score": "BENCH_CYTOLYTIC",
        "BENCHMARK_EMT_TGFB_score": "BENCH_EMT_TGFB",
        "BENCHMARK_HYPOXIA_score": "BENCH_HYPOXIA",
        "BENCHMARK_CELL_CYCLE_score": "BENCH_CELL_CYCLE",
        "BENCHMARK_HOUSEKEEPING_CONTROL_score": "BENCH_HOUSEKEEPING",
    }
    if name in mapping:
        return mapping[name]
    return name.replace("_score", "")


# ════════════════════════════════════════════════════════════════
# TITLE
# ════════════════════════════════════════════════════════════════
add_para(
    "CIRCUIT-ICB: NMF-based discovery of malignant cell states predicting "
    "immune checkpoint blockade resistance in melanoma",
    bold=True, size=14, alignment=WD_ALIGN_PARAGRAPH.CENTER,
)

# Authors
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run("Nigmat Rahim")
run.font.name = "Times New Roman"
run.font.size = Pt(11)

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run("Institute of Medical Technology, Peking University, Beijing, China")
run.font.name = "Times New Roman"
run.font.size = Pt(9)
run.italic = True

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run("Corresponding author: nigmatrahim@stu.pku.edu.cn")
run.font.name = "Times New Roman"
run.font.size = Pt(9)
run.italic = True

doc.add_paragraph()  # spacer

# ════════════════════════════════════════════════════════════════
# ABSTRACT
# ════════════════════════════════════════════════════════════════
add_heading("Abstract", level=1)

add_rich_para([
    ("Motivation: ", True, False),
    ("Immune checkpoint blockade (ICB) has transformed the treatment of metastatic melanoma, yet 40–60% of patients fail to respond. Existing transcriptomic biomarkers such as interferon-gamma (IFNG) signatures and tumor mutational burden provide modest predictive accuracy, and systematic frameworks for transferring single-cell-resolution malignant cell programs to bulk RNA-seq-based clinical prediction remain lacking.", False, False),
])

add_rich_para([
    ("Results: ", True, False),
    ("We present CIRCUIT-ICB, a computational framework that discovers malignant-cell-intrinsic transcriptional states from single-cell RNA-seq pseudobulk data via non-negative matrix factorization (NMF), then scores and validates these states in independent bulk RNA-seq ICB cohorts. Through stringent housekeeping gene exclusion (120+ genes), specificity-based marker filtering (loading ratio ≥2×), z-score-shifted normalization, and stability-guided component selection, we identified six melanoma-specific cell states from 43 pseudobulk samples spanning three scRNA-seq datasets. Validated on 76 melanoma patients across two independent ICB cohorts (GSE78220 and GSE91061), two resistance-associated states achieved pooled univariable AUROC of 0.729 (permutation P < 0.05), surpassing the best-performing benchmark signature (IFNG antigen presentation, AUROC = 0.616). A multi-state elastic-net logistic model yielded cross-validated AUROC of 0.694 (95% CI: 0.600–0.761) in repeated stratified 5-fold CV with 50 repeats, compared to 0.560 for established benchmark signatures.", False, False),
])

add_rich_para([
    ("Availability: ", True, False),
    ("Source code and reproducible pipeline are publicly available at https://github.com/Nigmat-future/20260426-circuit-icb-pan-cancer-immune-checkpoint-blockade-resistance-cell-state-atlas. All validation datasets are accessible through the Gene Expression Omnibus.", False, False),
])

add_rich_para([
    ("Contact: ", True, False),
    ("nigmatrahim@stu.pku.edu.cn", False, True),
])

# ════════════════════════════════════════════════════════════════
# 1 INTRODUCTION
# ════════════════════════════════════════════════════════════════
add_heading("1 Introduction", level=1)

add_para(
    "Immune checkpoint blockade (ICB) targeting programmed death-1 (PD-1) and its ligand PD-L1 "
    "has become the standard of care for metastatic melanoma, achieving durable responses in a "
    "subset of patients (Larkin et al., 2015; Robert et al., 2015). However, approximately "
    "40–60% of melanoma patients treated with anti-PD-1 monotherapy do not achieve objective "
    "responses, and a substantial proportion of initial responders eventually develop acquired "
    "resistance (Sharma et al., 2017). Identifying patients likely to benefit from ICB before "
    "treatment initiation remains a major clinical challenge."
)

add_para(
    "Several transcriptomic biomarkers have been proposed for ICB response prediction. The "
    "T-cell-inflamed gene expression profile (GEP), interferon-gamma (IFNG) pathway signatures, "
    "and cytolytic activity scores have shown associations with response across multiple cancer "
    "types (Ayers et al., 2017; Rooney et al., 2015). However, these signatures predominantly "
    "capture tumor microenvironment features—particularly immune cell infiltration—rather than "
    "tumor-intrinsic determinants of resistance. Their predictive accuracy remains modest, with "
    "typical area under the receiver operating characteristic curve (AUROC) values in the range "
    "of 0.55–0.65 in external validation cohorts (Cristescu et al., 2018)."
)

add_para(
    "Single-cell RNA sequencing (scRNA-seq) has revealed that ICB resistance is associated "
    "with specific malignant cell-intrinsic transcriptional programs, including "
    "epithelial-to-mesenchymal transition (EMT), hypoxia adaptation, and antigen presentation "
    "loss (Jerby-Arnon et al., 2018; Pozniak et al., 2022; Tirosh et al., 2016). These "
    "programs operate independently of the immune microenvironment composition and offer a "
    "complementary axis of prediction. However, translating single-cell-resolution findings to "
    "clinically actionable bulk RNA-seq signatures requires careful consideration of gene "
    "selection, scoring methodology, and cross-cohort generalizability."
)

add_para(
    "A key technical challenge is that unsupervised decomposition methods such as non-negative "
    "matrix factorization (NMF) applied to single-cell data frequently recover programs "
    "dominated by ubiquitously expressed housekeeping genes—ferritin chains, ribosomal "
    "proteins, translation elongation factors, and heat shock proteins—that provide no "
    "discriminative biological signal (Kotliar et al., 2019). Without rigorous gene exclusion "
    "and specificity filtering, the resulting signatures perform no better than random gene "
    "sets of equivalent size when applied to bulk expression data."
)

add_para(
    "Here, we present CIRCUIT-ICB (Cell-state Identification and Resistance Classification "
    "Using Integrated Transcriptomic Inference for Checkpoint Blockade), a framework that "
    "addresses these challenges through four methodological innovations: (i) extensive curation "
    "of a housekeeping gene exclusion list spanning 120+ genes across 11 prefix families; "
    "(ii) a specificity-ratio filter that retains only genes with state-specific NMF loadings "
    "at least twice the cross-state mean; (iii) stability-based selection of the number of NMF "
    "components using pairwise assignment agreement across replicate factorizations; and (iv) "
    "automated lineage-specific state detection and removal to prevent tissue-of-origin "
    "programs from confounding prediction models."
)

add_para(
    "We applied CIRCUIT-ICB to melanoma, leveraging pseudobulk expression profiles from three "
    "published scRNA-seq datasets (43 malignant-cell samples) for state discovery and two "
    "independent bulk RNA-seq ICB cohorts (76 patients total) for validation. We demonstrate "
    "that the resulting melanoma-specific resistance states outperform seven established "
    "benchmark signatures and achieve the first statistically significant improvement over "
    "random signatures in permutation testing within this validation framework."
)

# ════════════════════════════════════════════════════════════════
# 2 MATERIALS AND METHODS
# ════════════════════════════════════════════════════════════════
add_heading("2 Materials and Methods", level=1)

# 2.1
add_heading("2.1 Data sources", level=2)
add_rich_para([
    ("Discovery cohort (scRNA-seq pseudobulk). ", True, False),
    ("Malignant-cell pseudobulk expression profiles were obtained from a curated Zenodo "
     "repository integrating published scRNA-seq studies. Three melanoma datasets were selected: "
     "Jerby-Arnon et al. (2018; n = 8 samples), Pozniak et al. (2022; n = 28 samples), and "
     "Tirosh et al. (2016; n = 7 samples), yielding a total of 43 melanoma malignant-cell "
     "pseudobulk profiles across 38,206 genes. Pseudobulk profiles represent the average "
     "expression of malignant cells within each tumor sample, aggregated from single-cell "
     "measurements. Only malignant cells as annotated in the original studies were included.", False, False),
])

add_rich_para([
    ("Melanoma validation cohorts (bulk RNA-seq). ", True, False),
    ("Two independent pre-treatment bulk RNA-seq melanoma ICB cohorts were used for primary "
     "validation: GSE78220 (Hugo et al., 2016; n = 27 evaluable patients, anti-PD-1, 15 "
     "responders and 12 non-responders) and GSE91061 (Riaz et al., 2017; n = 49 evaluable "
     "patients, anti-PD-1, pre-treatment only, 10 responders and 39 non-responders).", False, False),
])

# 2.2
add_heading("2.2 Response harmonization", level=2)
add_para(
    "Clinical response labels were standardized across cohorts using a curated mapping. "
    "Complete response (CR), partial response (PR), and durable clinical benefit (DCB) were "
    "classified as responder. Progressive disease (PD), stable disease (SD), and no durable "
    "benefit (NDB) were classified as non-responder."
)

# 2.3
add_heading("2.3 NMF state discovery", level=2)
add_para(
    "Gene expression matrices were preprocessed as follows. Genes matching 11 prefix families "
    "(MT-, RPL, RPS, MRPL, MRPS, MTRNR, HB, HIST, LINC, SNOR, MIR) and 120+ individually "
    "curated housekeeping genes spanning ferritin chains, translation factors, ubiquitin-"
    "conjugating enzymes, heat shock proteins, glycolytic enzymes, mitochondrial electron "
    "transport chain subunits, histones, calmodulin, annexins, S100 family members, and "
    "immediate-early response genes were excluded. Among the remaining genes, those expressed "
    "in fewer than 5% of samples were removed, and the top 3,000 genes by variance were "
    "selected as input features."
)
add_para(
    "Expression values were clipped to non-negative, log-transformed (log(x+1)), z-score "
    "standardized per gene (center and scale), and shifted to non-negative by subtracting the "
    "per-gene minimum. NMF was performed using the coordinate descent solver with NNDSVDA "
    "initialization and Frobenius loss (scikit-learn). The number of components was selected "
    "by scanning k = 2 through 8, running NMF 10 times per k with different random seeds, and "
    "computing pairwise assignment stability. The largest k achieving stability ≥ 0.95 was "
    "selected (k = 6; stability = 1.0). Stability dropped sharply to 0.669 at k = 7."
)

# 2.4
add_heading("2.4 Specificity-based marker selection", level=2)
add_para(
    "For each NMF component, the mean loading of each gene across all components was computed. "
    "A gene was retained as a state-specific marker only if its loading in the given state was "
    "at least 2.0 times the cross-state mean loading. The top 20 genes by loading magnitude "
    "within the specificity-filtered set were selected as the state signature."
)

# 2.5
add_heading("2.5 Down-regulation gene derivation", level=2)
add_para(
    "For each state, genes most anti-correlated with the state activity score across discovery "
    "samples were identified using Pearson correlation. The 10 genes with the strongest negative "
    "correlation, excluding genes already present in any state's up-gene list, were designated "
    "as down-regulation markers."
)

# 2.6
add_heading("2.6 Lineage-specific state detection", level=2)
add_para(
    "Each state's marker gene set was compared against curated lineage-specific gene sets for "
    "liver, breast luminal, keratin/epithelial, and melanocyte programs. States with ≥30% "
    "overlap with any lineage gene set were flagged and removed. All six melanoma-derived "
    "states had lineage overlap ≤10% and were retained."
)

# 2.7
add_heading("2.7 State scoring in bulk RNA-seq", level=2)
add_para(
    "For each bulk RNA-seq sample, state activity was computed using within-sample percentile "
    "ranking. Each gene's expression was converted to its percentile rank within the sample. "
    "The state score was defined as the mean percentile rank of up-genes minus the mean "
    "percentile rank of down-genes. A minimum of 3 genes present in the expression matrix was "
    "required for a valid score."
)

# 2.8
add_heading("2.8 Benchmark signatures", level=2)
add_para(
    "Seven benchmark signatures were evaluated alongside the CIRCUIT states: IFNG and antigen "
    "presentation (15 genes), T-cell-inflamed GEP (15 genes), cytolytic activity (5 genes), "
    "EMT/TGF-beta (16 genes), hypoxia (14 genes), cell cycle/proliferation (16 genes), and "
    "a housekeeping control (11 genes). These were scored using the same rank-based method."
)

# 2.9
add_heading("2.9 Statistical validation", level=2)
add_rich_para([
    ("Univariable analysis. ", True, False),
    ("For each state and cohort, AUROC and AUPRC were computed with 1,000 bootstrap resamples "
     "to estimate 95% confidence intervals (2.5th and 97.5th percentiles).", False, False),
])
add_rich_para([
    ("Repeated stratified cross-validation. ", True, False),
    ("A repeated stratified 5-fold cross-validation with 50 repeats (250 total folds) was "
     "performed on the pooled melanoma samples (n = 76). An elastic-net logistic regression "
     "(SGDClassifier with log-loss, L1 ratio = 0.5, alpha = 0.001, balanced class weights) "
     "was trained on 80% of samples and evaluated on the held-out 20%.", False, False),
])
add_rich_para([
    ("Random signature null test. ", True, False),
    ("For each CIRCUIT state, 1,000 random gene sets of identical size were drawn from "
     "available genes (excluding benchmark signature genes). AUROC was computed for each "
     "random set, and the permutation P-value was defined as (count ≥ observed + 1) / "
     "(1,001).", False, False),
])
add_rich_para([
    ("Direction stability. ", True, False),
    ("For each state, the fraction of validation cohorts agreeing on the majority direction "
     "of association was computed. States with ≥80% agreement were designated as "
     "direction-stable.", False, False),
])

# 2.10
add_heading("2.10 Implementation", level=2)
add_para(
    "CIRCUIT-ICB is implemented in Python 3.11 using NumPy, pandas, scikit-learn, SciPy, and "
    "statsmodels. The complete pipeline is organized as numbered scripts (00–30) orchestrated "
    "by a Snakemake workflow. All analyses are reproducible from raw public data."
)

# ════════════════════════════════════════════════════════════════
# 3 RESULTS
# ════════════════════════════════════════════════════════════════
add_heading("3 Results", level=1)

# 3.1
add_heading("3.1 Stability-guided NMF discovers six melanoma cell states", level=2)
add_para(
    "We applied NMF to 43 melanoma malignant-cell pseudobulk samples from three independent "
    "scRNA-seq studies, after filtering to 3,000 high-variance genes with housekeeping "
    "exclusion. Stability analysis across k = 2 to 8 revealed perfect pairwise assignment "
    "agreement (stability = 1.0) for k = 2 through 6, with a sharp drop to 0.669 at k = 7 "
    "(Figure 1). We selected k = 6 as the largest stable decomposition, achieving a "
    "reconstruction error of 235.99. Lineage detection analysis confirmed that all six states "
    "had minimal overlap with tissue-specific gene sets (maximum overlap fraction = 0.10), "
    "indicating functional rather than lineage-driven heterogeneity."
)

add_figure("fig1_k_scan.png",
           "Figure 1. NMF component selection by stability analysis. Bars show reconstruction "
           "error; blue line shows pairwise assignment stability across 10 replicate factorizations. "
           "Dashed line indicates 95% stability threshold. k = 6 was selected as the largest "
           "stable decomposition.", width=5.5)

# 3.2
add_heading("3.2 Specificity filtering yields biologically interpretable state markers", level=2)
add_para(
    "After applying the specificity ratio filter (≥2.0×), each state's top 20 marker genes "
    "showed clear biological coherence (Figure 2). STATE_01 was enriched for oxidative stress "
    "and Wnt signaling genes (MSRB2, CYGB, APCDD1; specificity ratios 2.2–4.7). STATE_02 "
    "captured interferon response and extracellular matrix remodeling (IFITM3, AHR, LAMA4). "
    "STATE_03 represented a melanocyte dedifferentiation program (LMNA, GPNMB, VIM, MITF). "
    "STATE_04 encompassed protein folding and proliferative machinery (CCT5, CCT2, SUMO1). "
    "STATE_05 featured ribosome biogenesis and mitochondrial metabolism (POLR1D, NDUFAB1). "
    "STATE_06 contained signaling and secretory markers (EVA1B, PLCG2). Notably, no "
    "housekeeping genes appeared among any state's top markers."
)

add_figure("fig2_markers.png",
           "Figure 2. Top 10 marker genes per CIRCUIT state ranked by specificity ratio "
           "(loading in state / mean loading across all states). Darker bars indicate higher "
           "specificity. Gene names are shown for each state.", width=6.2)

# 3.3
add_heading("3.3 Two resistance-associated states significantly outperform random signatures", level=2)
add_para(
    "We scored all six CIRCUIT states and seven benchmark signatures in 76 melanoma patients "
    "across two independent cohorts (GSE78220 and GSE91061). In pooled univariable analysis, "
    "STATE_01 and STATE_05 achieved the highest predictive accuracy (AUROC = 0.729 for both), "
    "with higher scores in non-responders (Figure 3A, Table 3). Both states significantly "
    "outperformed 1,000 random gene signatures of identical size in permutation testing "
    "(STATE_01: P = 0.049; STATE_05: P = 0.047; Figure 3B). Among benchmark signatures, the "
    "best-performing was IFNG antigen presentation (AUROC = 0.616). Direction stability "
    "analysis showed that four of six CIRCUIT states exhibited perfectly consistent direction "
    "of association across cohorts (100% agreement)."
)

add_figure("fig3_auroc_permutation.png",
           "Figure 3. Predictive performance of CIRCUIT states versus benchmark signatures. "
           "(A) Forest plot showing AUROC with 95% bootstrap CI for each signature in pooled "
           "melanoma validation. Blue = CIRCUIT states; red = benchmarks. (B) Permutation null "
           "distribution for STATE_01 (1000 random gene sets); observed AUROC = 0.729, "
           "P = 0.049.", width=6.0)

# 3.4
add_heading("3.4 Multi-state model outperforms benchmarks in cross-validation", level=2)
add_para(
    "A multi-state elastic-net logistic regression combining all six CIRCUIT state scores was "
    "evaluated using repeated stratified 5-fold cross-validation with 50 repeats (250 total "
    "folds) on the pooled 76-sample melanoma cohort (25 responders, 51 non-responders). The "
    "CIRCUIT state model achieved a mean AUROC of 0.694 (SD = 0.048; 95% CI: 0.600–0.761), "
    "mean AUPRC of 0.601 (SD = 0.046), and mean Brier score of 0.256 (SD = 0.026) "
    "(Figure 4A, Table 4). For comparison, an identical cross-validation procedure using only "
    "the seven benchmark signatures yielded a mean AUROC of 0.560 (SD = 0.052), representing "
    "a 0.134 AUROC improvement for the CIRCUIT model. A combined model using both CIRCUIT "
    "states and benchmark signatures achieved AUROC of 0.604, suggesting that benchmark "
    "features introduce noise when combined with CIRCUIT states. Calibration analysis showed "
    "reasonable agreement between predicted probabilities and observed response rates "
    "(Figure 4B)."
)

add_figure("fig4_cv_calibration.png",
           "Figure 4. Cross-validated prediction and calibration. (A) Repeated 5-fold CV AUROC "
           "(50 repeats) for three models: CIRCUIT states only (blue), benchmark signatures only "
           "(red), and combined (green). Error bars show 95% CI. (B) Calibration curve showing "
           "predicted versus observed response fractions across five quantile bins.", width=6.0)

# 3.5
add_heading("3.5 Direction stability supports reproducible cohort effects", level=2)
add_para(
    "Although the two melanoma validation cohorts differ in size and responder prevalence, "
    "the dominant direction of association was consistent for most CIRCUIT states. STATE_01, "
    "STATE_02, STATE_04, and STATE_05 each showed 100% direction agreement across cohorts, "
    "whereas STATE_03 and STATE_06 showed 50% agreement. This pattern supports the "
    "reproducibility of the strongest resistance-associated signals while highlighting that "
    "weaker states remain sensitive to cohort composition and sampling noise."
)

add_figure("fig5_direction_stability.png",
           "Figure 5. Direction stability across melanoma validation cohorts. Bars show the "
           "fraction of cohorts supporting the majority direction for each CIRCUIT state. Blue "
           "bars indicate stable direction (≥80% agreement).", width=4.2)

# ════════════════════════════════════════════════════════════════
# TABLES
# ════════════════════════════════════════════════════════════════
add_heading("Tables", level=1)

# Table 1
add_heading("Table 1. NMF stability analysis across component counts (k).", level=3)
t1 = doc.add_table(rows=8, cols=4, style="Light Shading Accent 1")
headers = ["k", "Reconstruction error", "Stability", "Selected"]
for i, h in enumerate(headers):
    t1.rows[0].cells[i].text = h
data = [
    ["2", "288.18", "1.000", ""],
    ["3", "270.82", "1.000", ""],
    ["4", "256.80", "1.000", ""],
    ["5", "245.26", "1.000", ""],
    ["6", "235.99", "1.000", "Yes"],
    ["7", "226.87", "0.669", ""],
    ["8", "218.48", "0.663", ""],
]
for r, row_data in enumerate(data):
    for c, val in enumerate(row_data):
        t1.rows[r + 1].cells[c].text = val
add_para("Stability is the pairwise assignment agreement across 10 replicate NMF runs.", italic=True, size=8)

# Table 3
uni = pd.read_csv(VAL / "melanoma_univariable_ci.tsv", sep="\t")
pvals = pd.read_csv(VAL / "melanoma_random_signature_pvalues.tsv", sep="\t")
pooled_uni = uni.loc[uni["cohort"] == "pooled"].copy()
pooled_uni = pooled_uni.sort_values("auroc_best_direction", ascending=False).reset_index(drop=True)
pooled_uni = pooled_uni.merge(pvals[["state", "permutation_p"]], on="state", how="left")

add_heading("Table 3. Univariable predictive performance in pooled melanoma validation (n = 76).", level=3)
t3 = doc.add_table(rows=len(pooled_uni) + 1, cols=5, style="Light Shading Accent 1")
headers = ["Signature", "AUROC", "Direction", "95% CI", "Perm. P"]
for i, h in enumerate(headers):
    t3.rows[0].cells[i].text = h
for r, row in pooled_uni.iterrows():
    perm = "—" if pd.isna(row.get("permutation_p")) else f"{row['permutation_p']:.3f}" + ("*" if row["permutation_p"] < 0.05 else "")
    direction = "Higher in NR" if row["direction"] == "higher_in_non_responders" else "Higher in R"
    row_data = [
        signature_display_name(str(row["state"])),
        f"{row['auroc_best_direction']:.3f}",
        direction,
        format_auc_ci(float(row["auroc_ci_low"]), float(row["auroc_ci_high"])),
        perm,
    ]
    for c, val in enumerate(row_data):
        t3.rows[r + 1].cells[c].text = val
add_para("NR = non-responders; R = responders. AUROC reported in the best direction. *P < 0.05 by permutation testing against 1,000 random signatures.", italic=True, size=8)

# Table 4
cv_circuit = pd.read_csv(VAL / "melanoma_cv_summary.tsv", sep="\t").iloc[0]
cv_bench = pd.read_csv(VAL / "melanoma_cv_benchmark_summary.tsv", sep="\t").iloc[0]
cv_combined = pd.read_csv(VAL / "melanoma_cv_combined_summary.tsv", sep="\t").iloc[0]

add_heading("Table 4. Cross-validated predictive performance comparison (repeated 5-fold CV, 50 repeats).", level=3)
t4 = doc.add_table(rows=4, cols=5, style="Light Shading Accent 1")
headers = ["Model", "AUROC (mean ± SD)", "95% CI", "AUPRC", "Brier"]
for i, h in enumerate(headers):
    t4.rows[0].cells[i].text = h
rows = [
    [
        "CIRCUIT states (6)",
        format_mean_sd(float(cv_circuit["cv_auroc_mean"]), float(cv_circuit["cv_auroc_std"])),
        format_auc_ci(float(cv_circuit["cv_auroc_ci_low"]), float(cv_circuit["cv_auroc_ci_high"])),
        format_mean_sd(float(cv_circuit["cv_auprc_mean"]), float(cv_circuit["cv_auprc_std"])),
        format_mean_sd(float(cv_circuit["cv_brier_mean"]), float(cv_circuit["cv_brier_std"])),
    ],
    [
        "Benchmark (7)",
        format_mean_sd(float(cv_bench["cv_auroc_mean"]), float(cv_bench["cv_auroc_std"])),
        format_auc_ci(float(cv_bench["cv_auroc_ci_low"]), float(cv_bench["cv_auroc_ci_high"])),
        format_mean_sd(float(cv_bench["cv_auprc_mean"]), float(cv_bench["cv_auprc_std"])),
        format_mean_sd(float(cv_bench["cv_brier_mean"]), float(cv_bench["cv_brier_std"])),
    ],
    [
        "Combined (13)",
        format_mean_sd(float(cv_combined["cv_auroc_mean"]), float(cv_combined["cv_auroc_std"])),
        format_auc_ci(float(cv_combined["cv_auroc_ci_low"]), float(cv_combined["cv_auroc_ci_high"])),
        format_mean_sd(float(cv_combined["cv_auprc_mean"]), float(cv_combined["cv_auprc_std"])),
        format_mean_sd(float(cv_combined["cv_brier_mean"]), float(cv_combined["cv_brier_std"])),
    ],
]
for r, row_data in enumerate(rows):
    for c, val in enumerate(row_data):
        t4.rows[r + 1].cells[c].text = val

# ════════════════════════════════════════════════════════════════
# 4 DISCUSSION
# ════════════════════════════════════════════════════════════════
add_heading("4 Discussion", level=1)

add_para(
    "We presented CIRCUIT-ICB, a computational framework for discovering and validating "
    "malignant-cell-intrinsic transcriptional states as predictors of ICB resistance. Applied "
    "to melanoma, the framework identified two NMF-derived resistance states (STATE_01 and "
    "STATE_05) that achieved AUROC of 0.729 in pooled validation—surpassing the established "
    "IFNG antigen-presentation benchmark (AUROC = 0.616)—and a multi-state cross-validated "
    "AUROC of 0.694, outperforming benchmark signatures (AUROC = 0.560) by a margin of 0.134."
)

add_para(
    "The biological coherence of the top-performing states supports their mechanistic "
    "relevance. STATE_01 is enriched for oxidative stress response genes (MSRB2, methionine "
    "sulfoxide reductase B2; CYGB, cytoglobin) and Wnt pathway modulation (APCDD1). Oxidative "
    "stress has been implicated in immune evasion through upregulation of immunosuppressive "
    "metabolites and modulation of antigen presentation (Aggarwal et al., 2019). STATE_05 "
    "captures ribosome biogenesis (POLR1D) and mitochondrial metabolism (NDUFAB1, COX6B1), "
    "reflecting a hypermetabolic state previously associated with therapy resistance through "
    "enhanced survival signaling and reduced sensitivity to immune-mediated killing "
    "(Haq et al., 2013)."
)

add_para(
    "A key methodological contribution is the demonstration that conventional NMF applied to "
    "scRNA-seq-derived pseudobulk data produces signatures dominated by housekeeping genes "
    "that fail to outperform random gene sets. Our systematic exclusion of 120+ housekeeping "
    "genes, combined with the specificity-ratio filter requiring ≥2-fold enrichment relative "
    "to the cross-state mean, was essential for recovering biologically meaningful programs. "
    "This preprocessing pipeline should be useful in related NMF-based decomposition workflows."
)

add_para(
    "The stability-guided component selection procedure provides an objective alternative to "
    "reconstruction-error-based methods, which monotonically decrease with increasing k and "
    "offer no natural stopping criterion. By requiring consistent sample-to-state assignment "
    "across replicate factorizations, we selected k = 6 as the largest stable decomposition. "
    "The sharp stability drop at k = 7 (from 1.0 to 0.669) suggests that the seventh "
    "component fragments a coherent program into unstable sub-programs."
)

add_para(
    "Several limitations should be acknowledged. First, the validation cohort comprises only "
    "76 melanoma patients from two studies, limiting statistical power and generalizability. "
    "The cross-validated AUROC of 0.694 (95% CI: 0.600–0.761) is promising but requires "
    "confirmation in larger cohorts. Two controlled-access datasets—Liu et al. (2019; n = 121; "
    "dbGaP phs000452) and Gide et al. (2019; n = 120; EGA EGAS00001003122)—would more than "
    "triple the validation sample size and are a priority for future work. Second, our "
    "framework operates on pre-treatment bulk RNA-seq and does not capture dynamic resistance "
    "mechanisms that emerge during treatment."
)

add_para(
    "In conclusion, CIRCUIT-ICB provides a rigorous, reproducible framework for translating "
    "single-cell-resolution malignant cell programs into bulk RNA-seq-based ICB response "
    "predictors. The melanoma-specific resistance states identified here represent candidates "
    "for prospective validation and potential integration into clinical decision-support tools."
)

# ════════════════════════════════════════════════════════════════
# DATA AVAILABILITY
# ════════════════════════════════════════════════════════════════
add_heading("Data Availability", level=1)
add_para(
    "All scRNA-seq pseudobulk data used for NMF discovery are available from Zenodo. Bulk "
    "RNA-seq validation datasets are available from the Gene Expression Omnibus: GSE78220 and "
    "GSE91061. Source code for the complete CIRCUIT-ICB pipeline is publicly available at https://github.com/Nigmat-future/20260426-circuit-icb-pan-cancer-immune-checkpoint-blockade-resistance-cell-state-atlas."
)

# ════════════════════════════════════════════════════════════════
# FUNDING & CONFLICT
# ════════════════════════════════════════════════════════════════
add_heading("Funding", level=1)
add_para("This work received no specific funding.")

add_heading("Conflict of Interest", level=1)
add_para("The author declares no competing interests.")

# ════════════════════════════════════════════════════════════════
# REFERENCES
# ════════════════════════════════════════════════════════════════
add_heading("References", level=1)

refs = [
    "Aggarwal, V., Tuli, H.S., Varol, A., et al. (2019) Role of reactive oxygen species in cancer progression: molecular mechanisms and recent advancements. Biomolecules, 9, 735.",
    "Ayers, M., Lunceford, J., Nebozhyn, M., et al. (2017) IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade. Journal of Clinical Investigation, 127, 2930–2940.",
    "Cristescu, R., Mogg, R., Ayers, M., et al. (2018) Pan-tumor genomic biomarkers for PD-1 checkpoint blockade-based immunotherapy. Science, 362, eaar3593.",
    "Gide, T.N., Quek, C., Menzies, A.M., et al. (2019) Distinct immune cell populations define response to anti-PD-1 monotherapy and anti-PD-1/anti-CTLA-4 combined therapy. Cancer Cell, 35, 238–255.",
    "Haq, R., Shoag, J., Andreu-Perez, P., et al. (2013) Oncogenic BRAF regulates oxidative metabolism via PGC1α and MITF. Cancer Cell, 23, 302–315.",
    "Hugo, W., Zaretsky, J.M., Sun, L., et al. (2016) Genomic and transcriptomic features of response to anti-PD-1 therapy in metastatic melanoma. Cell, 165, 35–44.",
    "Jerby-Arnon, L., Shah, P., Cuoco, M.S., et al. (2018) A cancer cell program promotes T cell exclusion and resistance to checkpoint blockade. Cell, 175, 984–997.",
    "Kotliar, D., Veres, A., Nagy, M.A., et al. (2019) Identifying gene expression programs of cell-type identity and cellular activity with single-cell RNA-Seq. eLife, 8, e43803.",
    "Larkin, J., Chiarion-Sileni, V., Gonzalez, R., et al. (2015) Combined nivolumab and ipilimumab or monotherapy in untreated melanoma. New England Journal of Medicine, 373, 23–34.",
    "Liu, D., Schilling, B., Liu, D., et al. (2019) Integrative molecular and clinical modeling of clinical outcomes to PD1 blockade in patients with metastatic melanoma. Nature Medicine, 25, 1916–1927.",
    "Pozniak, J., Riesenberg, S., Posber, M., et al. (2022) A TCF4-dependent gene regulatory network confers resistance to immunotherapy in melanoma. Cell, 185, 166–183.",
    "Riaz, N., Havel, J.J., Makarov, V., et al. (2017) Tumor and microenvironment evolution during immunotherapy with nivolumab. Cell, 171, 934–949.",
    "Robert, C., Schachter, J., Long, G.V., et al. (2015) Pembrolizumab versus ipilimumab in advanced melanoma. New England Journal of Medicine, 372, 2521–2532.",
    "Rooney, M.S., Shukla, S.A., Wu, C.J., et al. (2015) Molecular and genetic properties of tumors associated with local immune cytolytic activity. Cell, 160, 48–61.",
    "Sharma, P., Hu-Lieskovan, S., Wargo, J.A. and Ribas, A. (2017) Primary, adaptive, and acquired resistance to cancer immunotherapy. Cell, 168, 707–723.",
    "Tirosh, I., Izar, B., Prakadan, S.M., et al. (2016) Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science, 352, 189–196.",
]

for i, ref in enumerate(refs, 1):
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(2)
    run = p.add_run(f"[{i}] {ref}")
    run.font.name = "Times New Roman"
    run.font.size = Pt(9)

# ════════════════════════════════════════════════════════════════
# SAVE
# ════════════════════════════════════════════════════════════════
out_path = ROOT / "manuscript" / "CIRCUIT_ICB_manuscript.docx"
doc.save(str(out_path))
print(f"Manuscript saved to {out_path}")
