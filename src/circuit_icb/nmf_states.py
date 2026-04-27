from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.decomposition import NMF


DEFAULT_EXCLUDE_PREFIXES = (
    "MT-",
    "RPL",
    "RPS",
    "MRPL",
    "MRPS",
    "MTRNR",
    "HB",
    "HIST",
    "LINC",
    "SNOR",
    "MIR",
)

DEFAULT_EXCLUDE_GENES = {
    # lncRNA / antisense / novel transcripts
    "MALAT1", "NEAT1", "XIST", "KCNQ1OT1", "SNHG16", "UGDH-AS1",
    "C17orf76-AS1", "FP236383.2", "FP236383.3", "FP671120.3",
    # housekeeping / cytoskeleton
    "GAPDH", "ACTB", "ACTG1", "TUBB", "TUBA1B", "TUBA1A",
    # translation / ribosome-associated
    "EEF1A1", "EEF1G", "EEF1B2", "EEF1D", "EEF2",
    "FAU", "UBA52", "RACK1", "GNB2L1",
    "NACA", "NAP1L1", "NPM1", "PABPC1",
    "EIF4A1", "EIF3C", "EIF3L", "BTF3",
    # ubiquitin / thymosin
    "UBB", "UBC", "TMSB4X", "TMSB10", "PTMA",
    # ferritin
    "FTH1", "FTL",
    # chaperones / heat shock
    "HSP90AA1", "HSP90AB1", "HSPA8", "HSPA1A", "HSPA1B",
    "HSPB1", "HSPE1", "DNAJB1",
    # small structural
    "PFN1", "CFL1", "TPT1", "MIF", "SERF2", "TOMM7", "ELOB",
    # antigen presentation (scored separately in benchmarks)
    "B2M", "HLA-A", "HLA-B", "HLA-C",
    # glycolysis / metabolic housekeeping
    "ALDOA", "ALDOB", "TPI1", "ENO1", "PGK1", "LDHA", "LDHB",
    "PKM", "PGAM1",
    # mitochondrial electron transport
    "ATP5F1E", "ATP5F1B", "ATP5MC2", "ATP5B", "ATP5A1", "ATP5E",
    "UQCRQ", "UQCRB", "NDUFA4", "COX7C", "COX6C", "SLC25A3", "SLC25A6",
    # histones
    "H3F3A", "H3F3B", "H2AFZ",
    # calmodulin / signalling scaffolds
    "PPIA", "PPIB", "CALM1", "CALM2", "CALM3", "YBX1", "DDX5",
    "MYL6", "MYL12A", "MYL12B",
    "YWHAZ", "YWHAB", "YWHAE",
    # redox / detox
    "GSTP1", "PRDX1", "PRDX2", "SOD1",
    # annexins / S100 generic
    "ANXA1", "ANXA2", "ANXA5",
    "S100A6", "S100A10", "S100A11",
    # immediate-early / stress-response (non-specific)
    "FOS", "JUN", "JUNB", "ZFP36", "IER3", "NFKBIA", "RHOB",
    "GNAS", "GNB1", "GNB2",
    # other pan-tissue high-expression
    "TAGLN2", "PPDPF", "CD63", "CD81",
}

METADATA_START_MARKERS = [
    "sample_id_pre_post",
    "sample_id",
    "cell_id",
    "orig.ident",
    "nCount_RNA",
    "nFeature_RNA",
]


@dataclass(frozen=True)
class NMFResult:
    activity: pd.DataFrame
    loadings: pd.DataFrame
    markers: pd.DataFrame
    metrics: pd.DataFrame
    feature_genes: list[str]


def dataset_id_from_path(path: Path) -> str:
    return path.stem.removesuffix("_data")


def find_gene_columns(columns: list[str]) -> list[str]:
    """Return gene-expression columns using metadata boundary markers."""

    marker_positions = [columns.index(marker) for marker in METADATA_START_MARKERS if marker in columns]
    if not marker_positions:
        raise ValueError("Could not find metadata boundary markers in expression table.")
    start_meta = min(marker_positions)
    return [column for column in columns[1:start_meta] if column and not str(column).startswith("Unnamed")]


def load_malignant_pseudobulk(directory: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load malignant-only pseudobulk CSV files as expression and metadata tables."""

    directory = Path(directory)
    files = sorted(path for path in directory.rglob("*.csv") if path.name != "Combined_data.csv")
    if not files:
        raise FileNotFoundError(f"No malignant pseudobulk CSV files found under {directory}")

    expression_parts = []
    metadata_parts = []
    gene_columns_reference: list[str] | None = None

    for path in files:
        dataset = dataset_id_from_path(path)
        df = pd.read_csv(path, low_memory=False)
        header = df.columns.tolist()
        gene_columns = find_gene_columns(header)
        start_meta = min(header.index(marker) for marker in METADATA_START_MARKERS if marker in header)

        if gene_columns_reference is None:
            gene_columns_reference = gene_columns
        else:
            if gene_columns_reference != gene_columns:
                raise ValueError(f"{path.name} gene columns differ from the reference file.")

        sample_ids = [f"{dataset}__{idx}" for idx in df[header[0]].astype(str)]
        expr = df.iloc[:, 1:start_meta].astype("float32", copy=False)
        expr.columns = gene_columns_reference
        expr.index = sample_ids
        expression_parts.append(expr)

        meta = pd.concat([df.iloc[:, [0]], df.iloc[:, start_meta:]], axis=1).copy()
        rename = {}
        if "sample" in meta.columns:
            rename["sample"] = "source_sample"
        if "dataset" in meta.columns:
            rename["dataset"] = "source_dataset"
        if rename:
            meta = meta.rename(columns=rename)
        meta.insert(0, "sample", sample_ids)
        meta.insert(1, "dataset", dataset)
        metadata_parts.append(meta)

    expression = pd.concat(expression_parts, axis=0)
    metadata = pd.concat(metadata_parts, axis=0, ignore_index=True)
    return expression, metadata


def select_variable_genes(
    expression: pd.DataFrame,
    *,
    n_top_genes: int = 3000,
    min_positive_fraction: float = 0.05,
    exclude_prefixes: tuple[str, ...] = DEFAULT_EXCLUDE_PREFIXES,
    exclude_genes: set[str] = DEFAULT_EXCLUDE_GENES,
) -> list[str]:
    """Select high-variance genes after basic presence filtering."""

    expr = expression.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    positive_fraction = (expr > 0).mean(axis=0)
    keep = positive_fraction >= min_positive_fraction
    gene_names = pd.Index(expr.columns.astype(str))
    excluded = gene_names.str.upper().isin(exclude_genes)
    for prefix in exclude_prefixes:
        excluded = excluded | gene_names.str.upper().str.startswith(prefix)
    keep = keep & ~pd.Series(excluded, index=expr.columns)
    filtered = expr.loc[:, keep]
    variances = filtered.var(axis=0).sort_values(ascending=False)
    return variances.head(min(n_top_genes, len(variances))).index.tolist()


def fit_nmf_states(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    n_components: int = 7,
    n_top_genes: int = 3000,
    marker_top_n: int = 20,
    random_state: int = 7,
    max_iter: int = 1000,
    specificity_ratio: float = 2.0,
) -> NMFResult:
    """Fit a malignant-cell state NMF model and return activity and markers."""

    feature_genes = select_variable_genes(expression, n_top_genes=n_top_genes)
    matrix = _prepare_nmf_matrix(expression, feature_genes)

    model = NMF(
        n_components=n_components,
        init="nndsvda",
        solver="cd",
        beta_loss="frobenius",
        max_iter=max_iter,
        random_state=random_state,
    )
    w = model.fit_transform(matrix)
    h = model.components_
    state_names = [f"CIRCUIT_STATE_{idx:02d}" for idx in range(1, n_components + 1)]

    activity = pd.DataFrame(w, columns=[f"{state}_activity" for state in state_names])
    activity.insert(0, "sample", expression.index)
    activity = activity.merge(metadata, on="sample", how="left")

    loadings = pd.DataFrame(h, index=state_names, columns=feature_genes)
    markers = markers_from_loadings(
        loadings, top_n=marker_top_n, specificity_ratio=specificity_ratio,
    )
    metrics = pd.DataFrame(
        [
            {
                "n_components": n_components,
                "n_samples": expression.shape[0],
                "n_features_selected": len(feature_genes),
                "reconstruction_error": model.reconstruction_err_,
                "n_iter": model.n_iter_,
                "random_state": random_state,
            }
        ]
    )
    return NMFResult(activity=activity, loadings=loadings, markers=markers, metrics=metrics, feature_genes=feature_genes)


def _prepare_nmf_matrix(expression: pd.DataFrame, feature_genes: list[str]) -> np.ndarray:
    """Shared preprocessing: log1p then per-gene z-score shifted to non-negative."""
    matrix = expression[feature_genes].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float)
    matrix = np.clip(matrix, a_min=0.0, a_max=None)
    matrix = np.log1p(matrix)
    gene_mean = matrix.mean(axis=0)
    gene_std = matrix.std(axis=0)
    gene_std[gene_std == 0] = 1.0
    matrix = (matrix - gene_mean) / gene_std
    matrix = matrix - matrix.min(axis=0)
    return matrix


def scan_nmf_k(
    expression: pd.DataFrame,
    *,
    k_values: list[int],
    n_top_genes: int = 3000,
    random_state: int = 7,
    max_iter: int = 700,
    n_runs: int = 5,
) -> pd.DataFrame:
    """Scan NMF component counts with stability analysis.

    Runs NMF *n_runs* times per k with different seeds and reports mean
    reconstruction error plus cophenetic-style stability (pairwise
    agreement of sample-to-dominant-state assignments across runs).
    """

    feature_genes = select_variable_genes(expression, n_top_genes=n_top_genes)
    matrix = _prepare_nmf_matrix(expression, feature_genes)

    rows = []
    for k in k_values:
        assignments_list = []
        errors = []
        iters = []
        for run in range(n_runs):
            seed = random_state + run
            model = NMF(
                n_components=k,
                init="nndsvda",
                solver="cd",
                beta_loss="frobenius",
                max_iter=max_iter,
                random_state=seed,
            )
            w = model.fit_transform(matrix)
            errors.append(model.reconstruction_err_)
            iters.append(model.n_iter_)
            assignments_list.append(w.argmax(axis=1))

        stability = _pairwise_assignment_stability(assignments_list)
        rows.append(
            {
                "n_components": k,
                "n_samples": expression.shape[0],
                "n_features_selected": len(feature_genes),
                "reconstruction_error": float(np.mean(errors)),
                "reconstruction_error_std": float(np.std(errors)),
                "stability": stability,
                "n_iter": int(np.mean(iters)),
                "n_runs": n_runs,
                "random_state": random_state,
            }
        )
    return pd.DataFrame(rows)


def _pairwise_assignment_stability(assignments: list[np.ndarray]) -> float:
    """Fraction of sample pairs with identical dominant-state across runs."""
    n_runs = len(assignments)
    if n_runs < 2:
        return 1.0
    n_samples = len(assignments[0])
    agree = 0
    total = 0
    for i in range(n_runs):
        for j in range(i + 1, n_runs):
            agree += int(np.sum(assignments[i] == assignments[j]))
            total += n_samples
    return agree / total if total else 1.0


def markers_from_loadings(
    loadings: pd.DataFrame,
    top_n: int = 20,
    *,
    specificity_ratio: float = 2.0,
) -> pd.DataFrame:
    """Extract top marker genes per state with a specificity filter.

    A gene is kept only if its loading in this state is at least
    *specificity_ratio* times the mean loading for that gene across all
    states. This eliminates ubiquitous genes shared across programs.
    """
    mean_loading = loadings.mean(axis=0)
    rows = []
    for state, values in loadings.iterrows():
        if specificity_ratio > 0:
            specificity = values / mean_loading.replace(0, np.nan)
            specific = specificity.loc[specificity >= specificity_ratio].index
            candidates = values.loc[values.index.isin(specific)]
        else:
            candidates = values
        top = candidates.sort_values(ascending=False).head(top_n)
        for rank, (gene, weight) in enumerate(top.items(), start=1):
            rows.append({
                "state": state,
                "rank": rank,
                "gene": gene,
                "weight": weight,
                "specificity": values[gene] / mean_loading[gene] if mean_loading[gene] > 0 else float("inf"),
            })
    return pd.DataFrame(rows)


LINEAGE_GENE_SETS: dict[str, set[str]] = {
    "liver": {"ALB", "SERPINA1", "APOC1", "APOA2", "AMBP", "FGG", "FGA", "VTN",
              "ORM1", "AGT", "TTR", "RBP4", "SEPP1", "APOE", "C3", "RARRES2"},
    "breast_luminal": {"SCGB2A2", "GATA3", "MUCL1", "AGR2", "TFF3", "AZGP1",
                       "SPINT2", "ADIRF", "XBP1"},
    "keratin_epithelial": {"KRT14", "KRT5", "KRT17", "KRT15", "KRT6A", "KRT6B",
                           "DSP", "SFN", "PERP", "S100A2", "S100A7", "S100A8",
                           "S100A9", "S100A14"},
    "melanocyte": {"PMEL", "MLANA", "PRAME", "GPNMB", "DCT", "TYR", "TYRP1",
                   "MITF", "SOX10"},
}


def detect_lineage_states(
    markers: pd.DataFrame,
    *,
    lineage_gene_sets: dict[str, set[str]] | None = None,
    overlap_threshold: float = 0.3,
) -> pd.DataFrame:
    """Flag states where >*overlap_threshold* of markers are lineage-specific genes."""
    if lineage_gene_sets is None:
        lineage_gene_sets = LINEAGE_GENE_SETS
    rows = []
    for state, group in markers.groupby("state"):
        genes = set(group["gene"].str.upper())
        n = len(genes)
        best_lineage = ""
        best_frac = 0.0
        for lineage, gene_set in lineage_gene_sets.items():
            overlap = len(genes & gene_set)
            frac = overlap / n if n else 0.0
            if frac > best_frac:
                best_frac = frac
                best_lineage = lineage
        rows.append({
            "state": state,
            "lineage_class": best_lineage if best_frac >= overlap_threshold else "pan_cancer",
            "lineage_overlap_fraction": best_frac,
            "is_lineage_specific": best_frac >= overlap_threshold,
        })
    return pd.DataFrame(rows)


def filter_lineage_states(
    result: NMFResult,
    *,
    overlap_threshold: float = 0.3,
) -> NMFResult:
    """Remove lineage-specific states from an NMFResult."""
    lineage_info = detect_lineage_states(
        result.markers, overlap_threshold=overlap_threshold,
    )
    keep_states = lineage_info.loc[~lineage_info["is_lineage_specific"], "state"].tolist()
    if not keep_states:
        return result
    activity_cols = ["sample"] + [col for col in result.activity.columns
                                  if col != "sample" and any(col.startswith(s) for s in keep_states)]
    other_cols = [col for col in result.activity.columns
                  if col != "sample" and not any(col.startswith(s + "_") for s in
                  lineage_info["state"].tolist())]
    keep_cols = list(dict.fromkeys(["sample"] + activity_cols[1:] + other_cols))
    return NMFResult(
        activity=result.activity[keep_cols],
        loadings=result.loadings.loc[result.loadings.index.isin(keep_states)],
        markers=result.markers.loc[result.markers["state"].isin(keep_states)].reset_index(drop=True),
        metrics=result.metrics,
        feature_genes=result.feature_genes,
    )


def derive_down_genes(
    expression: pd.DataFrame,
    activity: pd.DataFrame,
    state_col: str,
    *,
    n_down: int = 10,
    exclude_genes: set[str] | None = None,
) -> list[str]:
    """Find genes most anti-correlated with a state's activity score."""
    expr = expression.select_dtypes(include="number")
    if exclude_genes:
        expr = expr.drop(columns=[g for g in exclude_genes if g in expr.columns], errors="ignore")
    act = activity[state_col].astype(float)
    common = expr.index.intersection(act.index)
    if len(common) < 10:
        return []
    correlations = expr.loc[common].corrwith(act.loc[common])
    bottom = correlations.dropna().sort_values().head(n_down)
    return bottom.index.tolist()


def write_signature_yaml(
    markers: pd.DataFrame,
    path: str | Path,
    top_n: int = 20,
    *,
    down_genes: dict[str, list[str]] | None = None,
) -> None:
    """Write NMF marker genes as state signatures compatible with scoring code."""

    states: dict[str, dict[str, object]] = {}
    for state, group in markers.groupby("state", sort=False):
        top = group.sort_values("rank").head(top_n)
        down = down_genes.get(state, []) if down_genes else []
        states[state] = {
            "description": f"NMF-derived malignant-cell program {state}",
            "up": top["gene"].astype(str).tolist(),
            "down": down,
        }
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        yaml.safe_dump({"states": states}, handle, sort_keys=False)
