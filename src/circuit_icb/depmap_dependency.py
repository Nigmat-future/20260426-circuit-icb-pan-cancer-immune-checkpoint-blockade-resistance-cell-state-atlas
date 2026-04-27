from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


GENE_COLUMN_RE = re.compile(r"^(?P<symbol>.+?)\s+\(\d+\)$")


def clean_depmap_gene_column(column: str) -> str:
    """Convert DepMap matrix columns such as 'BRAF (673)' to gene symbols."""

    match = GENE_COLUMN_RE.match(str(column))
    if match:
        return match.group("symbol")
    return str(column)


def read_common_essential_genes(path: str | Path | None) -> set[str]:
    if not path:
        return set()
    path = Path(path)
    if not path.exists():
        return set()
    table = pd.read_csv(path)
    if table.empty:
        return set()
    for column in ["gene", "Gene", "gene_name", "Gene Name", "symbol", "Symbol"]:
        if column in table.columns:
            return {clean_depmap_gene_column(value).upper() for value in table[column].dropna().astype(str)}
    first_col = table.columns[0]
    return {clean_depmap_gene_column(value).upper() for value in table[first_col].dropna().astype(str)}


def dependency_columns_for_genes(columns: list[str], genes: list[str]) -> dict[str, list[str]]:
    requested = {gene.upper() for gene in genes}
    mapping: dict[str, list[str]] = {gene.upper(): [] for gene in genes}
    for column in columns:
        clean = clean_depmap_gene_column(column).upper()
        if clean in requested:
            mapping.setdefault(clean, []).append(column)
    return mapping


def build_depmap_dependency_table(
    dependency_matrix: pd.DataFrame,
    model_metadata: pd.DataFrame,
    genes: list[str],
    *,
    model_id_col: str = "ModelID",
    lineage_col: str = "OncotreeLineage",
    dependency_threshold: float = 0.5,
    common_essential_genes: set[str] | None = None,
) -> pd.DataFrame:
    """Summarize DepMap gene dependency probabilities for candidate genes."""

    dependency_matrix = dependency_matrix.copy()
    if model_id_col not in dependency_matrix.columns:
        first_col = dependency_matrix.columns[0]
        dependency_matrix = dependency_matrix.rename(columns={first_col: model_id_col})
    if model_id_col not in model_metadata.columns:
        metadata = pd.DataFrame({model_id_col: dependency_matrix[model_id_col].astype(str)})
    else:
        metadata = model_metadata.copy()

    if lineage_col not in metadata.columns:
        metadata[lineage_col] = "unknown"
    metadata = metadata[[model_id_col, lineage_col]].drop_duplicates(model_id_col)
    merged = dependency_matrix.merge(metadata, on=model_id_col, how="left")
    merged[lineage_col] = merged[lineage_col].fillna("unknown").astype(str)

    column_map = dependency_columns_for_genes(list(dependency_matrix.columns), genes)
    common = common_essential_genes or set()
    rows = []
    for gene in pd.Series(genes, dtype=str).drop_duplicates():
        gene_key = gene.upper()
        columns = column_map.get(gene_key, [])
        if not columns:
            rows.append(
                {
                    "gene": gene,
                    "dependency": pd.NA,
                    "depmap_mean_dependency_probability": pd.NA,
                    "depmap_fraction_dependent": pd.NA,
                    "depmap_max_lineage_median_dependency": pd.NA,
                    "depmap_top_lineage": "",
                    "depmap_models_evaluated": 0,
                    "depmap_common_essential": gene_key in common,
                }
            )
            continue

        values = merged[columns].apply(pd.to_numeric, errors="coerce")
        if values.shape[1] > 1:
            gene_values = values.max(axis=1)
        else:
            gene_values = values.iloc[:, 0]
        valid = gene_values.dropna()
        if valid.empty:
            rows.append(
                {
                    "gene": gene,
                    "dependency": pd.NA,
                    "depmap_mean_dependency_probability": pd.NA,
                    "depmap_fraction_dependent": pd.NA,
                    "depmap_max_lineage_median_dependency": pd.NA,
                    "depmap_top_lineage": "",
                    "depmap_models_evaluated": 0,
                    "depmap_common_essential": gene_key in common,
                }
            )
            continue

        lineage_median = pd.DataFrame({"dependency": gene_values, "lineage": merged[lineage_col]}).dropna()
        lineage_median = lineage_median.groupby("lineage")["dependency"].median().sort_values(ascending=False)
        top_lineage = str(lineage_median.index[0]) if not lineage_median.empty else ""
        max_lineage = float(lineage_median.iloc[0]) if not lineage_median.empty else float("nan")
        fraction_dependent = float((valid >= dependency_threshold).mean())
        mean_dependency = float(valid.mean())
        raw_dependency_score = max(fraction_dependent, max_lineage)
        lineage_selectivity_delta = max(0.0, max_lineage - float(valid.median())) if pd.notna(max_lineage) else 0.0
        dependency_score = raw_dependency_score
        if gene_key in common:
            dependency_score = min(dependency_score, 0.25)

        rows.append(
            {
                "gene": gene,
                "dependency": max(0.0, min(1.0, dependency_score)),
                "depmap_dependency_raw": raw_dependency_score,
                "depmap_mean_dependency_probability": mean_dependency,
                "depmap_fraction_dependent": fraction_dependent,
                "depmap_max_lineage_median_dependency": max_lineage,
                "depmap_lineage_selectivity_delta": lineage_selectivity_delta,
                "depmap_top_lineage": top_lineage,
                "depmap_models_evaluated": int(valid.shape[0]),
                "depmap_common_essential": gene_key in common,
            }
        )
    return pd.DataFrame(rows)
