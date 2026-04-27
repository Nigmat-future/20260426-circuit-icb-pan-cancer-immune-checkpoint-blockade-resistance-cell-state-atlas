from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


def read_hpa_tissue_consensus(path: str | Path) -> pd.DataFrame:
    """Read HPA RNA tissue consensus data from a TSV or zipped TSV file."""

    path = Path(path)
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as zf:
            names = [name for name in zf.namelist() if name.endswith(".tsv")]
            if not names:
                raise ValueError(f"No TSV file found in {path}")
            with zf.open(names[0]) as handle:
                table = pd.read_csv(handle, sep="\t")
    else:
        table = pd.read_csv(path, sep="\t")

    required = {"Gene name", "Tissue", "nTPM"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"HPA table is missing columns: {sorted(missing)}")
    table = table.rename(columns={"Gene name": "gene", "Tissue": "tissue", "nTPM": "nTPM"})
    table["gene"] = table["gene"].astype(str)
    table["tissue"] = table["tissue"].astype(str)
    table["nTPM"] = pd.to_numeric(table["nTPM"], errors="coerce")
    return table[["gene", "tissue", "nTPM"]].dropna(subset=["gene", "tissue"])


def normal_expression_safety_score(
    max_ntpm: pd.Series,
    broad_fraction: pd.Series,
    *,
    max_ntpm_cap: float = 100.0,
    max_expression_weight: float = 0.7,
) -> pd.Series:
    """Convert normal tissue expression burden into a 0-1 safety proxy.

    This is a conservative screening score, not a clinical toxicity estimate.
    High maximum normal expression and broad tissue expression both reduce the
    score because they increase the risk that a target is less tumor-selective.
    """

    if max_ntpm_cap <= 0:
        raise ValueError("max_ntpm_cap must be positive")
    if not 0 <= max_expression_weight <= 1:
        raise ValueError("max_expression_weight must be between 0 and 1")

    max_component = np.log1p(max_ntpm.clip(lower=0)) / np.log1p(max_ntpm_cap)
    max_component = max_component.clip(lower=0, upper=1)
    broad_component = broad_fraction.fillna(1.0).clip(lower=0, upper=1)
    burden = max_expression_weight * max_component + (1 - max_expression_weight) * broad_component
    return (1.0 - burden).clip(lower=0, upper=1)


def build_hpa_safety_table(
    hpa: pd.DataFrame,
    genes: list[str],
    *,
    expression_threshold: float = 1.0,
    max_ntpm_cap: float = 100.0,
) -> pd.DataFrame:
    """Summarize HPA normal expression and produce one safety row per gene."""

    requested = pd.DataFrame({"gene": pd.Series(genes, dtype=str).drop_duplicates()})
    subset = hpa.loc[hpa["gene"].isin(requested["gene"])].copy()
    if subset.empty:
        requested["safety"] = np.nan
        requested["hpa_max_nTPM"] = np.nan
        requested["hpa_median_nTPM"] = np.nan
        requested["hpa_broad_expression_fraction"] = np.nan
        requested["hpa_highest_tissue"] = ""
        return requested

    subset["nTPM"] = pd.to_numeric(subset["nTPM"], errors="coerce")
    grouped = subset.groupby("gene", sort=False)
    max_ntpm = grouped["nTPM"].max()
    median_ntpm = grouped["nTPM"].median()
    broad_fraction = grouped["nTPM"].apply(lambda values: float((values >= expression_threshold).mean()))

    highest_rows = (
        subset.sort_values(["gene", "nTPM", "tissue"], ascending=[True, False, True])
        .drop_duplicates("gene")
        .set_index("gene")["tissue"]
    )

    summary = pd.DataFrame(
        {
            "gene": max_ntpm.index,
            "hpa_max_nTPM": max_ntpm.values,
            "hpa_median_nTPM": max_ntpm.index.map(median_ntpm).astype(float),
            "hpa_broad_expression_fraction": max_ntpm.index.map(broad_fraction).astype(float),
            "hpa_highest_tissue": max_ntpm.index.map(highest_rows).fillna(""),
        }
    )
    summary["safety"] = normal_expression_safety_score(
        summary["hpa_max_nTPM"],
        summary["hpa_broad_expression_fraction"],
        max_ntpm_cap=max_ntpm_cap,
    )

    return requested.merge(summary, on="gene", how="left")
