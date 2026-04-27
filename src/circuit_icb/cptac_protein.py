from __future__ import annotations

import numpy as np
import pandas as pd


def build_cptac_protein_evidence(
    protein: pd.DataFrame,
    genes: list[str],
    *,
    gene_col: str = "gene",
    value_col: str = "protein_evidence",
) -> pd.DataFrame:
    """Normalize a CPTAC/PDC protein evidence table into target-priority input.

    The preferred input is one row per gene with a precomputed protein_evidence
    score in [0, 1]. If a quantitative abundance/fold-change column is supplied,
    values are min-max scaled over the candidate set as a transparent adapter.
    """

    if gene_col not in protein.columns:
        raise ValueError(f"Protein table must contain {gene_col}")
    if value_col not in protein.columns:
        raise ValueError(f"Protein table must contain {value_col}")

    requested = pd.DataFrame({"gene": pd.Series(genes, dtype=str).drop_duplicates()})
    subset = protein[[gene_col, value_col]].rename(columns={gene_col: "gene", value_col: "raw_protein_evidence"}).copy()
    subset["raw_protein_evidence"] = pd.to_numeric(subset["raw_protein_evidence"], errors="coerce")
    subset = subset.dropna(subset=["gene"]).groupby("gene", as_index=False)["raw_protein_evidence"].max()
    out = requested.merge(subset, on="gene", how="left")

    values = out["raw_protein_evidence"]
    if values.dropna().empty:
        out["protein_evidence"] = np.nan
    elif values.min() >= 0 and values.max() <= 1:
        out["protein_evidence"] = values
    elif values.max() == values.min():
        out["protein_evidence"] = np.where(values.notna(), 1.0, np.nan)
    else:
        out["protein_evidence"] = (values - values.min()) / (values.max() - values.min())
    return out
