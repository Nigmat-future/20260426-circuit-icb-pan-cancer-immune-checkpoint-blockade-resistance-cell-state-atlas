from __future__ import annotations

import re

import numpy as np
import pandas as pd


GENE_TOKEN_RE = re.compile(r"[A-Z0-9][A-Z0-9.-]*")


def target_tokens(target: object) -> set[str]:
    if pd.isna(target):
        return set()
    return {token.upper() for token in GENE_TOKEN_RE.findall(str(target).upper())}


def summarize_prism_compounds(prism: pd.DataFrame) -> pd.DataFrame:
    required = {"broad_id", "auc", "target"}
    missing = required.difference(prism.columns)
    if missing:
        raise ValueError(f"PRISM table is missing columns: {sorted(missing)}")

    table = prism.copy()
    table["auc"] = pd.to_numeric(table["auc"], errors="coerce")
    table = table.dropna(subset=["broad_id", "auc", "target"])
    for column in ["name", "moa", "phase"]:
        if column not in table.columns:
            table[column] = ""
    grouped = (
        table.groupby(["broad_id", "name", "moa", "target", "phase"], dropna=False)["auc"]
        .agg(prism_median_auc="median", prism_min_auc="min", prism_n_profiles="count")
        .reset_index()
    )
    grouped = grouped.loc[grouped["prism_n_profiles"] > 0].copy()
    if grouped.empty:
        grouped["prism_compound_sensitivity_score"] = []
        return grouped

    ranks = grouped["prism_median_auc"].rank(method="average", ascending=True)
    denom = max(len(grouped) - 1, 1)
    grouped["prism_compound_sensitivity_score"] = 1.0 - ((ranks - 1.0) / denom)
    grouped["prism_compound_sensitivity_score"] = grouped["prism_compound_sensitivity_score"].clip(0.0, 1.0)
    grouped["target_tokens"] = grouped["target"].map(target_tokens)
    return grouped


def build_prism_sensitivity_table(prism: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    compounds = summarize_prism_compounds(prism)
    rows = []
    for gene in pd.Series(genes, dtype=str).drop_duplicates():
        gene_key = gene.upper()
        if compounds.empty:
            matches = compounds
        else:
            matches = compounds.loc[compounds["target_tokens"].map(lambda tokens: gene_key in tokens)].copy()
        if matches.empty:
            rows.append(
                {
                    "gene": gene,
                    "drug_sensitivity": np.nan,
                    "prism_best_compound": "",
                    "prism_best_broad_id": "",
                    "prism_best_moa": "",
                    "prism_best_target_annotation": "",
                    "prism_best_phase": "",
                    "prism_best_median_auc": np.nan,
                    "prism_n_targeted_compounds": 0,
                }
            )
            continue
        best = matches.sort_values(["prism_compound_sensitivity_score", "prism_n_profiles"], ascending=[False, False]).iloc[0]
        rows.append(
            {
                "gene": gene,
                "drug_sensitivity": float(best["prism_compound_sensitivity_score"]),
                "prism_best_compound": best.get("name", ""),
                "prism_best_broad_id": best.get("broad_id", ""),
                "prism_best_moa": best.get("moa", ""),
                "prism_best_target_annotation": best.get("target", ""),
                "prism_best_phase": best.get("phase", ""),
                "prism_best_median_auc": float(best["prism_median_auc"]),
                "prism_n_targeted_compounds": int(matches["broad_id"].nunique()),
            }
        )
    return pd.DataFrame(rows)
