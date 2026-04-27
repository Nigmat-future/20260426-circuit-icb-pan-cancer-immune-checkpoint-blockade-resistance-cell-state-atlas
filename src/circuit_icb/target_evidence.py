from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from circuit_icb.state_validation import harmonize_response


EVIDENCE_COLUMNS = [
    "gene",
    "state_specificity",
    "clinical_association",
    "protein_evidence",
    "dependency",
    "drug_sensitivity",
    "perturbation_reversal",
    "safety",
    "druggability",
]


def candidate_genes_from_markers(markers: pd.DataFrame, top_n_per_state: int = 25) -> pd.DataFrame:
    """Create candidate target seed genes from NMF state markers."""

    top = markers.sort_values(["state", "rank"]).groupby("state").head(top_n_per_state).copy()
    max_weight = top.groupby("gene")["weight"].max()
    if max_weight.max() == max_weight.min():
        specificity = max_weight / max_weight.max()
    else:
        specificity = (max_weight - max_weight.min()) / (max_weight.max() - max_weight.min())
    best_state = top.sort_values("weight", ascending=False).drop_duplicates("gene").set_index("gene")["state"]
    evidence = pd.DataFrame(
        {
            "gene": max_weight.index,
            "top_state": max_weight.index.map(best_state),
            "state_marker_weight": max_weight.values,
            "state_specificity": max_weight.index.map(specificity).astype(float),
        }
    )
    return evidence.sort_values("state_specificity", ascending=False).reset_index(drop=True)


def clinical_gene_association(
    expression: pd.DataFrame,
    clinical: pd.DataFrame,
    genes: list[str],
    *,
    outcome_col: str = "response_original",
) -> pd.DataFrame:
    """Score candidate genes by response association in a prepared bulk cohort."""

    merged = expression.merge(clinical[["sample", outcome_col]], on="sample", how="left")
    merged["response_binary"] = merged[outcome_col].map(harmonize_response)
    merged = merged.loc[merged["response_binary"].isin(["responder", "non_responder"])].copy()
    if merged.empty:
        return pd.DataFrame(columns=["gene", "clinical_association", "clinical_direction", "clinical_auroc_best"])
    y = (merged["response_binary"] == "responder").astype(int)

    rows = []
    for gene in genes:
        if gene not in merged.columns or y.nunique() < 2:
            rows.append({"gene": gene, "clinical_association": np.nan, "clinical_direction": "", "clinical_auroc_best": np.nan})
            continue
        values = pd.to_numeric(merged[gene], errors="coerce")
        mask = values.notna()
        if mask.sum() < 4 or y.loc[mask].nunique() < 2:
            rows.append({"gene": gene, "clinical_association": np.nan, "clinical_direction": "", "clinical_auroc_best": np.nan})
            continue
        auc = roc_auc_score(y.loc[mask], values.loc[mask])
        best_auc = max(auc, 1.0 - auc)
        rows.append(
            {
                "gene": gene,
                "clinical_association": max(0.0, (best_auc - 0.5) * 2.0),
                "clinical_direction": "higher_in_responders" if auc >= 0.5 else "higher_in_non_responders",
                "clinical_auroc_best": best_auc,
            }
        )
    return pd.DataFrame(rows)


def pooled_clinical_association(cohort_tables: list[tuple[pd.DataFrame, pd.DataFrame]], genes: list[str]) -> pd.DataFrame:
    """Average clinical association over prepared external bulk cohorts."""

    cohort_results = []
    for idx, (expression, clinical) in enumerate(cohort_tables, start=1):
        result = clinical_gene_association(expression, clinical, genes)
        result["cohort_index"] = idx
        cohort_results.append(result)
    combined = pd.concat(cohort_results, ignore_index=True)
    rows = []
    for gene, group in combined.groupby("gene"):
        valid = group.dropna(subset=["clinical_association"])
        if valid.empty:
            rows.append(
                {
                    "gene": gene,
                    "clinical_association": np.nan,
                    "clinical_auroc_best_mean": np.nan,
                    "clinical_n_cohorts": 0,
                    "clinical_direction": "",
                }
            )
        else:
            direction = valid["clinical_direction"].mode().iloc[0] if not valid["clinical_direction"].mode().empty else ""
            rows.append(
                {
                    "gene": gene,
                    "clinical_association": float(valid["clinical_association"].mean()),
                    "clinical_auroc_best_mean": float(valid["clinical_auroc_best"].mean()),
                    "clinical_n_cohorts": int(valid.shape[0]),
                    "clinical_direction": direction,
                }
            )
    return pd.DataFrame(rows)


def build_target_evidence_matrix(
    markers: pd.DataFrame,
    cohort_tables: list[tuple[pd.DataFrame, pd.DataFrame]],
    *,
    top_n_per_state: int = 25,
) -> pd.DataFrame:
    candidates = candidate_genes_from_markers(markers, top_n_per_state=top_n_per_state)
    clinical = pooled_clinical_association(cohort_tables, candidates["gene"].astype(str).tolist())
    evidence = candidates.merge(clinical, on="gene", how="left")
    for col in ["protein_evidence", "dependency", "drug_sensitivity", "perturbation_reversal", "safety", "druggability"]:
        evidence[col] = np.nan
    return evidence
