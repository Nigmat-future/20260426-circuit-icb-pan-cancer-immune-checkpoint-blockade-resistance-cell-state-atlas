from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml


def load_weight_config(path: str | Path) -> tuple[dict[str, float], list[str]]:
    """Load target-priority weights and required evidence columns."""

    with Path(path).open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}
    weights = {str(k): float(v) for k, v in payload.get("weights", {}).items()}
    required = [str(col) for col in payload.get("required_columns", ["gene", *weights.keys()])]
    return weights, required


def score_targets(
    evidence: pd.DataFrame,
    weights: dict[str, float],
    *,
    gene_column: str = "gene",
) -> pd.DataFrame:
    """Compute a transparent weighted target priority score.

    Evidence columns are expected to be directionally aligned, where higher is
    better. Values are clipped to [0, 1]. Missing values contribute zero and are
    reported through completeness columns.
    """

    missing = [col for col in [gene_column, *weights.keys()] if col not in evidence.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    out = evidence.copy()
    total_weight = sum(weights.values())
    weighted_sum = pd.Series(0.0, index=out.index)
    evidence_present = pd.Series(0, index=out.index, dtype="int64")

    for column, weight in weights.items():
        values = pd.to_numeric(out[column], errors="coerce")
        evidence_present += values.notna().astype(int)
        clipped = values.clip(lower=0.0, upper=1.0).fillna(0.0)
        out[f"{column}_weighted"] = clipped * weight
        weighted_sum += out[f"{column}_weighted"]

    out["target_priority_score"] = weighted_sum / total_weight if total_weight else 0.0
    out["n_evidence_present"] = evidence_present
    out["n_evidence_total"] = len(weights)
    out["evidence_completeness"] = evidence_present / len(weights) if weights else 0.0
    out["conservative_priority_score"] = out["target_priority_score"] * out["evidence_completeness"]
    out = out.sort_values(["conservative_priority_score", "target_priority_score"], ascending=False)
    return out.reset_index(drop=True)


def apply_publish_v1_target_penalties(priority: pd.DataFrame) -> pd.DataFrame:
    """Apply manuscript-facing down-weighting for risky or context-specific targets."""

    required = {"state_relevance_class", "conservative_priority_score"}
    missing = required - set(priority.columns)
    if missing:
        raise ValueError(f"Missing required columns for publish_v1 target penalties: {sorted(missing)}")

    out = priority.copy()
    out["publish_v1_role"] = out["state_relevance_class"].map(
        {
            "core_resistance_state": "main_text_target",
            "immune_response_context": "response_context_supplement",
            "lineage_context_state": "lineage_context_supplement",
        }
    ).fillna("review_supplement")
    out["publish_v1_adjustment_factor"] = 1.0
    if "depmap_common_essential" in out.columns:
        out.loc[out["depmap_common_essential"].fillna(False).astype(bool), "publish_v1_adjustment_factor"] *= 0.75
    if "safety" in out.columns:
        safety = pd.to_numeric(out["safety"], errors="coerce").fillna(0.0)
        out.loc[safety < 0.05, "publish_v1_adjustment_factor"] *= 0.75
    out.loc[out["state_relevance_class"].eq("lineage_context_state"), "publish_v1_adjustment_factor"] *= 0.50
    out["publish_v1_priority_score"] = out["conservative_priority_score"] * out["publish_v1_adjustment_factor"]
    return out.sort_values(["publish_v1_role", "publish_v1_priority_score"], ascending=[True, False]).reset_index(drop=True)
