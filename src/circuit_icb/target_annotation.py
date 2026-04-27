from __future__ import annotations

import pandas as pd


LINEAGE_CONTEXT_LABELS = {"BREAST_LUMINAL", "LIVER_SECRETORY_COMPLEMENT", "EPITHELIAL_KERATIN"}
CORE_RESISTANCE_LABELS = {"HYPOXIA_STRESS", "EMT_TGFB_EXCLUSION"}
IMMUNE_CONTEXT_LABELS = {"IFNG_ANTIGEN_PRESENTATION"}


def classify_state_label(label: str, decision: str = "keep") -> str:
    label = str(label)
    decision = str(decision)
    if decision == "drop":
        return "dropped_state"
    if label in CORE_RESISTANCE_LABELS:
        return "core_resistance_state"
    if label in IMMUNE_CONTEXT_LABELS:
        return "immune_response_context"
    if label in LINEAGE_CONTEXT_LABELS:
        return "lineage_context_state"
    return "review_state"


def annotate_target_priority(priority: pd.DataFrame, state_annotation: pd.DataFrame) -> pd.DataFrame:
    required_priority = {"gene", "top_state"}
    required_annotation = {"state", "proposed_label", "decision"}
    missing_priority = required_priority.difference(priority.columns)
    missing_annotation = required_annotation.difference(state_annotation.columns)
    if missing_priority:
        raise ValueError(f"Priority table is missing columns: {sorted(missing_priority)}")
    if missing_annotation:
        raise ValueError(f"State annotation table is missing columns: {sorted(missing_annotation)}")

    annotation = state_annotation.rename(columns={"state": "top_state", "decision": "state_decision"})
    keep_cols = [
        "top_state",
        "proposed_label",
        "state_decision",
        "reason",
        "technical_fraction_top_genes",
        "best_pathway",
        "best_fdr",
    ]
    keep_cols = [col for col in keep_cols if col in annotation.columns]
    out = priority.merge(annotation[keep_cols].drop_duplicates("top_state"), on="top_state", how="left")
    out["state_relevance_class"] = [
        classify_state_label(label, decision)
        for label, decision in zip(out["proposed_label"].fillna(""), out["state_decision"].fillna("review"))
    ]
    out["primary_resistance_candidate"] = out["state_relevance_class"].eq("core_resistance_state")
    return out
