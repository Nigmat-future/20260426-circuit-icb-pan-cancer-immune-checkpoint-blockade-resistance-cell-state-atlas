from __future__ import annotations

import pandas as pd

from circuit_icb.signature import StateSignature
from circuit_icb.tcga_xena import summarize_state_landscape, zscore_signature_scores


def test_zscore_signature_scores_returns_sample_level_scores() -> None:
    expression = pd.DataFrame(
        {
            "S1": [1.0, 4.0, 2.0],
            "S2": [2.0, 3.0, 2.0],
            "S3": [3.0, 2.0, 2.0],
            "S4": [4.0, 1.0, 2.0],
        },
        index=["G1", "G2", "G3"],
    )
    signatures = {"STATE_A": StateSignature("STATE_A", up=("G1",), down=("G2",))}

    scores = zscore_signature_scores(expression, signatures, min_genes_present=2)

    assert scores.shape[0] == 4
    assert "STATE_A_score" in scores.columns
    assert scores.loc[scores["sample"].eq("S4"), "STATE_A_score"].iloc[0] > scores.loc[scores["sample"].eq("S1"), "STATE_A_score"].iloc[0]


def test_summarize_state_landscape_groups_by_cancer() -> None:
    scores = pd.DataFrame(
        {
            "sample": ["S1", "S2", "S3"],
            "cancer_type": ["A", "A", "B"],
            "STATE_A_score": [1.0, 3.0, 5.0],
        }
    )

    out = summarize_state_landscape(scores)

    assert out.loc[out["cancer_type"].eq("A"), "STATE_A_score_mean"].iloc[0] == 2.0
