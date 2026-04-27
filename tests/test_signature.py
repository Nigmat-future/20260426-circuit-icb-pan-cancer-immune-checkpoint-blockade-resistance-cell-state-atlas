from __future__ import annotations

import pandas as pd

from circuit_icb.signature import StateSignature, rank_state_scores, ssgsea_state_scores


def test_rank_state_score_orders_expected_sample() -> None:
    expression = pd.DataFrame(
        {
            "sample": ["low", "high"],
            "A": [1.0, 10.0],
            "B": [2.0, 11.0],
            "C": [10.0, 1.0],
            "D": [9.0, 2.0],
        }
    )
    signatures = {"STATE": StateSignature(name="STATE", up=("A", "B"), down=("C", "D"))}

    scores = rank_state_scores(expression, signatures, min_genes_present=2)
    low = scores.loc[scores["sample"] == "low", "STATE_score"].item()
    high = scores.loc[scores["sample"] == "high", "STATE_score"].item()

    assert high > low
    assert scores["STATE_n_genes_present"].iloc[0] == 4


def test_rank_state_score_reports_missing_fraction() -> None:
    expression = pd.DataFrame({"sample": ["S1"], "A": [1.0], "B": [2.0]})
    signatures = {"STATE": StateSignature(name="STATE", up=("A", "B", "C"), down=())}

    scores = rank_state_scores(expression, signatures, min_genes_present=2)

    assert scores["STATE_n_genes_present"].iloc[0] == 2
    assert round(scores["STATE_missing_fraction"].iloc[0], 3) == 0.333


def test_weighted_rank_scoring_respects_weights() -> None:
    expression = pd.DataFrame(
        {
            "sample": ["S1", "S2"],
            "A": [10.0, 1.0],
            "B": [1.0, 10.0],
        }
    )
    signatures = {"STATE": StateSignature(name="STATE", up=("A", "B"), down=())}
    weights = {"STATE": {"A": 10.0, "B": 1.0}}

    scores = rank_state_scores(expression, signatures, min_genes_present=2, weights=weights)
    s1 = scores.loc[scores["sample"] == "S1", "STATE_score"].item()
    s2 = scores.loc[scores["sample"] == "S2", "STATE_score"].item()

    assert s1 > s2


def test_ssgsea_scores_direction() -> None:
    expression = pd.DataFrame(
        {
            "sample": ["low", "high"],
            "A": [1.0, 10.0],
            "B": [2.0, 9.0],
            "C": [8.0, 3.0],
            "D": [9.0, 2.0],
        }
    )
    signatures = {"STATE": StateSignature(name="STATE", up=("A", "B"), down=())}

    scores = ssgsea_state_scores(expression, signatures, min_genes_present=2)
    low = scores.loc[scores["sample"] == "low", "STATE_score"].item()
    high = scores.loc[scores["sample"] == "high", "STATE_score"].item()

    assert high > low
