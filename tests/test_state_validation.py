from __future__ import annotations

import pandas as pd

from circuit_icb.state_validation import (
    evaluate_univariable_states,
    harmonize_response,
    leave_one_cohort_logistic,
)


def test_harmonize_response_maps_common_labels() -> None:
    assert harmonize_response("PR") == "responder"
    assert harmonize_response("NE") == "non_responder"
    assert harmonize_response("SD") == "non_responder"
    assert harmonize_response("Untreated") == "unknown"


def test_evaluate_univariable_states_reports_auc() -> None:
    df = pd.DataFrame(
        {
            "Combined_outcome": ["Favourable", "Favourable", "Unfavourable", "Unfavourable"],
            "dataset": ["A", "A", "A", "A"],
            "STATE_activity": [0.9, 0.8, 0.2, 0.1],
        }
    )

    result = evaluate_univariable_states(df, outcome_col="Combined_outcome", cohort_col="dataset")

    assert result.loc[result["cohort"] == "pooled", "auroc_high_responder"].iloc[0] == 1.0


def test_leave_one_cohort_logistic_skips_small_cohorts() -> None:
    df = pd.DataFrame(
        {
            "sample": [f"S{i}" for i in range(16)],
            "cohort": ["small"] * 4 + ["big"] * 12,
            "outcome": ["R", "NR", "R", "NR"] + ["R", "R", "R", "NR", "NR", "NR"] * 2,
            "STATE_score": [0.9, 0.1, 0.8, 0.2] + [0.9, 0.8, 0.7, 0.1, 0.2, 0.3] * 2,
        }
    )

    result = leave_one_cohort_logistic(
        df, outcome_col="outcome", cohort_col="cohort", min_test_size=10,
    )

    small_row = result.loc[result["heldout_cohort"] == "small"]
    assert small_row["status"].iloc[0] == "skipped"


def test_leave_one_cohort_logistic_elastic_net() -> None:
    df = pd.DataFrame(
        {
            "sample": [f"S{i}" for i in range(24)],
            "cohort": ["A"] * 12 + ["B"] * 12,
            "outcome": ["R", "R", "R", "NR", "NR", "NR"] * 4,
            "STATE_score": [0.9, 0.8, 0.7, 0.1, 0.2, 0.3] * 4,
        }
    )

    result = leave_one_cohort_logistic(
        df, outcome_col="outcome", cohort_col="cohort",
        min_test_size=4, use_elastic_net=True,
    )

    evaluated = result.loc[result["status"] == "evaluated"]
    assert len(evaluated) >= 1
    assert "model_type" in result.columns
