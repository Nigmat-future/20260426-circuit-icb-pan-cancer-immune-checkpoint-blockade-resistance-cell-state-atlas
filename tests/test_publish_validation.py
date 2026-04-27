from __future__ import annotations

import pandas as pd

from circuit_icb.publish_validation import (
    calibration_table,
    circuit_vs_benchmark_comparison,
    direction_stability,
    evaluate_states_with_ci,
    fixed_or_random_effects_meta,
    leave_one_group_logistic,
    random_signature_null,
)


def test_evaluate_states_with_ci_returns_bootstrap_columns() -> None:
    scores = pd.DataFrame(
        {
            "sample": [f"S{i}" for i in range(8)],
            "cohort_id": ["A"] * 4 + ["B"] * 4,
            "response_original": ["R", "R", "NR", "NR", "R", "NR", "R", "NR"],
            "CIRCUIT_STATE_02_score": [0.9, 0.8, 0.2, 0.1, 0.7, 0.3, 0.6, 0.4],
        }
    )

    out = evaluate_states_with_ci(scores, n_bootstrap=20)

    assert {"auroc_ci_low", "auroc_ci_high", "signature_class"} <= set(out.columns)
    meta = fixed_or_random_effects_meta(out)
    assert "random_effect_auroc" in meta.columns


def test_leave_one_group_logistic_and_calibration() -> None:
    scores = pd.DataFrame(
        {
            "sample": [f"S{i}" for i in range(24)],
            "cohort_id": ["A"] * 12 + ["B"] * 12,
            "response_original": ["R", "R", "R", "NR", "NR", "NR"] * 4,
            "CIRCUIT_STATE_02_score": [0.9, 0.8, 0.7, 0.1, 0.2, 0.3] * 4,
            "CIRCUIT_STATE_04_score": [0.8, 0.7, 0.6, 0.2, 0.1, 0.3] * 4,
        }
    )

    stats, predictions = leave_one_group_logistic(scores, group_col="cohort_id", min_test_size=4)
    cal = calibration_table(predictions, n_bins=3)

    assert (stats["status"] == "evaluated").all()
    assert not predictions.empty
    assert not cal.empty


def test_leave_one_group_logistic_skips_small_groups() -> None:
    scores = pd.DataFrame(
        {
            "sample": [f"S{i}" for i in range(16)],
            "cohort_id": ["A"] * 4 + ["B"] * 12,
            "response_original": ["R", "NR", "R", "NR"] + ["R", "R", "R", "NR", "NR", "NR"] * 2,
            "CIRCUIT_STATE_02_score": [0.9, 0.1, 0.8, 0.2] + [0.9, 0.8, 0.7, 0.1, 0.2, 0.3] * 2,
        }
    )

    stats, _ = leave_one_group_logistic(scores, group_col="cohort_id", min_test_size=10)

    a_row = stats.loc[stats["heldout_group"] == "A"]
    assert a_row["status"].iloc[0] == "skipped"


def test_random_signature_null_skips_missing_gene_scores() -> None:
    expression = pd.DataFrame(
        {
            "sample": ["S1", "S2", "S3", "S4"],
            "G1": [1.0, 2.0, None, None],
            "G2": [2.0, 1.0, None, None],
            "G3": [None, None, 3.0, 4.0],
            "G4": [None, None, 4.0, 3.0],
        }
    )
    clinical = pd.DataFrame({"sample": ["S1", "S2", "S3", "S4"], "response_original": ["R", "NR", "R", "NR"]})

    out = random_signature_null(
        expression,
        clinical,
        "CIRCUIT_STATE_02_score",
        signature_size=2,
        excluded_genes=set(),
        n_permutations=5,
    )

    assert {"state", "random_auroc_best_direction"} <= set(out.columns)


def test_direction_stability_detects_flip() -> None:
    per_cohort = pd.DataFrame(
        {
            "cohort": ["A", "B", "C", "A", "B", "C"],
            "state": ["S1"] * 3 + ["S2"] * 3,
            "direction": [
                "higher_in_responders",
                "higher_in_non_responders",
                "higher_in_non_responders",
                "higher_in_responders",
                "higher_in_responders",
                "higher_in_responders",
            ],
            "auroc_best_direction": [0.6] * 6,
        }
    )

    result = direction_stability(per_cohort)
    s1 = result.loc[result["state"] == "S1"].iloc[0]
    s2 = result.loc[result["state"] == "S2"].iloc[0]

    assert not s1["is_stable"]
    assert s2["is_stable"]


def test_circuit_vs_benchmark_comparison() -> None:
    per_cohort = pd.DataFrame(
        {
            "cohort": ["A", "A"],
            "state": ["CIRCUIT_STATE_01_score", "BENCHMARK_IFNG_score"],
            "auroc_best_direction": [0.55, 0.70],
            "signature_class": ["circuit_state", "benchmark"],
            "n": [50, 50],
        }
    )

    result = circuit_vs_benchmark_comparison(per_cohort)

    assert len(result) == 1
    assert result["gap_best"].iloc[0] > 0
