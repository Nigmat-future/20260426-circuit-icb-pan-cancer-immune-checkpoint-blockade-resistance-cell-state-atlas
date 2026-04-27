from __future__ import annotations

import pandas as pd

from circuit_icb.target_priority import score_targets


def test_target_priority_ranks_multievidence_target_first() -> None:
    evidence = pd.DataFrame(
        {
            "gene": ["A", "B"],
            "state_specificity": [0.9, 1.0],
            "clinical_association": [0.9, 0.2],
            "dependency": [0.8, 0.1],
        }
    )
    weights = {"state_specificity": 1.0, "clinical_association": 1.0, "dependency": 1.0}

    ranked = score_targets(evidence, weights)

    assert ranked.loc[0, "gene"] == "A"
    assert "conservative_priority_score" in ranked.columns


def test_target_priority_penalizes_missing_evidence() -> None:
    evidence = pd.DataFrame(
        {
            "gene": ["A", "B"],
            "state_specificity": [1.0, 0.8],
            "clinical_association": [None, 0.8],
        }
    )
    weights = {"state_specificity": 1.0, "clinical_association": 1.0}

    ranked = score_targets(evidence, weights)

    assert ranked.loc[ranked["gene"] == "A", "evidence_completeness"].item() == 0.5
    assert ranked.loc[ranked["gene"] == "B", "evidence_completeness"].item() == 1.0
