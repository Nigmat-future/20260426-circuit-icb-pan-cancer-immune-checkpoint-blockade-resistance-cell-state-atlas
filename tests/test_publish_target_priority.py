from __future__ import annotations

import pandas as pd

from circuit_icb.target_priority import apply_publish_v1_target_penalties


def test_publish_v1_target_penalties_downweight_common_essential_safety_and_lineage() -> None:
    priority = pd.DataFrame(
        {
            "gene": ["A", "B", "C"],
            "state_relevance_class": ["core_resistance_state", "lineage_context_state", "core_resistance_state"],
            "conservative_priority_score": [1.0, 1.0, 1.0],
            "depmap_common_essential": [True, False, False],
            "safety": [0.10, 0.10, 0.01],
        }
    )

    out = apply_publish_v1_target_penalties(priority).set_index("gene")

    assert out.loc["A", "publish_v1_adjustment_factor"] == 0.75
    assert out.loc["B", "publish_v1_adjustment_factor"] == 0.50
    assert out.loc["C", "publish_v1_adjustment_factor"] == 0.75
    assert out.loc["A", "publish_v1_role"] == "main_text_target"
