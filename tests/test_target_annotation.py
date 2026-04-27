from __future__ import annotations

import pandas as pd

from circuit_icb.target_annotation import annotate_target_priority, classify_state_label


def test_classify_state_label_separates_resistance_from_lineage() -> None:
    assert classify_state_label("HYPOXIA_STRESS") == "core_resistance_state"
    assert classify_state_label("BREAST_LUMINAL") == "lineage_context_state"
    assert classify_state_label("IFNG_ANTIGEN_PRESENTATION") == "immune_response_context"


def test_annotate_target_priority_adds_primary_flag() -> None:
    priority = pd.DataFrame({"gene": ["A", "B"], "top_state": ["S1", "S2"]})
    annotation = pd.DataFrame(
        {
            "state": ["S1", "S2"],
            "proposed_label": ["EMT_TGFB_EXCLUSION", "BREAST_LUMINAL"],
            "decision": ["keep", "keep"],
        }
    )

    out = annotate_target_priority(priority, annotation)

    assert out.loc[out["gene"] == "A", "primary_resistance_candidate"].item() is True
    assert out.loc[out["gene"] == "B", "primary_resistance_candidate"].item() is False
