from __future__ import annotations

import pandas as pd

from circuit_icb.cohort_schema import PUBLISH_V1_CLINICAL_COLUMNS, standardize_publish_clinical, validate_publish_clinical
from circuit_icb.state_validation import harmonize_response


def test_standardize_publish_clinical_adds_required_columns() -> None:
    clinical = pd.DataFrame({"sample": ["S1", "S2"], "response_original": ["R", "NR"], "timing": ["pre", "pre"]})

    out = standardize_publish_clinical(
        clinical,
        cohort_id="C1",
        cancer_type="melanoma",
        therapy_class="anti-PD-1",
        sample_timing_col="timing",
        source_url="https://example.org",
    )

    assert out.columns.tolist() == PUBLISH_V1_CLINICAL_COLUMNS
    assert out["response_binary"].tolist() == ["responder", "non_responder"]
    summary = validate_publish_clinical(out)
    assert summary.loc[0, "n_binary_response"] == 2


def test_objective_response_maps_stable_disease_to_non_responder() -> None:
    assert harmonize_response("CR") == "responder"
    assert harmonize_response("SD") == "non_responder"
    assert harmonize_response("non-responder") == "non_responder"
    assert harmonize_response("stable_or_mixed") == "stable_or_mixed"
