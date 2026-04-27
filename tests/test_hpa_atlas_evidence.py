from __future__ import annotations

import pandas as pd

from circuit_icb.hpa_atlas_evidence import build_hpa_atlas_evidence, score_druggability


def test_score_druggability_prioritizes_fda_targets() -> None:
    assert score_druggability("FDA approved drug targets", "Evidence at protein level") == 1.0
    assert score_druggability("Kinases", "Evidence at protein level") == 0.8
    assert score_druggability("", "Evidence at protein level") == 0.25


def test_build_hpa_atlas_evidence_counts_tcga_unfavorable_calls() -> None:
    atlas = pd.DataFrame(
        {
            "gene": ["A"],
            "Protein class": ["Kinases"],
            "Evidence": ["Evidence at protein level"],
            "Cancer prognostics - Cancer 1 (TCGA)": ["prognostic unfavorable (1.0e-4)"],
            "Cancer prognostics - Cancer 2 (TCGA)": ["prognostic favorable (1.0e-3)"],
            "Cancer prognostics - Cancer 3 (validation)": ["prognostic unfavorable (1.0e-4)"],
        }
    )

    table = build_hpa_atlas_evidence(atlas, ["A", "B"])

    a = table.loc[table["gene"] == "A"].iloc[0]
    b = table.loc[table["gene"] == "B"].iloc[0]
    assert a["tcga_unfavorable_count"] == 1
    assert a["tcga_favorable_count"] == 1
    assert a["druggability"] == 0.8
    assert pd.isna(b["tcga_prognosis"])
