from __future__ import annotations

import pandas as pd

from circuit_icb.cptac_protein import build_cptac_protein_evidence


def test_build_cptac_protein_evidence_scales_quantitative_values() -> None:
    protein = pd.DataFrame({"gene": ["A", "B"], "tumor_normal_log2_fc": [2.0, 4.0]})

    out = build_cptac_protein_evidence(protein, ["B", "A", "C"], value_col="tumor_normal_log2_fc")

    assert out.loc[out["gene"] == "B", "protein_evidence"].item() == 1.0
    assert out.loc[out["gene"] == "A", "protein_evidence"].item() == 0.0
    assert pd.isna(out.loc[out["gene"] == "C", "protein_evidence"].item())
