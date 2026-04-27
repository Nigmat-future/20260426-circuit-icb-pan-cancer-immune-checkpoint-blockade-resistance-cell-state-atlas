from __future__ import annotations

import pandas as pd

from circuit_icb.prism_sensitivity import build_prism_sensitivity_table, target_tokens


def test_target_tokens_handles_multi_target_annotations() -> None:
    assert {"POLA1", "POLB", "POLD1", "POLE"} <= target_tokens("POLA1, POLB, POLD1, POLE")


def test_build_prism_sensitivity_table_matches_exact_gene_targets() -> None:
    prism = pd.DataFrame(
        {
            "broad_id": ["D1", "D1", "D2", "D3"],
            "name": ["drug1", "drug1", "drug2", "drug3"],
            "moa": ["m1", "m1", "m2", "m3"],
            "target": ["A, B", "A, B", "C", "ABCD"],
            "phase": ["Launched", "Launched", "Phase 2", ""],
            "auc": [0.2, 0.4, 2.0, 0.1],
        }
    )

    out = build_prism_sensitivity_table(prism, ["A", "AB"])

    assert out.loc[out["gene"] == "A", "prism_best_compound"].item() == "drug1"
    assert out.loc[out["gene"] == "A", "drug_sensitivity"].notna().item()
    assert out.loc[out["gene"] == "AB", "prism_n_targeted_compounds"].item() == 0
