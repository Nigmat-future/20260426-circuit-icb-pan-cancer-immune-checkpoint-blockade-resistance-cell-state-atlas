from __future__ import annotations

import pandas as pd

from circuit_icb.hpa_safety import build_hpa_safety_table, normal_expression_safety_score


def test_normal_expression_safety_penalizes_high_and_broad_expression() -> None:
    score = normal_expression_safety_score(
        pd.Series([1.0, 1000.0]),
        pd.Series([0.1, 1.0]),
        max_ntpm_cap=100.0,
    )

    assert score.iloc[0] > score.iloc[1]
    assert 0 <= score.iloc[0] <= 1
    assert 0 <= score.iloc[1] <= 1


def test_build_hpa_safety_table_returns_requested_gene_order() -> None:
    hpa = pd.DataFrame(
        {
            "gene": ["A", "A", "B", "B"],
            "tissue": ["heart", "lung", "heart", "lung"],
            "nTPM": [0.2, 0.4, 50.0, 80.0],
        }
    )

    table = build_hpa_safety_table(hpa, ["B", "A", "C"], expression_threshold=1.0)

    assert table["gene"].tolist() == ["B", "A", "C"]
    assert table.loc[table["gene"] == "A", "safety"].item() > table.loc[table["gene"] == "B", "safety"].item()
    assert pd.isna(table.loc[table["gene"] == "C", "safety"].item())
