from __future__ import annotations

import numpy as np
import pandas as pd

from circuit_icb.nmf_states import (
    detect_lineage_states,
    find_gene_columns,
    markers_from_loadings,
)


def test_find_gene_columns_stops_at_metadata_boundary() -> None:
    columns = ["Unnamed: 0", "A", "B", "sample_id_pre_post", "cell_id", "Outcome"]

    assert find_gene_columns(columns) == ["A", "B"]


def test_markers_from_loadings_returns_ranked_genes() -> None:
    loadings = pd.DataFrame({"A": [0.2], "B": [0.8], "C": [0.1]}, index=["STATE_01"])

    markers = markers_from_loadings(loadings, top_n=2, specificity_ratio=0)

    assert markers["gene"].tolist() == ["B", "A"]
    assert markers["rank"].tolist() == [1, 2]


def test_markers_specificity_filter_removes_shared_genes() -> None:
    loadings = pd.DataFrame(
        {
            "SHARED": [0.5, 0.5],
            "SPECIFIC_A": [0.9, 0.05],
            "SPECIFIC_B": [0.05, 0.9],
        },
        index=["STATE_01", "STATE_02"],
    )

    markers = markers_from_loadings(loadings, top_n=3, specificity_ratio=1.5)

    state1_genes = markers.loc[markers["state"] == "STATE_01", "gene"].tolist()
    state2_genes = markers.loc[markers["state"] == "STATE_02", "gene"].tolist()
    assert "SPECIFIC_A" in state1_genes
    assert "SPECIFIC_B" in state2_genes
    assert "SHARED" not in state1_genes
    assert "SHARED" not in state2_genes


def test_markers_specificity_column_present() -> None:
    loadings = pd.DataFrame({"A": [0.8], "B": [0.2]}, index=["STATE_01"])

    markers = markers_from_loadings(loadings, top_n=2, specificity_ratio=0)

    assert "specificity" in markers.columns


def test_detect_lineage_states_flags_liver() -> None:
    markers = pd.DataFrame({
        "state": ["S1"] * 5 + ["S2"] * 5,
        "gene": ["ALB", "SERPINA1", "APOC1", "VTN", "AMBP",
                 "BRCA1", "TP53", "MYC", "EGFR", "KRAS"],
        "rank": list(range(1, 6)) * 2,
        "weight": [1.0] * 10,
    })

    result = detect_lineage_states(markers, overlap_threshold=0.3)

    s1 = result.loc[result["state"] == "S1"].iloc[0]
    s2 = result.loc[result["state"] == "S2"].iloc[0]
    assert s1["is_lineage_specific"] is True or s1["is_lineage_specific"] == True
    assert s1["lineage_class"] == "liver"
    assert s2["lineage_class"] == "pan_cancer"
