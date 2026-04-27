from __future__ import annotations

import pandas as pd

from circuit_icb.depmap_dependency import build_depmap_dependency_table, clean_depmap_gene_column


def test_clean_depmap_gene_column_removes_entrez_suffix() -> None:
    assert clean_depmap_gene_column("BRAF (673)") == "BRAF"
    assert clean_depmap_gene_column("MIR21") == "MIR21"


def test_build_depmap_dependency_table_scores_lineage_signal() -> None:
    matrix = pd.DataFrame(
        {
            "ModelID": ["M1", "M2", "M3"],
            "A (1)": [0.8, 0.7, 0.1],
            "B (2)": [0.1, 0.2, 0.1],
        }
    )
    model = pd.DataFrame({"ModelID": ["M1", "M2", "M3"], "OncotreeLineage": ["skin", "skin", "lung"]})

    table = build_depmap_dependency_table(matrix, model, ["A", "B", "C"], common_essential_genes={"A"})

    a = table.loc[table["gene"] == "A"].iloc[0]
    b = table.loc[table["gene"] == "B"].iloc[0]
    c = table.loc[table["gene"] == "C"].iloc[0]
    assert a["dependency"] > b["dependency"]
    assert a["depmap_top_lineage"] == "skin"
    assert bool(a["depmap_common_essential"]) is True
    assert c["depmap_models_evaluated"] == 0


def test_build_depmap_dependency_table_renames_unnamed_model_column() -> None:
    matrix = pd.DataFrame(
        {
            "Unnamed: 0": ["M1", "M2"],
            "A (1)": [0.9, 0.1],
        }
    )
    model = pd.DataFrame({"ModelID": ["M1", "M2"], "OncotreeLineage": ["skin", "lung"]})

    table = build_depmap_dependency_table(matrix, model, ["A"])

    assert table.loc[0, "depmap_top_lineage"] == "skin"
