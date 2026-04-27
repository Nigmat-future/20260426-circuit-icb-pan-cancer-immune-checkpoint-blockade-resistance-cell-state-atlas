from __future__ import annotations

import pandas as pd

from circuit_icb.lincs_reversal import build_lincs_gene_reversal_table, flatten_l1000_response


def test_flatten_l1000_response_extracts_overlap_fields() -> None:
    response = {
        "shareId": "abc",
        "topMeta": [
            {
                "pert_desc": "drug-a",
                "pert_id": "BRD1",
                "cell_id": "A549",
                "score": 0.7,
                "overlap": {"up/dn": ["VIM"], "dn/up": ["CDH1"]},
            }
        ],
    }

    table = flatten_l1000_response("S1", response)

    assert table.loc[0, "lincs_share_id"] == "abc"
    assert table.loc[0, "overlap_up_down"] == "VIM"


def test_build_lincs_gene_reversal_table_bridges_prism_targets() -> None:
    lincs = pd.DataFrame(
        {
            "state": ["S1", "S2"],
            "lincs_rank": [1, 2],
            "pert_desc": ["drug-a", "drug-b"],
            "lincs_reversal_score": [0.7, 0.9],
        }
    )
    prism = pd.DataFrame(
        {
            "broad_id": ["D1", "D2"],
            "name": ["drug-a", "drug-c"],
            "moa": ["m1", "m2"],
            "target": ["A", "B"],
            "phase": ["", ""],
            "auc": [0.2, 0.3],
        }
    )

    table = build_lincs_gene_reversal_table(lincs, prism, ["A", "B"])

    assert table.loc[table["gene"] == "A", "perturbation_reversal"].item() == 0.7
    assert pd.isna(table.loc[table["gene"] == "B", "perturbation_reversal"].item())
