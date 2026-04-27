from __future__ import annotations

from circuit_icb.pathway_annotation import technical_fraction


def test_technical_fraction_flags_prefixes_and_genes() -> None:
    technical = {"x": {"prefixes": ["MT-", "RPL"], "genes": ["MALAT1"]}}
    genes = ["MT-CO1", "RPL3", "MALAT1", "CXCL10"]

    assert technical_fraction(genes, technical) == 0.75
