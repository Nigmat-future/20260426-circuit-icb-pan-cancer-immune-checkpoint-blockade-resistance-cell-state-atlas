from __future__ import annotations

from circuit_icb.benchmark_signatures import BENCHMARK_SIGNATURES, benchmark_summary, signature_gene_sets


def test_benchmark_signatures_include_required_controls() -> None:
    assert "BENCHMARK_HOUSEKEEPING_CONTROL" in BENCHMARK_SIGNATURES
    assert "BENCHMARK_TCELL_INFLAMED_GEP" in BENCHMARK_SIGNATURES
    summary = benchmark_summary()
    assert summary["n_up"].min() > 0
    assert "GAPDH" in signature_gene_sets()["BENCHMARK_HOUSEKEEPING_CONTROL"]
