from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.hpa_safety import build_hpa_safety_table, read_hpa_tissue_consensus
from circuit_icb.io import read_table, write_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Build HPA normal tissue expression safety evidence for candidate targets.")
    parser.add_argument("--hpa", default="data/raw/hpa/rna_tissue_consensus.tsv.zip")
    parser.add_argument("--candidates", default="results/tables/target_evidence/circuit_icb_target_evidence.tsv")
    parser.add_argument("--gene-col", default="gene")
    parser.add_argument("--expression-threshold", type=float, default=1.0)
    parser.add_argument("--max-ntpm-cap", type=float, default=100.0)
    parser.add_argument("--out", default="results/tables/target_evidence/hpa_safety.tsv")
    args = parser.parse_args()

    candidates = read_table(args.candidates)
    if args.gene_col not in candidates.columns:
        raise ValueError(f"{args.candidates} must contain column {args.gene_col}")
    genes = candidates[args.gene_col].astype(str).dropna().drop_duplicates().tolist()
    hpa = read_hpa_tissue_consensus(args.hpa)
    safety = build_hpa_safety_table(
        hpa,
        genes,
        expression_threshold=args.expression_threshold,
        max_ntpm_cap=args.max_ntpm_cap,
    )
    write_table(safety, args.out)
    observed = safety["safety"].notna().sum()
    print(f"Wrote HPA safety table to {args.out} ({observed}/{len(safety)} genes matched)")


if __name__ == "__main__":
    main()
