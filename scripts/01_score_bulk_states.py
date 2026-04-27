from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table
from circuit_icb.signature import load_signatures, rank_state_scores


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute bulk-compatible malignant-cell state scores.")
    parser.add_argument("--expression", required=True, help="Sample-by-gene expression table with a sample column.")
    parser.add_argument("--signatures", default="config/state_signatures.example.yml")
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-genes-present", type=int, default=3)
    args = parser.parse_args()

    expression = read_table(args.expression)
    signatures = load_signatures(args.signatures)
    scores = rank_state_scores(expression, signatures, min_genes_present=args.min_genes_present)
    write_table(scores, args.out)
    print(f"Wrote {scores.shape[0]} samples x {scores.shape[1]} columns to {args.out}")


if __name__ == "__main__":
    main()
