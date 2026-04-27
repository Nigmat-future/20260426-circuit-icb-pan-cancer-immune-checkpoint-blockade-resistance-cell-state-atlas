from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table
from circuit_icb.target_priority import load_weight_config, score_targets


def main() -> None:
    parser = argparse.ArgumentParser(description="Score and rank candidate CIRCUIT-ICB targets.")
    parser.add_argument("--evidence", required=True, help="Target evidence table.")
    parser.add_argument("--weights", default="config/target_weights.yml")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    evidence = read_table(args.evidence)
    weights, _ = load_weight_config(args.weights)
    ranked = score_targets(evidence, weights)
    write_table(ranked, args.out)
    print(f"Wrote {ranked.shape[0]} ranked targets to {args.out}")


if __name__ == "__main__":
    main()
