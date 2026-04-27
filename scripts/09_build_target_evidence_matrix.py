from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table
from circuit_icb.target_evidence import build_target_evidence_matrix
from circuit_icb.target_priority import load_weight_config, score_targets


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a CIRCUIT-ICB target evidence matrix from NMF states and external bulk cohorts.")
    parser.add_argument("--markers", default="results/tables/nmf_states/nmf_state_markers.tsv")
    parser.add_argument("--weights", default="config/target_weights.yml")
    parser.add_argument("--out-evidence", default="results/tables/target_evidence/circuit_icb_target_evidence.tsv")
    parser.add_argument("--out-priority", default="results/tables/target_evidence/circuit_icb_target_priority.tsv")
    parser.add_argument("--top-n-per-state", type=int, default=25)
    args = parser.parse_args()

    cohort_tables = [
        (
            read_table("data/processed/bulk_icb/GSE78220/expression.tsv"),
            read_table("data/processed/bulk_icb/GSE78220/clinical.tsv"),
        ),
        (
            read_table("data/processed/bulk_icb/GSE91061/expression.tsv"),
            read_table("data/processed/bulk_icb/GSE91061/clinical.tsv"),
        ),
    ]
    markers = read_table(args.markers)
    evidence = build_target_evidence_matrix(markers, cohort_tables, top_n_per_state=args.top_n_per_state)
    write_table(evidence, args.out_evidence)

    weights, _ = load_weight_config(args.weights)
    priority = score_targets(evidence, weights)
    write_table(priority, args.out_priority)
    print(f"Wrote evidence matrix to {args.out_evidence}")
    print(f"Wrote target priority table to {args.out_priority}")


if __name__ == "__main__":
    main()
