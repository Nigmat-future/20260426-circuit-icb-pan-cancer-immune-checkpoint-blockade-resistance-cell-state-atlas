from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table
from circuit_icb.target_annotation import annotate_target_priority


def main() -> None:
    parser = argparse.ArgumentParser(description="Annotate target priority rows with state interpretation classes.")
    parser.add_argument("--priority", default="results/tables/target_evidence/circuit_icb_target_priority.merged.tsv")
    parser.add_argument("--state-annotation", default="results/tables/state_annotation/state_annotation.tsv")
    parser.add_argument("--out", default="results/tables/target_evidence/circuit_icb_target_priority.annotated.tsv")
    parser.add_argument("--out-core", default="results/tables/target_evidence/core_resistance_target_priority.tsv")
    args = parser.parse_args()

    priority = read_table(args.priority)
    annotation = read_table(args.state_annotation)
    annotated = annotate_target_priority(priority, annotation)
    write_table(annotated, args.out)
    core = annotated.loc[annotated["primary_resistance_candidate"]].copy()
    write_table(core, args.out_core)
    print(f"Wrote annotated target priority to {args.out}")
    print(f"Wrote core resistance target priority to {args.out_core} ({core.shape[0]} rows)")


if __name__ == "__main__":
    main()
