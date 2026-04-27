from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table
from circuit_icb.state_validation import evaluate_univariable_states, leave_one_cohort_logistic


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate CIRCUIT-ICB state scores against response labels.")
    parser.add_argument("--scores", required=True)
    parser.add_argument("--clinical", default=None)
    parser.add_argument("--merge-on", default="sample")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--outcome-col", default="Combined_outcome")
    parser.add_argument("--cohort-col", default="dataset")
    args = parser.parse_args()

    scores = read_table(args.scores)
    if args.clinical:
        clinical = read_table(args.clinical)
        scores = scores.merge(clinical, on=args.merge_on, how="left", suffixes=("", "_clinical"))
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    univariable = evaluate_univariable_states(scores, outcome_col=args.outcome_col, cohort_col=args.cohort_col)
    write_table(univariable, f"{out_prefix}.univariable.tsv")

    loco = leave_one_cohort_logistic(scores, outcome_col=args.outcome_col, cohort_col=args.cohort_col)
    write_table(loco, f"{out_prefix}.leave_one_cohort.tsv")

    print(f"Wrote validation tables with prefix {out_prefix}")


if __name__ == "__main__":
    main()
