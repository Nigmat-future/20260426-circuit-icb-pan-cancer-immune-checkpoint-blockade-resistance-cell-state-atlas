from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table
from circuit_icb.signature import load_signatures, rank_state_scores


def score_one(expression: pd.DataFrame, signatures_path: str, min_genes_present: int) -> pd.DataFrame:
    signatures = load_signatures(signatures_path)
    return rank_state_scores(expression, signatures, min_genes_present=min_genes_present)


def main() -> None:
    parser = argparse.ArgumentParser(description="Score CIRCUIT and benchmark signatures in all publish_v1 cohorts.")
    parser.add_argument("--registry", default="data/manifests/publish_v1_bulk_icb_cohorts.tsv")
    parser.add_argument("--state-signatures", default="results/tables/state_annotation/final_state_signatures.yml")
    parser.add_argument("--benchmark-signatures", default="config/benchmark_signatures.yml")
    parser.add_argument("--out-dir", default="results/tables/publish_v1")
    parser.add_argument("--min-genes-present", type=int, default=5)
    args = parser.parse_args()

    registry = read_table(args.registry)
    out_dir = Path(args.out_dir)
    combined_rows = []
    expression_rows = []
    for _, row in registry.loc[registry["status"].eq("included")].iterrows():
        cohort_id = str(row["cohort_id"])
        expression = read_table(row["expression_file"])
        clinical = read_table(row["clinical_file"])
        state_scores = score_one(expression, args.state_signatures, args.min_genes_present)
        benchmark_scores = score_one(expression, args.benchmark_signatures, 3)
        benchmark_scores = benchmark_scores.drop(columns=["sample"])
        scores = pd.concat([state_scores, benchmark_scores], axis=1)
        merged = scores.merge(clinical, on="sample", how="left")
        merged["sample_uid"] = cohort_id + "::" + merged["sample"].astype(str)
        write_table(merged, out_dir / f"{cohort_id}.scores.tsv")
        combined_rows.append(merged)
        expr = expression.copy()
        expr["cohort_id"] = cohort_id
        expr["sample_uid"] = cohort_id + "::" + expr["sample"].astype(str)
        expr["sample"] = expr["sample_uid"]
        expression_rows.append(expr.drop(columns=["sample_uid", "cohort_id"]))
    combined = pd.concat(combined_rows, ignore_index=True) if combined_rows else pd.DataFrame()
    combined_expression = pd.concat(expression_rows, ignore_index=True, sort=False) if expression_rows else pd.DataFrame()
    if not combined.empty:
        combined["sample"] = combined["sample_uid"]
    write_table(combined, out_dir / "publish_v1_all_scores.tsv")
    write_table(combined_expression, out_dir / "publish_v1_all_expression.tsv")
    print(f"Wrote publish_v1 combined scores: {combined.shape}")


if __name__ == "__main__":
    main()
