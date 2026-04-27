from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.benchmark_signatures import signature_gene_sets
from circuit_icb.io import read_table, write_table
from circuit_icb.publish_validation import (
    calibration_table,
    circuit_vs_benchmark_comparison,
    direction_stability,
    evaluate_states_with_ci,
    fixed_or_random_effects_meta,
    leave_one_group_logistic,
    random_signature_null,
)


def load_signature_sizes(path: str) -> dict[str, int]:
    with Path(path).open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}
    sizes = {}
    for name, spec in payload.get("states", {}).items():
        genes = {str(g).upper() for g in spec.get("up", [])} | {str(g).upper() for g in spec.get("down", [])}
        sizes[f"{name}_score"] = len(genes)
    return sizes


def main() -> None:
    parser = argparse.ArgumentParser(description="Run publish_v1 validation statistics.")
    parser.add_argument("--scores", default="results/tables/publish_v1/publish_v1_all_scores.tsv")
    parser.add_argument("--expression", default="results/tables/publish_v1/publish_v1_all_expression.tsv")
    parser.add_argument("--state-signatures", default="results/tables/state_annotation/final_state_signatures.yml")
    parser.add_argument("--out-dir", default="results/tables/publish_v1/validation")
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--n-random", type=int, default=1000)
    parser.add_argument("--min-test-size", type=int, default=10)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    scores = read_table(args.scores)
    expression = read_table(args.expression)

    univariable = evaluate_states_with_ci(scores, n_bootstrap=args.n_bootstrap)
    write_table(univariable, out_dir / "publish_v1_univariable_ci.tsv")
    meta = fixed_or_random_effects_meta(univariable)
    write_table(meta, out_dir / "publish_v1_random_effects_meta.tsv")

    stability = direction_stability(univariable)
    write_table(stability, out_dir / "publish_v1_direction_stability.tsv")
    print("Direction stability:")
    for _, row in stability.iterrows():
        flag = "STABLE" if row["is_stable"] else "UNSTABLE"
        print(f"  {row['state']}: {row['majority_direction']} ({row['direction_agreement']:.0%}) [{flag}]")

    comparison = circuit_vs_benchmark_comparison(univariable)
    write_table(comparison, out_dir / "publish_v1_circuit_vs_benchmark.tsv")
    if not comparison.empty:
        print(f"\nCircuit vs Benchmark gap (mean AUROC):")
        for _, row in comparison.iterrows():
            print(f"  {row['cohort']}: circuit={row['circuit_mean_auroc']:.3f} bench={row['benchmark_mean_auroc']:.3f} gap={row['gap_mean']:+.3f}")

    loco, cohort_predictions = leave_one_group_logistic(
        scores, group_col="cohort_id", min_test_size=args.min_test_size,
    )
    write_table(loco, out_dir / "publish_v1_leave_one_cohort.tsv")
    locancer, cancer_predictions = leave_one_group_logistic(
        scores, group_col="cancer_type", min_test_size=args.min_test_size,
    )
    write_table(locancer, out_dir / "publish_v1_leave_one_cancer_type.tsv")
    predictions = pd.concat([cohort_predictions.assign(split_type="cohort"), cancer_predictions.assign(split_type="cancer_type")], ignore_index=True)
    write_table(predictions, out_dir / "publish_v1_logistic_predictions.tsv")
    write_table(calibration_table(cohort_predictions), out_dir / "publish_v1_calibration.tsv")

    print(f"\nLeave-one-cohort-out results:")
    for _, row in loco.iterrows():
        if row["status"] == "evaluated":
            print(f"  {row['heldout_group']}: AUROC={row['auroc']:.3f} AUPRC={row['auprc']:.3f} (n={row['n_test']})")
        else:
            print(f"  {row['heldout_group']}: SKIPPED ({row['reason']})")

    sizes = load_signature_sizes(args.state_signatures)
    excluded = set()
    for genes in signature_gene_sets().values():
        excluded.update(genes)

    circuit_states = [col for col in scores.columns if col.startswith("CIRCUIT_STATE_") and col.endswith("_score")]
    primary_states = [s for s in circuit_states if s in sizes]

    random_tables = []
    for state in primary_states:
        if state not in scores.columns:
            continue
        random_tables.append(
            random_signature_null(
                expression,
                scores[["sample", "response_original"]].copy(),
                state,
                signature_size=sizes[state],
                excluded_genes=excluded,
                n_permutations=args.n_random,
            )
        )
    random_null = pd.concat(random_tables, ignore_index=True) if random_tables else pd.DataFrame()
    write_table(random_null, out_dir / "publish_v1_random_signature_null.tsv")
    if not random_null.empty:
        observed = univariable.loc[univariable["cohort"].eq("pooled"), ["state", "auroc_best_direction"]]
        observed = observed.loc[observed["state"].isin(primary_states)].copy()
        merged = observed.merge(random_null, on="state", how="left")
        rows = []
        for (state, observed_auc), group in merged.groupby(["state", "auroc_best_direction"], dropna=False):
            values = group["random_auroc_best_direction"].dropna()
            rows.append(
                {
                    "state": state,
                    "auroc_best_direction": observed_auc,
                    "n_random": int(values.shape[0]),
                    "permutation_p": float(((values >= observed_auc).sum() + 1) / (values.shape[0] + 1)),
                }
            )
        pvals = pd.DataFrame(rows)
    else:
        pvals = pd.DataFrame(columns=["state", "auroc_best_direction", "n_random", "permutation_p"])
    write_table(pvals, out_dir / "publish_v1_random_signature_pvalues.tsv")

    print(f"\nRandom signature null test:")
    for _, row in pvals.iterrows():
        sig = "***" if row["permutation_p"] < 0.001 else "**" if row["permutation_p"] < 0.01 else "*" if row["permutation_p"] < 0.05 else "ns"
        print(f"  {row['state']}: AUROC={row['auroc_best_direction']:.3f} p={row['permutation_p']:.4f} [{sig}]")

    print(f"\nWrote publish_v1 validation outputs to {out_dir}")


if __name__ == "__main__":
    main()
