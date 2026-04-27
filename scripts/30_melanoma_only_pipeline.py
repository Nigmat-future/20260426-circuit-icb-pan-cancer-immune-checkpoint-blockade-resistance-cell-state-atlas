"""Melanoma-only pipeline: NMF discovery → signature scoring → validation."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.benchmark_signatures import signature_gene_sets
from circuit_icb.io import read_table, write_table
from circuit_icb.nmf_states import (
    _prepare_nmf_matrix,
    derive_down_genes,
    detect_lineage_states,
    filter_lineage_states,
    fit_nmf_states,
    load_malignant_pseudobulk,
    markers_from_loadings,
    scan_nmf_k,
    select_variable_genes,
    write_signature_yaml,
)
from circuit_icb.publish_validation import (
    calibration_table,
    circuit_vs_benchmark_comparison,
    direction_stability,
    evaluate_states_with_ci,
    fixed_or_random_effects_meta,
    leave_one_group_logistic,
    random_signature_null,
    repeated_cv_summary,
    repeated_stratified_cv,
)
from circuit_icb.signature import load_signatures, rank_state_scores

MELANOMA_DATASETS = {"Jerby_Arnon_Mel", "Pozniak_Mel", "Tirosh_Mel"}
MELANOMA_COHORTS = {"GSE78220", "GSE91061"}
GENERALIZATION_COHORTS = {"GSE126044", "GSE135222", "GSE176307"}


def filter_melanoma_pseudobulk(
    expression: pd.DataFrame, metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mask = metadata["dataset"].isin(MELANOMA_DATASETS)
    mel_samples = metadata.loc[mask, "sample"].tolist()
    return (
        expression.loc[expression.index.isin(mel_samples)],
        metadata.loc[mask].reset_index(drop=True),
    )


def main() -> None:
    out_dir = ROOT / "results" / "tables" / "melanoma_only"
    out_dir.mkdir(parents=True, exist_ok=True)
    val_dir = out_dir / "validation"
    val_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Load & filter melanoma pseudobulk ──
    print("Loading pseudobulk data...")
    expression, metadata = load_malignant_pseudobulk(
        ROOT / "data" / "raw" / "zenodo" / "malignant_pseudobulk",
    )
    mel_expr, mel_meta = filter_melanoma_pseudobulk(expression, metadata)
    print(f"Melanoma pseudobulk: {mel_expr.shape[0]} samples, {mel_expr.shape[1]} genes")
    print(f"Datasets: {mel_meta['dataset'].value_counts().to_dict()}")
    write_table(mel_meta, out_dir / "melanoma_pseudobulk_metadata.tsv")

    # ── 2. Scan k with stability ──
    print("\nScanning k=2..8 with stability analysis...")
    k_metrics = scan_nmf_k(
        mel_expr, k_values=list(range(2, 9)),
        n_top_genes=3000, random_state=7, max_iter=500, n_runs=10,
    )
    write_table(k_metrics, out_dir / "melanoma_k_scan.tsv")
    print(k_metrics[["n_components", "stability", "reconstruction_error"]].to_string(index=False))

    best_k_row = k_metrics.loc[k_metrics["stability"] >= 0.95].sort_values(
        ["stability", "reconstruction_error"], ascending=[False, True],
    )
    best_k = int(best_k_row.iloc[0]["n_components"]) if not best_k_row.empty else 4
    print(f"\nSelected k={best_k} (highest stability with lowest error)")

    # ── 3. Fit NMF ──
    print(f"\nFitting NMF with k={best_k}...")
    result = fit_nmf_states(
        mel_expr, mel_meta,
        n_components=best_k, n_top_genes=3000,
        marker_top_n=20, random_state=7, max_iter=1000,
        specificity_ratio=2.0,
    )

    lineage_info = detect_lineage_states(result.markers, overlap_threshold=0.3)
    write_table(lineage_info, out_dir / "melanoma_lineage_info.tsv")
    print(f"Lineage detection:\n{lineage_info[['state', 'lineage_class', 'lineage_overlap_fraction']].to_string(index=False)}")

    n_before = len(result.loadings)
    result = filter_lineage_states(result, overlap_threshold=0.3)
    n_after = len(result.loadings)
    if n_after < n_before:
        print(f"Filtered {n_before - n_after} lineage states, {n_after} remain")

    # ── 4. Derive down genes ──
    down_genes: dict[str, list[str]] = {}
    for state in result.loadings.index:
        act_col = f"{state}_activity"
        if act_col in result.activity.columns:
            up_genes = set(result.markers.loc[result.markers["state"] == state, "gene"].str.upper())
            down = derive_down_genes(
                mel_expr, result.activity.set_index("sample"),
                act_col, n_down=10, exclude_genes=up_genes,
            )
            if down:
                down_genes[state] = down

    # ── 5. Write outputs ──
    write_table(result.activity, out_dir / "melanoma_nmf_activity.tsv")
    result.loadings.to_csv(out_dir / "melanoma_nmf_loadings.tsv", sep="\t")
    write_table(result.markers, out_dir / "melanoma_nmf_markers.tsv")
    write_table(result.metrics, out_dir / "melanoma_nmf_metrics.tsv")
    write_table(
        result.markers.groupby("state").head(15).reset_index(drop=True),
        out_dir / "melanoma_top15_markers.tsv",
    )
    sig_path = out_dir / "melanoma_state_signatures.yml"
    write_signature_yaml(result.markers, sig_path, top_n=20, down_genes=down_genes)

    print(f"\nTop markers per state:")
    for state in result.loadings.index:
        genes = result.markers.loc[result.markers["state"] == state].head(5)["gene"].tolist()
        down = down_genes.get(state, [])[:3]
        print(f"  {state} UP: {genes}  DOWN: {down}")

    # ── 6. Score melanoma validation cohorts ──
    print("\n── Scoring melanoma validation cohorts ──")
    registry = read_table(ROOT / "data" / "manifests" / "publish_v1_bulk_icb_cohorts.tsv")
    bench_path = ROOT / "config" / "benchmark_signatures.yml"
    signatures = load_signatures(sig_path)
    bench_sigs = load_signatures(bench_path)

    combined_rows = []
    expression_rows = []
    for _, row in registry.loc[registry["status"].eq("included")].iterrows():
        cohort_id = str(row["cohort_id"])
        if cohort_id not in (MELANOMA_COHORTS | GENERALIZATION_COHORTS):
            continue
        expr = read_table(row["expression_file"])
        clin = read_table(row["clinical_file"])
        state_scores = rank_state_scores(expr, signatures, min_genes_present=3)
        bench_scores = rank_state_scores(expr, bench_sigs, min_genes_present=3)
        bench_scores = bench_scores.drop(columns=["sample"])
        scores = pd.concat([state_scores, bench_scores], axis=1)
        merged = scores.merge(clin, on="sample", how="left")
        merged["sample_uid"] = cohort_id + "::" + merged["sample"].astype(str)
        split = "melanoma" if cohort_id in MELANOMA_COHORTS else "generalization"
        merged["validation_split"] = split
        combined_rows.append(merged)
        e = expr.copy()
        e["sample"] = cohort_id + "::" + e["sample"].astype(str)
        expression_rows.append(e)
        print(f"  {cohort_id} ({split}): {len(merged)} samples scored")

    combined = pd.concat(combined_rows, ignore_index=True)
    combined["sample"] = combined["sample_uid"]
    combined_expr = pd.concat(expression_rows, ignore_index=True, sort=False)
    write_table(combined, out_dir / "melanoma_all_scores.tsv")
    write_table(combined_expr, out_dir / "melanoma_all_expression.tsv")

    # ── 7. Validate on melanoma cohorts ──
    mel_scores = combined.loc[combined["validation_split"] == "melanoma"].copy()
    print(f"\n── Melanoma validation ({len(mel_scores)} samples) ──")

    univariable = evaluate_states_with_ci(mel_scores, n_bootstrap=1000)
    write_table(univariable, val_dir / "melanoma_univariable_ci.tsv")

    pooled = univariable.loc[univariable["cohort"] == "pooled"].copy()
    print("\nPooled univariable AUROC:")
    for _, r in pooled.iterrows():
        cls = r["signature_class"]
        print(f"  {r['state']}: {r['auroc_best_direction']:.3f} ({r['direction']}) [{cls}]")

    meta = fixed_or_random_effects_meta(univariable)
    write_table(meta, val_dir / "melanoma_random_effects_meta.tsv")

    stab = direction_stability(univariable)
    write_table(stab, val_dir / "melanoma_direction_stability.tsv")

    comp = circuit_vs_benchmark_comparison(univariable)
    write_table(comp, val_dir / "melanoma_circuit_vs_benchmark.tsv")

    # Repeated stratified CV (multi-state logistic)
    cv_folds, cv_preds = repeated_stratified_cv(
        mel_scores, n_splits=5, n_repeats=50,
    )
    write_table(cv_folds, val_dir / "melanoma_cv_folds.tsv")
    write_table(cv_preds, val_dir / "melanoma_cv_predictions.tsv")
    write_table(calibration_table(cv_preds.rename(columns={"mean_predicted_probability": "predicted_probability"})), val_dir / "melanoma_calibration.tsv")
    cv_summary = repeated_cv_summary(cv_folds)
    write_table(pd.DataFrame([cv_summary]), val_dir / "melanoma_cv_summary.tsv")

    print(f"\nRepeated 5-fold CV (50 repeats, multi-state logistic):")
    print(f"  AUROC: {cv_summary['cv_auroc_mean']:.3f} +/- {cv_summary['cv_auroc_std']:.3f} (95% CI: {cv_summary['cv_auroc_ci_low']:.3f}-{cv_summary['cv_auroc_ci_high']:.3f})")
    print(f"  AUPRC: {cv_summary['cv_auprc_mean']:.3f} +/- {cv_summary['cv_auprc_std']:.3f}")
    print(f"  Brier: {cv_summary['cv_brier_mean']:.3f} +/- {cv_summary['cv_brier_std']:.3f}")

    # Also run CV with benchmark features for comparison
    cv_folds_bench, cv_preds_bench = repeated_stratified_cv(
        mel_scores, feature_prefixes=("BENCHMARK_",), n_splits=5, n_repeats=50,
    )
    cv_bench_summary = repeated_cv_summary(cv_folds_bench)
    if cv_bench_summary:
        write_table(pd.DataFrame([cv_bench_summary]), val_dir / "melanoma_cv_benchmark_summary.tsv")
        print(f"\nBenchmark signatures CV (for comparison):")
        print(f"  AUROC: {cv_bench_summary['cv_auroc_mean']:.3f} +/- {cv_bench_summary['cv_auroc_std']:.3f}")

    # Combined: CIRCUIT + BENCHMARK
    cv_folds_all, cv_preds_all = repeated_stratified_cv(
        mel_scores, feature_prefixes=("CIRCUIT_STATE_", "BENCHMARK_"), n_splits=5, n_repeats=50,
    )
    cv_all_summary = repeated_cv_summary(cv_folds_all)
    if cv_all_summary:
        write_table(pd.DataFrame([cv_all_summary]), val_dir / "melanoma_cv_combined_summary.tsv")
        print(f"\nCombined (CIRCUIT + BENCHMARK) CV:")
        print(f"  AUROC: {cv_all_summary['cv_auroc_mean']:.3f} +/- {cv_all_summary['cv_auroc_std']:.3f}")

    # LOCO kept for reference
    loco, cohort_preds = leave_one_group_logistic(
        mel_scores, group_col="cohort_id", min_test_size=10,
    )
    write_table(loco, val_dir / "melanoma_leave_one_cohort.tsv")
    print(f"\nLeave-one-cohort-out (reference):")
    for _, r in loco.iterrows():
        if r["status"] == "evaluated":
            print(f"  {r['heldout_group']}: AUROC={r['auroc']:.3f} (n={r['n_test']})")
        else:
            print(f"  {r['heldout_group']}: SKIPPED ({r['reason']})")

    # Random signature null
    sizes = {}
    with sig_path.open("r") as f:
        sig_yaml = yaml.safe_load(f) or {}
    for name, spec in sig_yaml.get("states", {}).items():
        genes = {str(g).upper() for g in spec.get("up", [])} | {str(g).upper() for g in spec.get("down", [])}
        sizes[f"{name}_score"] = len(genes)

    excluded = set()
    for genes in signature_gene_sets().values():
        excluded.update(genes)

    circuit_states = [c for c in mel_scores.columns if c.startswith("CIRCUIT_STATE_") and c.endswith("_score")]
    random_tables = []
    for state in circuit_states:
        if state not in sizes:
            continue
        mel_expr_val = combined_expr.loc[combined_expr["sample"].isin(mel_scores["sample"])]
        random_tables.append(random_signature_null(
            mel_expr_val, mel_scores[["sample", "response_original"]].copy(),
            state, signature_size=sizes[state], excluded_genes=excluded,
            n_permutations=1000,
        ))
    random_null = pd.concat(random_tables, ignore_index=True) if random_tables else pd.DataFrame()
    write_table(random_null, val_dir / "melanoma_random_signature_null.tsv")

    if not random_null.empty:
        observed = pooled.loc[pooled["state"].isin(circuit_states), ["state", "auroc_best_direction"]]
        merged_r = observed.merge(random_null, on="state", how="left")
        pval_rows = []
        for (state, obs_auc), grp in merged_r.groupby(["state", "auroc_best_direction"], dropna=False):
            vals = grp["random_auroc_best_direction"].dropna()
            pval_rows.append({
                "state": state, "auroc_best_direction": obs_auc,
                "n_random": int(vals.shape[0]),
                "permutation_p": float(((vals >= obs_auc).sum() + 1) / (vals.shape[0] + 1)),
            })
        pvals = pd.DataFrame(pval_rows)
    else:
        pvals = pd.DataFrame(columns=["state", "auroc_best_direction", "n_random", "permutation_p"])
    write_table(pvals, val_dir / "melanoma_random_signature_pvalues.tsv")

    print("\nRandom signature null test (melanoma):")
    for _, r in pvals.iterrows():
        sig = "***" if r["permutation_p"] < 0.001 else "**" if r["permutation_p"] < 0.01 else "*" if r["permutation_p"] < 0.05 else "ns"
        print(f"  {r['state']}: AUROC={r['auroc_best_direction']:.3f} p={r['permutation_p']:.4f} [{sig}]")

    # ── 8. Generalization test on NSCLC + urothelial ──
    gen_scores = combined.loc[combined["validation_split"] == "generalization"].copy()
    if not gen_scores.empty:
        print(f"\n── Generalization test ({len(gen_scores)} samples) ──")
        gen_uni = evaluate_states_with_ci(gen_scores, n_bootstrap=1000)
        write_table(gen_uni, val_dir / "generalization_univariable_ci.tsv")
        gen_pooled = gen_uni.loc[gen_uni["cohort"] == "pooled"]
        print("Pooled univariable AUROC (generalization):")
        for _, r in gen_pooled.iterrows():
            if r["state"].startswith("CIRCUIT_STATE_"):
                print(f"  {r['state']}: {r['auroc_best_direction']:.3f} ({r['direction']})")

    print(f"\nAll outputs written to {out_dir}")


if __name__ == "__main__":
    main()
