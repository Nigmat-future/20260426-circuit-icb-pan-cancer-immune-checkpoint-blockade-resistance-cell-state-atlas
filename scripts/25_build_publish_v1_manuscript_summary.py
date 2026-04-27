from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table


PRIMARY_STATES = {
    "CIRCUIT_STATE_01_score",
    "CIRCUIT_STATE_02_score",
    "CIRCUIT_STATE_04_score",
    "CIRCUIT_STATE_06_score",
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Build manuscript-ready publish_v1 summary tables.")
    parser.add_argument("--publish-dir", default="results/tables/publish_v1")
    parser.add_argument("--out", default="results/tables/publish_v1/manuscript_ready_summary.tsv")
    args = parser.parse_args()

    publish_dir = Path(args.publish_dir)
    cohort_summary = read_table(publish_dir / "cohort_summary.tsv")
    validation = read_table(publish_dir / "validation" / "publish_v1_univariable_ci.tsv")
    meta_path = publish_dir / "validation" / "publish_v1_random_effects_meta.tsv"
    meta = read_table(meta_path) if meta_path.exists() else pd.DataFrame()
    targets = read_table(publish_dir / "target_evidence" / "publish_v1_core_resistance_target_priority.tsv")
    curation = read_table(publish_dir / "cohort_curation_log.tsv")
    random_path = publish_dir / "validation" / "publish_v1_random_signature_pvalues.tsv"
    random_p = read_table(random_path) if random_path.exists() else pd.DataFrame()
    loco_path = publish_dir / "validation" / "publish_v1_leave_one_cohort.tsv"
    loco = read_table(loco_path) if loco_path.exists() else pd.DataFrame()
    tcga_path = publish_dir / "tcga" / "tcga_state_scores.tsv"
    tcga = read_table(tcga_path) if tcga_path.exists() else pd.DataFrame()
    tcga_survival_path = publish_dir / "tcga" / "tcga_state_survival_os.tsv"
    tcga_survival = read_table(tcga_survival_path) if tcga_survival_path.exists() else pd.DataFrame()
    cptac_path = publish_dir / "cptac" / "cptac_protein_evidence_status.tsv"
    cptac_status = read_table(cptac_path) if cptac_path.exists() else pd.DataFrame()

    pooled = validation.loc[validation["cohort"].eq("pooled")].copy() if not validation.empty else pd.DataFrame()
    pooled_primary = pooled.loc[pooled["state"].isin(PRIMARY_STATES)].copy() if not pooled.empty else pd.DataFrame()
    best_primary = pooled_primary.sort_values("auroc_best_direction", ascending=False).iloc[0] if not pooled_primary.empty else None
    acceptance = "met_minimum_bulk_validation_gate" if int(cohort_summary["n_binary_response"].sum()) >= 150 and int((curation["status"] == "included").sum()) >= 4 else "below_minimum_bulk_validation_gate"
    mean_loco_auc = float(pd.to_numeric(loco.loc[loco["status"].eq("evaluated"), "auroc"], errors="coerce").mean()) if not loco.empty else float("nan")
    classifier_claim = "association_resource_not_classifier" if not mean_loco_auc or pd.isna(mean_loco_auc) or mean_loco_auc < 0.60 else "classifier_performance_requires_external_review"

    rows = [
        {"section": "cohorts", "metric": "included_cohorts", "value": int((curation["status"] == "included").sum())},
        {"section": "cohorts", "metric": "binary_response_samples", "value": int(cohort_summary["n_binary_response"].sum()) if not cohort_summary.empty else 0},
        {"section": "cohorts", "metric": "minimum_gate_status", "value": acceptance},
        {"section": "validation", "metric": "best_primary_pooled_state", "value": best_primary["state"] if best_primary is not None else ""},
        {"section": "validation", "metric": "best_primary_pooled_auroc", "value": float(best_primary["auroc_best_direction"]) if best_primary is not None else ""},
        {"section": "validation", "metric": "meta_analyzed_states", "value": int(meta.shape[0]) if not meta.empty else 0},
        {"section": "validation", "metric": "mean_leave_one_cohort_auroc", "value": mean_loco_auc},
        {"section": "validation", "metric": "random_signature_tests", "value": int(random_p["n_random"].sum()) if not random_p.empty else 0},
        {"section": "tcga", "metric": "primary_tumors_scored", "value": int(tcga.shape[0]) if not tcga.empty else 0},
        {"section": "tcga", "metric": "os_cox_tests", "value": int(tcga_survival.shape[0]) if not tcga_survival.empty else 0},
        {"section": "cptac", "metric": "protein_evidence_status", "value": ";".join(cptac_status["status"].astype(str).tolist()) if not cptac_status.empty else "not_run"},
        {"section": "targets", "metric": "core_resistance_targets", "value": int(targets.shape[0]) if not targets.empty else 0},
        {"section": "claim_boundary", "metric": "classifier_claim", "value": classifier_claim},
        {"section": "claim_boundary", "metric": "primary_claim", "value": "Retrospective public-data state-resource and target-prioritization framework; not a clinical-use predictor."},
    ]
    write_table(pd.DataFrame(rows), args.out)
    print(f"Wrote manuscript-ready summary to {args.out}")


if __name__ == "__main__":
    main()
