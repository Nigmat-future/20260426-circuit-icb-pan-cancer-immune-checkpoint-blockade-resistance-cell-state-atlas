from __future__ import annotations

from pathlib import Path

import pandas as pd


def main() -> None:
    pairs = [
        (
            "GSE78220",
            Path("results/tables/bulk_icb/GSE78220.state_scores.tsv"),
            Path("data/processed/bulk_icb/GSE78220/clinical.tsv"),
        ),
        (
            "GSE91061",
            Path("results/tables/bulk_icb/GSE91061.state_scores.tsv"),
            Path("data/processed/bulk_icb/GSE91061/clinical.tsv"),
        ),
    ]
    rows = []
    for cohort, scores_path, clinical_path in pairs:
        scores = pd.read_csv(scores_path, sep="\t")
        clinical = pd.read_csv(clinical_path, sep="\t")
        merged = scores.merge(clinical, on="sample", how="left", suffixes=("", "_clinical"))
        merged["cohort_id"] = cohort
        rows.append(merged)

    out = pd.concat(rows, ignore_index=True, sort=False)
    output = Path("results/tables/bulk_icb/locked_external_bulk_state_scores.tsv")
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False)
    print(f"Wrote {out.shape[0]} rows to {output}")


if __name__ == "__main__":
    main()
