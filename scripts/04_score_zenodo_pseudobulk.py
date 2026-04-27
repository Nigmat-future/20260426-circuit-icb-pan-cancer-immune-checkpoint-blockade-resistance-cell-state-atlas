from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import write_table
from circuit_icb.signature import load_signatures, rank_state_scores


def _dataset_from_file(path: Path) -> str:
    name = path.stem
    return name.removesuffix("_data")


def score_file(path: Path, signatures_path: Path, min_genes_present: int) -> tuple[pd.DataFrame, dict[str, object]]:
    expression = pd.read_csv(path)
    first_col = expression.columns[0]
    dataset = _dataset_from_file(path)
    expression = expression.rename(columns={first_col: "sample"})
    expression["sample"] = expression["sample"].map(lambda value: f"{dataset}__{value}")
    signatures = load_signatures(signatures_path)
    scores = rank_state_scores(expression, signatures, min_genes_present=min_genes_present)
    scores.insert(1, "dataset", dataset)
    profile = {
        "dataset": dataset,
        "file": str(path),
        "rows": expression.shape[0],
        "genes_or_columns": expression.shape[1] - 1,
    }
    return scores, profile


def main() -> None:
    parser = argparse.ArgumentParser(description="Score downloaded Zenodo pseudobulk matrices with state signatures.")
    parser.add_argument("--pseudobulk-dir", default="data/raw/zenodo/pseudobulk")
    parser.add_argument("--signatures", default="config/state_signatures.example.yml")
    parser.add_argument("--out", default="results/tables/zenodo_pseudobulk_state_scores.tsv")
    parser.add_argument("--profile-out", default="data/manifests/zenodo_pseudobulk_profile.tsv")
    parser.add_argument("--min-genes-present", type=int, default=3)
    args = parser.parse_args()

    root = Path(args.pseudobulk_dir)
    files = sorted(root.rglob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files found under {root}")

    all_scores = []
    profiles = []
    for path in files:
        scores, profile = score_file(path, Path(args.signatures), args.min_genes_present)
        all_scores.append(scores)
        profiles.append(profile)

    combined = pd.concat(all_scores, ignore_index=True)
    profile_df = pd.DataFrame(profiles)
    write_table(combined, args.out)
    write_table(profile_df, args.profile_out)
    print(f"Wrote {combined.shape[0]} scored pseudobulk rows to {args.out}")
    print(f"Wrote pseudobulk profile to {args.profile_out}")


if __name__ == "__main__":
    main()
