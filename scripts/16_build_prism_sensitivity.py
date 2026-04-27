from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table
from circuit_icb.prism_sensitivity import build_prism_sensitivity_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Build PRISM target-compound sensitivity evidence for CIRCUIT-ICB candidates.")
    parser.add_argument("--prism", default="data/raw/prism/secondary-screen-dose-response-curve-parameters.csv")
    parser.add_argument("--candidates", default="results/tables/target_evidence/circuit_icb_target_evidence.tsv")
    parser.add_argument("--out", default="results/tables/target_evidence/prism_drug_sensitivity.tsv")
    args = parser.parse_args()

    candidates = read_table(args.candidates)
    if "gene" not in candidates.columns:
        raise ValueError(f"{args.candidates} must contain a gene column")
    genes = candidates["gene"].astype(str).drop_duplicates().tolist()
    prism = pd.read_csv(args.prism, usecols=["broad_id", "auc", "target", "name", "moa", "phase"])
    sensitivity = build_prism_sensitivity_table(prism, genes)
    write_table(sensitivity, args.out)
    observed = sensitivity["drug_sensitivity"].notna().sum()
    print(f"Wrote PRISM drug sensitivity table to {args.out} ({observed}/{len(sensitivity)} genes matched)")


if __name__ == "__main__":
    main()
