from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table
from circuit_icb.lincs_reversal import build_lincs_gene_reversal_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Bridge LINCS state-reversal perturbations back to candidate target genes through PRISM target annotations.")
    parser.add_argument("--lincs", default="results/tables/target_evidence/lincs_state_reversal.tsv")
    parser.add_argument("--prism", default="data/raw/prism/secondary-screen-dose-response-curve-parameters.csv")
    parser.add_argument("--candidates", default="results/tables/target_evidence/circuit_icb_target_evidence.tsv")
    parser.add_argument("--out", default="results/tables/target_evidence/lincs_gene_reversal.tsv")
    args = parser.parse_args()

    lincs = read_table(args.lincs)
    candidates = read_table(args.candidates)
    prism = pd.read_csv(args.prism, usecols=["broad_id", "auc", "target", "name", "moa", "phase"])
    genes = candidates["gene"].astype(str).drop_duplicates().tolist()
    reversal = build_lincs_gene_reversal_table(lincs, prism, genes)
    write_table(reversal, args.out)
    observed = reversal["perturbation_reversal"].notna().sum()
    print(f"Wrote LINCS gene reversal table to {args.out} ({observed}/{len(reversal)} genes matched)")


if __name__ == "__main__":
    main()
