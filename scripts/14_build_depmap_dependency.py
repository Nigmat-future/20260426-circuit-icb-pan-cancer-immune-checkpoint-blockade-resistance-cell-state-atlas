from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.depmap_dependency import build_depmap_dependency_table, read_common_essential_genes
from circuit_icb.io import read_table, write_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Build DepMap CRISPR dependency evidence for CIRCUIT-ICB target candidates.")
    parser.add_argument("--dependency-matrix", default="data/raw/depmap/26Q1/CRISPRGeneDependency.csv")
    parser.add_argument("--model", default="data/raw/depmap/26Q1/Model.csv")
    parser.add_argument("--common-essential", default="data/raw/depmap/26Q1/CRISPRInferredCommonEssentials.csv")
    parser.add_argument("--candidates", default="results/tables/target_evidence/circuit_icb_target_evidence.tsv")
    parser.add_argument("--dependency-threshold", type=float, default=0.5)
    parser.add_argument("--out", default="results/tables/target_evidence/depmap_dependency.tsv")
    args = parser.parse_args()

    candidates = read_table(args.candidates)
    if "gene" not in candidates.columns:
        raise ValueError(f"{args.candidates} must contain a gene column")
    genes = candidates["gene"].astype(str).drop_duplicates().tolist()
    matrix = pd.read_csv(args.dependency_matrix)
    model = pd.read_csv(args.model)
    common = read_common_essential_genes(args.common_essential)
    dependency = build_depmap_dependency_table(
        matrix,
        model,
        genes,
        dependency_threshold=args.dependency_threshold,
        common_essential_genes=common,
    )
    write_table(dependency, args.out)
    observed = dependency["dependency"].notna().sum()
    print(f"Wrote DepMap dependency table to {args.out} ({observed}/{len(dependency)} genes matched)")


if __name__ == "__main__":
    main()
