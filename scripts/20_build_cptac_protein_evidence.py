from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.cptac_protein import build_cptac_protein_evidence
from circuit_icb.io import read_table, write_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Adapt a CPTAC/PDC protein evidence table into CIRCUIT-ICB target-priority format.")
    parser.add_argument("--protein-table", required=True, help="TSV/CSV with gene and protein evidence/abundance column.")
    parser.add_argument("--candidates", default="results/tables/target_evidence/circuit_icb_target_evidence.tsv")
    parser.add_argument("--gene-col", default="gene")
    parser.add_argument("--value-col", default="protein_evidence")
    parser.add_argument("--out", default="results/tables/target_evidence/cptac_protein_evidence.tsv")
    args = parser.parse_args()

    candidates = read_table(args.candidates)
    protein = read_table(args.protein_table)
    genes = candidates["gene"].astype(str).drop_duplicates().tolist()
    evidence = build_cptac_protein_evidence(protein, genes, gene_col=args.gene_col, value_col=args.value_col)
    write_table(evidence, args.out)
    print(f"Wrote CPTAC protein evidence to {args.out} ({evidence['protein_evidence'].notna().sum()}/{len(evidence)} genes matched)")


if __name__ == "__main__":
    main()
