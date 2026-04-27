from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table
from circuit_icb.target_priority import load_weight_config, score_targets


def merge_optional(base, path: str | None, column: str):
    if not path:
        return base
    table = read_table(path)
    if "gene" not in table.columns or column not in table.columns:
        raise ValueError(f"{path} must contain columns: gene and {column}")
    table = table.drop_duplicates("gene")
    replace_columns = [col for col in table.columns if col != "gene" and col in base.columns]
    if replace_columns:
        base = base.drop(columns=replace_columns)
    return base.merge(table, on="gene", how="left")


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge optional TCGA/CPTAC/DepMap/PRISM/LINCS evidence into the target matrix.")
    parser.add_argument("--base", default="results/tables/target_evidence/circuit_icb_target_evidence.tsv")
    parser.add_argument("--weights", default="config/target_weights.yml")
    parser.add_argument("--tcga", default=None, help="Optional table with gene,clinical_association or additional columns; reserved for custom use.")
    parser.add_argument("--cptac", default=None, help="Optional table with gene,protein_evidence.")
    parser.add_argument("--depmap", default=None, help="Optional table with gene,dependency.")
    parser.add_argument("--prism", default=None, help="Optional table with gene,drug_sensitivity.")
    parser.add_argument("--lincs", default=None, help="Optional table with gene,perturbation_reversal.")
    parser.add_argument("--safety", default=None, help="Optional table with gene,safety.")
    parser.add_argument("--druggability", default=None, help="Optional table with gene,druggability.")
    parser.add_argument("--out-evidence", default="results/tables/target_evidence/circuit_icb_target_evidence.merged.tsv")
    parser.add_argument("--out-priority", default="results/tables/target_evidence/circuit_icb_target_priority.merged.tsv")
    args = parser.parse_args()

    evidence = read_table(args.base)
    evidence = merge_optional(evidence, args.tcga, "tcga_prognosis")
    evidence = merge_optional(evidence, args.cptac, "protein_evidence")
    evidence = merge_optional(evidence, args.depmap, "dependency")
    evidence = merge_optional(evidence, args.prism, "drug_sensitivity")
    evidence = merge_optional(evidence, args.lincs, "perturbation_reversal")
    evidence = merge_optional(evidence, args.safety, "safety")
    evidence = merge_optional(evidence, args.druggability, "druggability")

    write_table(evidence, args.out_evidence)
    weights, _ = load_weight_config(args.weights)
    priority = score_targets(evidence, weights)
    write_table(priority, args.out_priority)
    print(f"Wrote merged evidence to {args.out_evidence}")
    print(f"Wrote merged priority to {args.out_priority}")


if __name__ == "__main__":
    main()
