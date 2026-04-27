from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table
from circuit_icb.lincs_reversal import (
    build_state_contrast,
    flatten_l1000_response,
    load_state_signatures,
    query_l1000cds2,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Query L1000CDS2 for state-reversal perturbations.")
    parser.add_argument("--signatures", default="results/tables/state_annotation/final_state_signatures.yml")
    parser.add_argument("--markers", default="results/tables/nmf_states/nmf_state_markers.tsv")
    parser.add_argument("--state-annotation", default="results/tables/state_annotation/state_annotation.tsv")
    parser.add_argument("--core-only", action="store_true", default=True)
    parser.add_argument("--top-n-up", type=int, default=40)
    parser.add_argument("--top-n-down", type=int, default=40)
    parser.add_argument("--top-n-results", type=int, default=50)
    parser.add_argument("--raw-dir", default="data/raw/lincs/l1000cds2")
    parser.add_argument("--out", default="results/tables/target_evidence/lincs_state_reversal.tsv")
    args = parser.parse_args()

    signatures = load_state_signatures(args.signatures)
    markers = read_table(args.markers)
    annotation = read_table(args.state_annotation)
    states = list(signatures)
    if args.core_only:
        core_labels = {"HYPOXIA_STRESS", "EMT_TGFB_EXCLUSION"}
        states = annotation.loc[annotation["proposed_label"].isin(core_labels), "state"].astype(str).tolist()
        states = [state for state in states if state in signatures]

    raw_dir = Path(args.raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    tables = []
    for state in states:
        up, down = build_state_contrast(state, signatures, markers, top_n_up=args.top_n_up, top_n_down=args.top_n_down)
        response = query_l1000cds2(up, down, tag=f"CIRCUIT_ICB_{state}")
        raw_path = raw_dir / f"{state}.json"
        raw_path.write_text(json.dumps({"upGenes": up, "dnGenes": down, "response": response}, indent=2), encoding="utf-8")
        table = flatten_l1000_response(state, response, top_n=args.top_n_results)
        table["query_up_genes"] = ";".join(up)
        table["query_down_genes"] = ";".join(down)
        tables.append(table)
        print(f"Queried {state}: {table.shape[0]} perturbations")

    out = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    write_table(out, args.out)
    print(f"Wrote LINCS state reversal table to {args.out}")


if __name__ == "__main__":
    main()
