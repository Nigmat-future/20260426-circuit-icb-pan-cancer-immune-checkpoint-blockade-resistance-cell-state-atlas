from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table
from circuit_icb.nmf_states import detect_lineage_states
from circuit_icb.pathway_annotation import (
    annotate_states,
    enrich_state_markers,
    load_pathway_config,
    write_final_signatures,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Annotate NMF states with curated resistance/pathway gene sets.")
    parser.add_argument("--markers", default="results/tables/nmf_states/nmf_state_markers.tsv")
    parser.add_argument("--pathways", default="config/pathway_sets.yml")
    parser.add_argument("--out-dir", default="results/tables/state_annotation")
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--signature-top-n", type=int, default=20)
    parser.add_argument("--lineage-threshold", type=float, default=0.3)
    args = parser.parse_args()

    markers = read_table(args.markers)
    pathways, technical = load_pathway_config(args.pathways)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    lineage_info = detect_lineage_states(markers, overlap_threshold=args.lineage_threshold)
    write_table(lineage_info, out_dir / "state_lineage_info.tsv")

    enrichment = enrich_state_markers(markers, pathways, top_n=args.top_n)
    annotation = annotate_states(markers, enrichment, technical, top_n=args.top_n)

    annotation = annotation.merge(
        lineage_info[["state", "lineage_class", "is_lineage_specific"]],
        on="state", how="left",
    )
    annotation.loc[annotation["is_lineage_specific"].fillna(False).astype(bool), "decision"] = "drop"
    annotation.loc[annotation["is_lineage_specific"].fillna(False).astype(bool), "reason"] = "lineage-specific state"

    write_table(enrichment, out_dir / "state_pathway_enrichment.tsv")
    write_table(annotation, out_dir / "state_annotation.tsv")
    write_final_signatures(markers, annotation, out_dir / "final_state_signatures.yml", top_n=args.signature_top_n)

    kept = annotation.loc[annotation["decision"].isin({"keep", "review"})]
    print(f"Kept {len(kept)}/{len(annotation)} states:")
    for _, row in kept.iterrows():
        print(f"  {row['state']}: {row['proposed_label']} (lineage={row.get('lineage_class', 'n/a')})")
    print(f"Wrote state annotation outputs to {out_dir}")


if __name__ == "__main__":
    main()
