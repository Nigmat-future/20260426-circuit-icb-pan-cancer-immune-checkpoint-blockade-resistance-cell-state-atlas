from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import write_table
from circuit_icb.nmf_states import (
    derive_down_genes,
    detect_lineage_states,
    filter_lineage_states,
    fit_nmf_states,
    load_malignant_pseudobulk,
    scan_nmf_k,
    write_signature_yaml,
)


def parse_k_range(value: str) -> list[int]:
    if "-" in value:
        start, end = value.split("-", 1)
        return list(range(int(start), int(end) + 1))
    return [int(item.strip()) for item in value.split(",") if item.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract NMF malignant-cell states from malignant-only pseudobulk matrices.")
    parser.add_argument("--input-dir", default="data/raw/zenodo/malignant_pseudobulk")
    parser.add_argument("--out-dir", default="results/tables/nmf_states")
    parser.add_argument("--k", type=int, default=7)
    parser.add_argument("--scan-k", default="3-10")
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--marker-top-n", type=int, default=20)
    parser.add_argument("--specificity-ratio", type=float, default=2.0)
    parser.add_argument("--random-state", type=int, default=7)
    parser.add_argument("--max-iter", type=int, default=500)
    parser.add_argument("--scan-max-iter", type=int, default=300)
    parser.add_argument("--scan-n-runs", type=int, default=5)
    parser.add_argument("--n-down-genes", type=int, default=10)
    parser.add_argument("--lineage-threshold", type=float, default=0.3)
    parser.add_argument("--skip-lineage-filter", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    expression, metadata = load_malignant_pseudobulk(args.input_dir)
    write_table(metadata, out_dir / "malignant_pseudobulk_metadata.tsv")

    k_metrics = scan_nmf_k(
        expression,
        k_values=parse_k_range(args.scan_k),
        n_top_genes=args.n_top_genes,
        random_state=args.random_state,
        max_iter=args.scan_max_iter,
        n_runs=args.scan_n_runs,
    )
    write_table(k_metrics, out_dir / "nmf_k_scan_metrics.tsv")
    print(f"K-scan stability: {k_metrics[['n_components', 'stability', 'reconstruction_error']].to_string(index=False)}")

    result = fit_nmf_states(
        expression,
        metadata,
        n_components=args.k,
        n_top_genes=args.n_top_genes,
        marker_top_n=args.marker_top_n,
        random_state=args.random_state,
        max_iter=args.max_iter,
        specificity_ratio=args.specificity_ratio,
    )

    lineage_info = detect_lineage_states(result.markers, overlap_threshold=args.lineage_threshold)
    write_table(lineage_info, out_dir / "nmf_state_lineage_info.tsv")
    print(f"Lineage detection:\n{lineage_info[['state', 'lineage_class', 'lineage_overlap_fraction']].to_string(index=False)}")

    if not args.skip_lineage_filter:
        n_before = len(result.loadings)
        result = filter_lineage_states(result, overlap_threshold=args.lineage_threshold)
        n_after = len(result.loadings)
        if n_after < n_before:
            print(f"Filtered {n_before - n_after} lineage-specific states, {n_after} remain")

    down_genes: dict[str, list[str]] = {}
    for state in result.loadings.index:
        act_col = f"{state}_activity"
        if act_col in result.activity.columns:
            up_genes = set(result.markers.loc[result.markers["state"] == state, "gene"].str.upper())
            down = derive_down_genes(
                expression, result.activity.set_index("sample"),
                act_col, n_down=args.n_down_genes, exclude_genes=up_genes,
            )
            if down:
                down_genes[state] = down
                print(f"  {state} down genes: {down[:5]}...")

    write_table(result.activity, out_dir / "nmf_state_activity.tsv")
    result.loadings.to_csv(out_dir / "nmf_state_gene_loadings.tsv", sep="\t")
    write_table(result.markers, out_dir / "nmf_state_markers.tsv")
    write_table(result.metrics, out_dir / "nmf_fit_metrics.tsv")
    write_table(
        result.markers.groupby("state").head(15).reset_index(drop=True),
        out_dir / "nmf_state_top15_markers.tsv",
    )
    write_signature_yaml(
        result.markers, out_dir / "nmf_state_signatures.yml",
        top_n=args.marker_top_n, down_genes=down_genes,
    )
    print(f"Wrote NMF state outputs to {out_dir}")


if __name__ == "__main__":
    main()
