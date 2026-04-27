from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.hpa_atlas_evidence import build_hpa_atlas_evidence, read_hpa_proteinatlas
from circuit_icb.io import read_table, write_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Build TCGA prognostic and HPA druggability annotations from proteinatlas.tsv.")
    parser.add_argument("--atlas", default="data/raw/hpa/proteinatlas.tsv.zip")
    parser.add_argument("--candidates", default="results/tables/target_evidence/circuit_icb_target_evidence.tsv")
    parser.add_argument("--out", default="results/tables/target_evidence/hpa_atlas_tcga_druggability.tsv")
    args = parser.parse_args()

    candidates = read_table(args.candidates)
    genes = candidates["gene"].astype(str).drop_duplicates().tolist()
    atlas = read_hpa_proteinatlas(args.atlas)
    evidence = build_hpa_atlas_evidence(atlas, genes)
    write_table(evidence, args.out)
    print(f"Wrote HPA atlas evidence to {args.out} ({evidence['tcga_prognosis'].notna().sum()}/{len(evidence)} genes matched)")


if __name__ == "__main__":
    main()
