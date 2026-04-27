from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import read_table, write_table
from circuit_icb.target_priority import apply_publish_v1_target_penalties


def main() -> None:
    parser = argparse.ArgumentParser(description="Create publish_v1 target priority resources focused on core resistance states.")
    parser.add_argument("--annotated", default="results/tables/target_evidence/circuit_icb_target_priority.annotated.tsv")
    parser.add_argument("--out-dir", default="results/tables/publish_v1/target_evidence")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    annotated = apply_publish_v1_target_penalties(read_table(args.annotated))
    core = annotated.loc[annotated["publish_v1_role"].eq("main_text_target")].copy()
    write_table(annotated, out_dir / "publish_v1_target_priority_all.tsv")
    write_table(core, out_dir / "publish_v1_core_resistance_target_priority.tsv")
    print(f"Wrote publish_v1 target resources to {out_dir}")


if __name__ == "__main__":
    main()
