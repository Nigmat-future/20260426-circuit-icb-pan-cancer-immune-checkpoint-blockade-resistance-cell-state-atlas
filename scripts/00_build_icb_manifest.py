from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.manifest import write_manifest_set


def main() -> None:
    parser = argparse.ArgumentParser(description="Build normalized ICB scRNA-seq manifests from the Zenodo workbook.")
    parser.add_argument("--workbook", default="data/raw/zenodo/v4_ICB_datasets.xlsx")
    parser.add_argument("--out-dir", default="data/manifests")
    args = parser.parse_args()

    outputs = write_manifest_set(args.workbook, args.out_dir)
    for label, path in outputs.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
