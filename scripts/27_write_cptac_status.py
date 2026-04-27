from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from circuit_icb.io import write_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Write publish_v1 CPTAC/PDC protein evidence status.")
    parser.add_argument("--out", default="results/tables/publish_v1/cptac/cptac_protein_evidence_status.tsv")
    args = parser.parse_args()

    rows = [
        {
            "resource": "CPTAC/PDC protein evidence",
            "status": "adapter_available_no_matrix_loaded",
            "source_url": "https://pdc.cancer.gov/pdc/",
            "local_adapter": "scripts/20_build_cptac_protein_evidence.py",
            "reason": "No transparent CPTAC/PDC protein matrix was supplied or pre-downloaded for this run.",
        },
        {
            "resource": "Python cptac package",
            "status": "runtime_blocked_windows_build_dependency",
            "source_url": "https://pypi.org/project/cptac/",
            "local_adapter": "",
            "reason": "pip install cptac failed because ncls/sorted_nearest require Microsoft Visual C++ Build Tools in this Windows environment.",
        },
        {
            "resource": "protein_evidence field",
            "status": "left_missing_not_substituted",
            "source_url": "",
            "local_adapter": "results/tables/target_evidence/circuit_icb_target_priority.annotated.tsv",
            "reason": "HPA protein existence is retained in separate HPA fields and is not used as a CPTAC replacement.",
        },
    ]
    write_table(pd.DataFrame(rows), args.out)
    print(f"Wrote CPTAC status to {args.out}")


if __name__ == "__main__":
    main()
