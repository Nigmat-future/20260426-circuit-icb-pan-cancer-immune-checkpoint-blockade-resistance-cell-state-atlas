from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd


def table_shape(path: Path) -> tuple[int, int]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        if not header:
            return 0, 0
        columns = len(header.rstrip("\n").split("\t"))
        rows = sum(1 for _ in handle)
    return rows, columns


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a compact release index from generated result tables.")
    parser.add_argument("--results-dir", default="results/tables")
    parser.add_argument("--out", default="results/tables/release_index.tsv")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    rows = []
    for path in sorted(results_dir.rglob("*.tsv")):
        if path.name == Path(args.out).name:
            continue
        try:
            n_rows, n_cols = table_shape(path)
            rows.append({"file": str(path), "rows": n_rows, "columns": n_cols})
        except Exception as exc:
            rows.append({"file": str(path), "rows": "", "columns": "", "error": str(exc)})

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"Wrote release index to {out}")


if __name__ == "__main__":
    main()
