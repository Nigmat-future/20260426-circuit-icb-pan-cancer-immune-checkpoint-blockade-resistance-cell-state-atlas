from __future__ import annotations

from pathlib import Path

import pandas as pd


def read_table(path: str | Path, **kwargs) -> pd.DataFrame:
    """Read a delimited table using the extension as a sensible default."""

    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", **kwargs)
    if suffix == ".csv":
        return pd.read_csv(path, **kwargs)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, **kwargs)
    raise ValueError(f"Unsupported table extension: {path.suffix}")


def write_table(df: pd.DataFrame, path: str | Path, index: bool = False) -> None:
    """Write a table as TSV or CSV, creating parent directories."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        df.to_csv(path, sep="\t", index=index)
        return
    if suffix == ".csv":
        df.to_csv(path, index=index)
        return
    raise ValueError(f"Unsupported output extension: {path.suffix}")
