from __future__ import annotations

from pathlib import Path

import pandas as pd


def validate_columns(path: str | Path, required_columns: list[str], sep: str = "\t") -> list[str]:
    """Return missing columns for a delimited table."""

    df = pd.read_csv(path, sep=sep, nrows=0)
    return [column for column in required_columns if column not in df.columns]


def summarize_manifest(path: str | Path, sep: str = "\t") -> dict[str, int]:
    """Return a compact table summary for workflow smoke tests."""

    df = pd.read_csv(path, sep=sep)
    return {"rows": int(df.shape[0]), "columns": int(df.shape[1])}
