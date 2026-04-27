from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from io import StringIO
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def read_depmap_registry(path: str | Path) -> pd.DataFrame:
    text = Path(path).read_text(encoding="utf-8")
    payload = json.loads(text)
    if not isinstance(payload, str):
        raise ValueError("DepMap registry payload must be a JSON-encoded CSV string")
    return pd.read_csv(StringIO(payload))


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_url(url: str, out: Path, *, expected_md5: str | None = None, force: bool = False) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    start = 0
    if out.exists() and out.stat().st_size > 0 and not force:
        start = out.stat().st_size

    request = urllib.request.Request(url)
    if start:
        request.add_header("Range", f"bytes={start}-")

    with urllib.request.urlopen(request, timeout=600) as response:
        mode = "ab" if start and getattr(response, "status", 200) == 206 else "wb"
        with out.open(mode) as handle:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
    if expected_md5:
        observed = md5sum(out)
        if observed.lower() != expected_md5.lower():
            raise ValueError(f"MD5 mismatch for {out}: expected {expected_md5}, observed {observed}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch selected files from a DepMap public download registry.")
    parser.add_argument("--registry", default="data/raw/depmap/portal_download_files.json")
    parser.add_argument("--release", default="DepMap Public 26Q1")
    parser.add_argument(
        "--filenames",
        nargs="+",
        default=["CRISPRGeneDependency.csv", "Model.csv", "CRISPRInferredCommonEssentials.csv"],
    )
    parser.add_argument("--out-dir", default="data/raw/depmap/26Q1")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    registry = read_depmap_registry(args.registry)
    selected = registry.loc[registry["release"].eq(args.release) & registry["filename"].isin(args.filenames)].copy()
    missing = sorted(set(args.filenames).difference(selected["filename"]))
    if missing:
        raise ValueError(f"Files not found in registry for {args.release}: {missing}")

    out_dir = Path(args.out_dir)
    for _, row in selected.iterrows():
        out = out_dir / row["filename"]
        if out.exists() and out.stat().st_size > 0 and not args.force:
            existing_md5 = md5sum(out)
            if existing_md5.lower() == str(row["md5_hash"]).lower():
                print(f"Using existing {out}")
                continue
            print(f"Resuming or replacing incomplete {out}")
        print(f"Downloading {row['filename']} -> {out}")
        download_url(str(row["url"]), out, expected_md5=str(row["md5_hash"]), force=args.force)
        print(f"Wrote {out} ({out.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
