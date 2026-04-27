from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

import pandas as pd


def fetch_json(url: str) -> dict:
    with urllib.request.urlopen(url) as response:
        return json.loads(response.read().decode("utf-8"))


def download(url: str, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, out)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch Zenodo metadata and the small ICB workbook.")
    parser.add_argument("--record", default="14511579")
    parser.add_argument("--out-dir", default="data/raw/zenodo")
    parser.add_argument("--workbook-name", default="v4 ICB datasets.xlsx")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    record = fetch_json(f"https://zenodo.org/api/records/{args.record}")

    files = []
    workbook_url = None
    for item in record.get("files", []):
        key = item.get("key", "")
        link = item.get("links", {}).get("self", "")
        files.append(
            {
                "record": args.record,
                "key": key,
                "size": item.get("size", ""),
                "checksum": item.get("checksum", ""),
                "download_url": link,
            }
        )
        if key == args.workbook_name:
            workbook_url = link

    file_index = out_dir / f"zenodo_{args.record}_files.tsv"
    pd.DataFrame(files).to_csv(file_index, sep="\t", index=False)

    if workbook_url is None:
        raise RuntimeError(f"Could not find workbook {args.workbook_name!r} in Zenodo record {args.record}")

    workbook_out = out_dir / "v4_ICB_datasets.xlsx"
    if not workbook_out.exists():
        download(workbook_url, workbook_out)

    print(f"metadata: {file_index}")
    print(f"workbook: {workbook_out}")


if __name__ == "__main__":
    main()
