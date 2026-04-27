from __future__ import annotations

from pathlib import Path

import pandas as pd


DISCOVERY_COLUMNS = [
    "source_order",
    "pmid",
    "authors",
    "technology",
    "cancer_type",
    "cancer_subtype",
    "malignant_cells",
    "non_malignant_cells",
    "patients_malignant",
    "patients_other",
    "pre_or_untreated_malignant",
    "post_malignant",
    "pre_post_same_patient",
    "primary_or_metastatic",
    "original_data_location",
    "license_or_dua",
    "doi",
]


FILE_COLUMNS = [
    "object_group",
    "object_type",
    "description",
    "dataset",
    "file_name",
    "location",
]


def _clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def build_discovery_manifest(workbook: str | Path) -> pd.DataFrame:
    """Create a normalized discovery-study manifest from the Zenodo workbook."""

    raw = pd.read_excel(workbook, sheet_name="Metadata for datasets", header=1)
    raw["source_order"] = pd.to_numeric(raw["#"], errors="coerce")
    df = raw.loc[raw["source_order"].notna()].copy()

    # Rows without a PMID are subtype rows belonging to the previous study.
    parent_row = df["PMID"].notna().copy()
    carry_cols = ["PMID", "Authors", "Technology", "Cancer type", "Location of original data", "License/Data Use Agreement", "DOI's"]
    for col in carry_cols:
        last = None
        values = []
        for is_parent, (_, row) in zip(parent_row, df.iterrows()):
            current = row[col]
            if is_parent:
                last = current
                values.append(current)
            else:
                values.append(last if pd.isna(current) else current)
        df[col] = values

    out = pd.DataFrame(
        {
            "source_order": df["source_order"].astype(int),
            "pmid": df["PMID"].apply(_clean_text),
            "authors": df["Authors"].apply(_clean_text),
            "technology": df["Technology"].apply(_clean_text),
            "cancer_type": df["Cancer type"].apply(_clean_text),
            "cancer_subtype": df["Cancer subtype"].apply(_clean_text),
            "malignant_cells": pd.to_numeric(df["# of malignant cells"], errors="coerce").astype("Int64"),
            "non_malignant_cells": pd.to_numeric(df["# of non-malignant cells"], errors="coerce").astype("Int64"),
            "patients_malignant": pd.to_numeric(df["# of patients or samples (mal)"], errors="coerce"),
            "patients_other": pd.to_numeric(df["# of patients or samples (other)"], errors="coerce"),
            "pre_or_untreated_malignant": df["Pre/UT (mal)"].apply(_clean_text),
            "post_malignant": df["Post (mal)"].apply(_clean_text),
            "pre_post_same_patient": df["Pre/Post (same p)"].apply(_clean_text),
            "primary_or_metastatic": df["Primary/Meta"].apply(_clean_text),
            "original_data_location": df["Location of original data"].apply(_clean_text),
            "license_or_dua": df["License/Data Use Agreement"].apply(_clean_text),
            "doi": df["DOI's"].apply(_clean_text),
        }
    )
    return out[DISCOVERY_COLUMNS]


def build_file_manifest(workbook: str | Path) -> pd.DataFrame:
    """Create a normalized file-location manifest from the Zenodo workbook."""

    raw = pd.read_excel(workbook, sheet_name="Data location", header=1)
    raw = raw.dropna(how="all").copy()
    raw = raw.rename(
        columns={
            "Unnamed: 1": "object_group",
            "Objects": "object_type",
            "Description": "description",
            "Datasets": "dataset",
            "File name": "file_name",
            "Location": "location",
        }
    )
    keep = [col for col in FILE_COLUMNS if col in raw.columns]
    out = raw[keep].copy()
    for col in ["object_group", "object_type", "description", "location"]:
        if col in out.columns:
            out[col] = out[col].ffill()
    for col in FILE_COLUMNS:
        if col not in out.columns:
            out[col] = ""
    out = out[FILE_COLUMNS]
    for col in FILE_COLUMNS:
        out[col] = out[col].apply(_clean_text)
    out = out.loc[out["file_name"].ne("") | out["dataset"].ne("")]
    return out


def build_treatment_manifest(workbook: str | Path) -> pd.DataFrame:
    """Create a treatment metadata table from the Zenodo workbook."""

    raw = pd.read_excel(workbook, sheet_name="Treatment_information")
    raw = raw.dropna(how="all").copy()
    raw = raw.rename(
        columns={
            "PMID": "pmid",
            "Sample": "sample",
            "Sex": "sex",
            "Treatment or Mode of Action": "treatment_or_mode",
        }
    )
    for col in ["pmid", "sample", "sex", "treatment_or_mode"]:
        if col not in raw.columns:
            raw[col] = ""
        raw[col] = raw[col].ffill() if col == "pmid" else raw[col]
        raw[col] = raw[col].apply(_clean_text)
    return raw[["pmid", "sample", "sex", "treatment_or_mode"]]


def write_manifest_set(workbook: str | Path, out_dir: str | Path) -> dict[str, Path]:
    """Write all normalized manifests and return their paths."""

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "discovery": out_dir / "icb_sc_discovery_datasets.tsv",
        "files": out_dir / "icb_sc_discovery_files.tsv",
        "treatment": out_dir / "icb_sc_treatment_information.tsv",
    }
    build_discovery_manifest(workbook).to_csv(outputs["discovery"], sep="\t", index=False)
    build_file_manifest(workbook).to_csv(outputs["files"], sep="\t", index=False)
    build_treatment_manifest(workbook).to_csv(outputs["treatment"], sep="\t", index=False)
    return outputs
