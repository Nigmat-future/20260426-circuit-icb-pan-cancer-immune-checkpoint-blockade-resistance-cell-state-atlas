from __future__ import annotations

from pathlib import Path

import pandas as pd

from circuit_icb.state_validation import harmonize_response


PUBLISH_V1_CLINICAL_COLUMNS = [
    "sample",
    "cohort_id",
    "cancer_type",
    "therapy_class",
    "sample_timing",
    "response_original",
    "response_binary",
    "dc_benefit",
    "pfs_time",
    "pfs_event",
    "os_time",
    "os_event",
    "source_url",
    "processing_note",
]


def durable_clinical_benefit(value: object) -> str:
    response = harmonize_response(value)
    if response == "responder":
        return "benefit"
    if response == "non_responder":
        return "no_benefit"
    return "unknown"


def standardize_publish_clinical(
    clinical: pd.DataFrame,
    *,
    cohort_id: str,
    cancer_type: str,
    therapy_class: str,
    source_url: str,
    sample_timing_col: str | None = None,
    processing_note: str = "",
) -> pd.DataFrame:
    out = clinical.copy()
    if "sample" not in out.columns:
        raise ValueError("Clinical table must contain a sample column")
    if "response_original" not in out.columns:
        out["response_original"] = ""
    out["cohort_id"] = cohort_id
    out["cancer_type"] = cancer_type
    out["therapy_class"] = therapy_class
    out["sample_timing"] = out[sample_timing_col] if sample_timing_col and sample_timing_col in out.columns else "pre-treatment"
    out["response_binary"] = out["response_original"].map(harmonize_response)
    out["dc_benefit"] = out["response_original"].map(durable_clinical_benefit)
    out["source_url"] = source_url
    out["processing_note"] = processing_note
    for column in ["pfs_time", "pfs_event", "os_time", "os_event"]:
        if column not in out.columns:
            out[column] = pd.NA
    return out[PUBLISH_V1_CLINICAL_COLUMNS].copy()


def validate_publish_clinical(clinical: pd.DataFrame) -> pd.DataFrame:
    missing = [column for column in PUBLISH_V1_CLINICAL_COLUMNS if column not in clinical.columns]
    if missing:
        raise ValueError(f"Clinical table is missing publish_v1 columns: {missing}")
    rows = []
    for cohort_id, group in clinical.groupby("cohort_id", dropna=False):
        binary = group.loc[group["response_binary"].isin(["responder", "non_responder"])]
        rows.append(
            {
                "cohort_id": cohort_id,
                "n_samples": int(group.shape[0]),
                "n_binary_response": int(binary.shape[0]),
                "n_responders": int((binary["response_binary"] == "responder").sum()),
                "n_non_responders": int((binary["response_binary"] == "non_responder").sum()),
                "cancer_type": ";".join(sorted(group["cancer_type"].dropna().astype(str).unique())),
                "therapy_class": ";".join(sorted(group["therapy_class"].dropna().astype(str).unique())),
                "sample_timing": ";".join(sorted(group["sample_timing"].dropna().astype(str).unique())),
            }
        )
    return pd.DataFrame(rows)


def build_cohort_curation_row(
    *,
    cohort_id: str,
    accession: str,
    status: str,
    reason: str,
    expression_file: str = "",
    clinical_file: str = "",
    source_url: str = "",
    notes: str = "",
) -> dict[str, object]:
    return {
        "cohort_id": cohort_id,
        "accession": accession,
        "status": status,
        "reason": reason,
        "expression_file": expression_file,
        "clinical_file": clinical_file,
        "source_url": source_url,
        "notes": notes,
    }


def existing_file_or_empty(path: str | Path) -> str:
    path = Path(path)
    return str(path) if path.exists() else ""
