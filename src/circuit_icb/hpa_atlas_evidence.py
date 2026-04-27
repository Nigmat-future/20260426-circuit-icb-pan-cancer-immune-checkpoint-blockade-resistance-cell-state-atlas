from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


def read_hpa_proteinatlas(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as zf:
            names = [name for name in zf.namelist() if name.endswith(".tsv")]
            if not names:
                raise ValueError(f"No TSV file found in {path}")
            with zf.open(names[0]) as handle:
                table = pd.read_csv(handle, sep="\t")
    else:
        table = pd.read_csv(path, sep="\t")
    if "Gene" not in table.columns:
        raise ValueError("HPA protein atlas table must contain a Gene column")
    return table.rename(columns={"Gene": "gene"})


def score_druggability(protein_class: object, evidence: object) -> float:
    classes = str(protein_class).lower() if not pd.isna(protein_class) else ""
    evidence_text = str(evidence).lower() if not pd.isna(evidence) else ""
    if "fda approved drug targets" in classes:
        return 1.0
    if any(term in classes for term in ["kinases", "g-protein coupled receptors", "transporters", "enzymes"]):
        return 0.8
    if any(term in classes for term in ["predicted membrane proteins", "cd markers", "plasma proteins"]):
        return 0.6
    if any(term in classes for term in ["secreted proteins", "surface proteins"]):
        return 0.5
    if "protein level" in evidence_text:
        return 0.25
    return 0.0


def score_protein_existence(evidence: object) -> float:
    text = str(evidence).lower() if not pd.isna(evidence) else ""
    if "protein level" in text:
        return 1.0
    if "transcript level" in text:
        return 0.5
    if "uncertain" in text:
        return 0.2
    return np.nan


def summarize_tcga_prognostics(row: pd.Series, prognostic_columns: list[str]) -> dict[str, object]:
    unfavorable = 0
    potential_unfavorable = 0
    favorable = 0
    tested = 0
    adverse_cancers = []
    protective_cancers = []
    for column in prognostic_columns:
        value = row.get(column)
        if pd.isna(value):
            continue
        text = str(value).lower()
        if not text:
            continue
        tested += 1
        cancer = column.replace("Cancer prognostics - ", "").replace(" (TCGA)", "")
        if "unfavorable" in text:
            if "potential prognostic unfavorable" in text:
                potential_unfavorable += 1
            else:
                unfavorable += 1
            adverse_cancers.append(cancer)
        elif "favorable" in text:
            favorable += 1
            protective_cancers.append(cancer)
    weighted_adverse = unfavorable + 0.5 * potential_unfavorable
    score = weighted_adverse / tested if tested else np.nan
    return {
        "tcga_prognosis": score,
        "tcga_unfavorable_count": unfavorable,
        "tcga_potential_unfavorable_count": potential_unfavorable,
        "tcga_favorable_count": favorable,
        "tcga_prognostic_cancers_tested": tested,
        "tcga_adverse_cancers": ";".join(adverse_cancers),
        "tcga_protective_cancers": ";".join(protective_cancers),
    }


def build_hpa_atlas_evidence(table: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    requested = pd.DataFrame({"gene": pd.Series(genes, dtype=str).drop_duplicates()})
    subset = table.loc[table["gene"].astype(str).isin(requested["gene"])].copy()
    prognostic_columns = [
        col for col in table.columns if col.startswith("Cancer prognostics - ") and col.endswith("(TCGA)")
    ]
    rows = []
    for gene in requested["gene"]:
        hit = subset.loc[subset["gene"].astype(str).eq(gene)]
        if hit.empty:
            rows.append(
                {
                    "gene": gene,
                    "tcga_prognosis": np.nan,
                    "druggability": np.nan,
                    "hpa_protein_existence": np.nan,
                    "hpa_protein_class": "",
                    "hpa_evidence": "",
                    "tcga_unfavorable_count": 0,
                    "tcga_potential_unfavorable_count": 0,
                    "tcga_favorable_count": 0,
                    "tcga_prognostic_cancers_tested": 0,
                    "tcga_adverse_cancers": "",
                    "tcga_protective_cancers": "",
                }
            )
            continue
        row = hit.iloc[0]
        tcga = summarize_tcga_prognostics(row, prognostic_columns)
        evidence = row.get("Evidence", "")
        protein_class = row.get("Protein class", "")
        rows.append(
            {
                "gene": gene,
                **tcga,
                "druggability": score_druggability(protein_class, evidence),
                "hpa_protein_existence": score_protein_existence(evidence),
                "hpa_protein_class": protein_class if not pd.isna(protein_class) else "",
                "hpa_evidence": evidence if not pd.isna(evidence) else "",
            }
        )
    return pd.DataFrame(rows)
