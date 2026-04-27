from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from circuit_icb.prism_sensitivity import summarize_prism_compounds


L1000CDS2_URL = "https://maayanlab.cloud/L1000CDS2/query"


def load_state_signatures(path: str | Path) -> dict[str, dict[str, list[str]]]:
    with Path(path).open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}
    return payload.get("states", {})


def build_state_contrast(
    state: str,
    signatures: dict[str, dict[str, list[str]]],
    markers: pd.DataFrame,
    *,
    top_n_up: int = 40,
    top_n_down: int = 40,
) -> tuple[list[str], list[str]]:
    up = [str(gene) for gene in signatures[state].get("up", []) if str(gene)]
    up = list(dict.fromkeys(up))[:top_n_up]
    up_set = {gene.upper() for gene in up}
    anti_state = markers.loc[markers["state"].astype(str).ne(state)].sort_values(["rank", "state"])
    down = []
    for gene in anti_state["gene"].astype(str):
        if gene.upper() not in up_set and gene not in down:
            down.append(gene)
        if len(down) >= top_n_down:
            break
    return up, down


def query_l1000cds2(up_genes: list[str], down_genes: list[str], *, tag: str) -> dict[str, object]:
    if len(up_genes) < 3 or len(down_genes) < 3:
        raise ValueError("L1000CDS2 gene-set search requires at least 3 up and 3 down genes")
    payload = {
        "data": {"upGenes": up_genes, "dnGenes": down_genes},
        "config": {
            "aggravate": False,
            "searchMethod": "geneSet",
            "share": False,
            "combination": False,
            "db-version": "latest",
        },
        "metadata": [{"key": "tag", "value": tag}],
    }
    request = urllib.request.Request(
        L1000CDS2_URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def flatten_l1000_response(state: str, response: dict[str, object], *, top_n: int = 50) -> pd.DataFrame:
    rows = []
    share_id = response.get("shareId", "")
    for rank, item in enumerate(response.get("topMeta", [])[:top_n], start=1):
        overlap = item.get("overlap", {}) if isinstance(item, dict) else {}
        rows.append(
            {
                "state": state,
                "lincs_rank": rank,
                "lincs_share_id": share_id,
                "pert_desc": item.get("pert_desc", ""),
                "pert_id": item.get("pert_id", ""),
                "cell_id": item.get("cell_id", ""),
                "pert_dose": item.get("pert_dose", ""),
                "pert_dose_unit": item.get("pert_dose_unit", ""),
                "pert_time": item.get("pert_time", ""),
                "pert_time_unit": item.get("pert_time_unit", ""),
                "lincs_reversal_score": item.get("score", np.nan),
                "lincs_deg_count": item.get("DEGcount", np.nan),
                "pubchem_id": item.get("pubchem_id", ""),
                "overlap_up_down": ";".join(overlap.get("up/dn", []) if isinstance(overlap, dict) else []),
                "overlap_down_up": ";".join(overlap.get("dn/up", []) if isinstance(overlap, dict) else []),
            }
        )
    return pd.DataFrame(rows)


def build_lincs_gene_reversal_table(lincs: pd.DataFrame, prism: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    compounds = summarize_prism_compounds(prism)
    if not compounds.empty:
        compounds["compound_key"] = compounds["name"].astype(str).str.strip().str.lower()
    lincs = lincs.copy()
    lincs["compound_key"] = lincs["pert_desc"].astype(str).str.strip().str.lower()
    lincs["lincs_reversal_score"] = pd.to_numeric(lincs["lincs_reversal_score"], errors="coerce").clip(0.0, 1.0)

    rows = []
    for gene in pd.Series(genes, dtype=str).drop_duplicates():
        gene_key = gene.upper()
        if compounds.empty:
            targeted = compounds
        else:
            targeted = compounds.loc[compounds["target_tokens"].map(lambda tokens: gene_key in tokens)].copy()
        if targeted.empty:
            rows.append(
                {
                    "gene": gene,
                    "perturbation_reversal": np.nan,
                    "lincs_best_perturbation": "",
                    "lincs_best_state": "",
                    "lincs_best_score": np.nan,
                    "lincs_matched_prism_target": "",
                    "lincs_n_matching_compounds": 0,
                }
            )
            continue
        matches = lincs.merge(targeted[["compound_key", "target", "moa"]].drop_duplicates(), on="compound_key", how="inner")
        if matches.empty:
            rows.append(
                {
                    "gene": gene,
                    "perturbation_reversal": np.nan,
                    "lincs_best_perturbation": "",
                    "lincs_best_state": "",
                    "lincs_best_score": np.nan,
                    "lincs_matched_prism_target": "",
                    "lincs_n_matching_compounds": 0,
                }
            )
            continue
        best = matches.sort_values(["lincs_reversal_score", "lincs_rank"], ascending=[False, True]).iloc[0]
        rows.append(
            {
                "gene": gene,
                "perturbation_reversal": float(best["lincs_reversal_score"]),
                "lincs_best_perturbation": best["pert_desc"],
                "lincs_best_state": best["state"],
                "lincs_best_score": float(best["lincs_reversal_score"]),
                "lincs_matched_prism_target": best["target"],
                "lincs_n_matching_compounds": int(matches["compound_key"].nunique()),
            }
        )
    return pd.DataFrame(rows)
