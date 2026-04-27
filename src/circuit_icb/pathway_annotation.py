from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def load_pathway_config(path: str | Path) -> tuple[dict[str, set[str]], dict[str, object]]:
    with Path(path).open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}
    pathways = {
        str(name): {str(gene).upper() for gene in genes}
        for name, genes in payload.get("pathways", {}).items()
    }
    technical = payload.get("technical_gene_sets", {})
    return pathways, technical


def technical_fraction(genes: list[str], technical: dict[str, object]) -> float:
    if not genes:
        return 0.0
    upper = [gene.upper() for gene in genes]
    technical_genes: set[str] = set()
    prefixes: list[str] = []
    for spec in technical.values():
        technical_genes.update(str(gene).upper() for gene in spec.get("genes", []) or [])
        prefixes.extend(str(prefix).upper() for prefix in spec.get("prefixes", []) or [])
    flagged = 0
    for gene in upper:
        if gene in technical_genes or any(gene.startswith(prefix) for prefix in prefixes):
            flagged += 1
    return flagged / len(upper)


def enrich_state_markers(
    markers: pd.DataFrame,
    pathways: dict[str, set[str]],
    *,
    top_n: int = 50,
) -> pd.DataFrame:
    universe = set(markers["gene"].astype(str).str.upper())
    rows = []
    for state, group in markers.groupby("state", sort=False):
        genes = group.sort_values("rank").head(top_n)["gene"].astype(str).str.upper().tolist()
        gene_set = set(genes)
        n = len(gene_set)
        m = len(universe)
        for pathway, pathway_genes in pathways.items():
            pathway_in_universe = pathway_genes & universe
            k = len(pathway_in_universe)
            x = len(gene_set & pathway_in_universe)
            pvalue = hypergeom.sf(x - 1, m, k, n) if m and k and n and x else 1.0
            rows.append(
                {
                    "state": state,
                    "pathway": pathway,
                    "overlap": x,
                    "state_genes": n,
                    "pathway_genes_in_universe": k,
                    "pvalue": pvalue,
                    "overlap_genes": ";".join(sorted(gene_set & pathway_in_universe)),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr"] = multipletests(out["pvalue"], method="fdr_bh")[1]
        out["enrichment_score"] = -out["fdr"].clip(lower=1e-300).map(lambda value: __import__("math").log10(value))
    return out


def annotate_states(
    markers: pd.DataFrame,
    enrichment: pd.DataFrame,
    technical: dict[str, object],
    *,
    top_n: int = 50,
    keep_fdr: float = 0.25,
    max_technical_fraction: float = 0.35,
) -> pd.DataFrame:
    rows = []
    for state, group in markers.groupby("state", sort=False):
        genes = group.sort_values("rank").head(top_n)["gene"].astype(str).tolist()
        tech = technical_fraction(genes, technical)
        enr = enrichment.loc[enrichment["state"] == state].sort_values(["fdr", "pvalue", "pathway"])
        best = enr.iloc[0].to_dict() if not enr.empty else {}
        significant = enr.loc[enr["fdr"] <= keep_fdr]
        label = str(best.get("pathway", "UNANNOTATED"))
        decision = "keep"
        reason = "interpretable pathway support"
        if tech > max_technical_fraction and significant.empty:
            decision = "drop"
            reason = "high technical fraction and no pathway support"
        elif significant.empty:
            decision = "review"
            reason = "no pathway passes relaxed FDR threshold"
        rows.append(
            {
                "state": state,
                "proposed_label": label,
                "decision": decision,
                "reason": reason,
                "technical_fraction_top_genes": tech,
                "best_pathway": best.get("pathway", ""),
                "best_overlap": best.get("overlap", 0),
                "best_fdr": best.get("fdr", 1.0),
                "top_markers": ";".join(genes[:15]),
            }
        )
    return pd.DataFrame(rows)


def write_final_signatures(
    markers: pd.DataFrame,
    annotation: pd.DataFrame,
    path: str | Path,
    *,
    top_n: int = 40,
    include_review: bool = True,
) -> None:
    keep_decisions = {"keep", "review"} if include_review else {"keep"}
    selected = annotation.loc[annotation["decision"].isin(keep_decisions)]
    states: dict[str, dict[str, object]] = {}
    for _, row in selected.iterrows():
        state = row["state"]
        genes = (
            markers.loc[markers["state"] == state]
            .sort_values("rank")
            .head(top_n)["gene"]
            .astype(str)
            .tolist()
        )
        states[state] = {
            "description": f"{row['proposed_label']} ({state})",
            "up": genes,
            "down": [],
        }
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        yaml.safe_dump({"states": states}, handle, sort_keys=False)
