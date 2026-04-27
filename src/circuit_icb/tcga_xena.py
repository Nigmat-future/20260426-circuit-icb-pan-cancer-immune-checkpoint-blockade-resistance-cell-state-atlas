from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.duration.hazard_regression import PHReg

from circuit_icb.signature import StateSignature


def zscore_signature_scores(
    gene_by_sample: pd.DataFrame,
    signatures: dict[str, StateSignature],
    *,
    min_genes_present: int = 3,
) -> pd.DataFrame:
    """Score signatures on a single platform using gene-wise z-scored expression."""

    expr = gene_by_sample.copy()
    expr.index = expr.index.astype(str).str.upper()
    expr = expr.apply(pd.to_numeric, errors="coerce")
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0).mean()
    means = expr.mean(axis=1, skipna=True)
    stds = expr.std(axis=1, skipna=True).replace(0.0, np.nan)
    zexpr = expr.sub(means, axis=0).div(stds, axis=0)
    records = pd.DataFrame(index=expr.columns)

    for state, signature in signatures.items():
        up_present = [gene for gene in signature.up if gene in zexpr.index]
        down_present = [gene for gene in signature.down if gene in zexpr.index]
        present = len(set(up_present) | set(down_present))
        if present < min_genes_present:
            score = pd.Series(np.nan, index=expr.columns)
        else:
            up_score = zexpr.loc[up_present].mean(axis=0) if up_present else pd.Series(np.nan, index=expr.columns)
            if signature.down:
                down_score = zexpr.loc[down_present].mean(axis=0) if down_present else pd.Series(0.0, index=expr.columns)
                score = up_score - down_score
            else:
                score = up_score
        records[f"{state}_score"] = score
        records[f"{state}_n_genes_present"] = present
        records[f"{state}_n_genes_total"] = signature.total_genes
        records[f"{state}_missing_fraction"] = 1.0 - (present / signature.total_genes if signature.total_genes else 0.0)

    records.insert(0, "sample", records.index)
    return records.reset_index(drop=True)


def summarize_state_landscape(scores: pd.DataFrame, *, group_col: str = "cancer_type") -> pd.DataFrame:
    score_cols = [col for col in scores.columns if col.endswith("_score")]
    rows = []
    for cancer_type, group in scores.groupby(group_col, dropna=False):
        row = {"cancer_type": cancer_type, "n_samples": int(group.shape[0])}
        for col in score_cols:
            values = pd.to_numeric(group[col], errors="coerce")
            row[f"{col}_mean"] = float(values.mean(skipna=True))
            row[f"{col}_median"] = float(values.median(skipna=True))
        rows.append(row)
    return pd.DataFrame(rows)


def cox_by_cancer_state(
    scores: pd.DataFrame,
    *,
    group_col: str = "cancer_type",
    time_col: str = "OS.time",
    event_col: str = "OS",
    min_samples: int = 40,
    min_events: int = 5,
) -> pd.DataFrame:
    rows = []
    score_cols = [col for col in scores.columns if col.endswith("_score") and not col.startswith("BENCHMARK_")]
    for cancer_type, group in scores.groupby(group_col, dropna=False):
        for state_col in score_cols:
            work = group[[time_col, event_col, state_col]].copy()
            work[time_col] = pd.to_numeric(work[time_col], errors="coerce")
            work[event_col] = pd.to_numeric(work[event_col], errors="coerce")
            work[state_col] = pd.to_numeric(work[state_col], errors="coerce")
            work = work.dropna()
            work = work.loc[work[time_col] > 0].copy()
            if work.shape[0] < min_samples or int(work[event_col].sum()) < min_events:
                continue
            std = work[state_col].std()
            if not np.isfinite(std) or std == 0:
                continue
            x = ((work[state_col] - work[state_col].mean()) / std).to_frame("state_z")
            try:
                fit = PHReg(work[time_col], x, status=work[event_col]).fit(disp=False)
                beta = float(fit.params[0])
                se = float(fit.bse[0])
                pvalue = float(fit.pvalues[0])
            except Exception:
                continue
            rows.append(
                {
                    "cancer_type": cancer_type,
                    "state": state_col,
                    "endpoint": time_col.replace(".time", ""),
                    "n": int(work.shape[0]),
                    "n_events": int(work[event_col].sum()),
                    "hazard_ratio_per_sd": math.exp(beta),
                    "ci_low": math.exp(beta - 1.96 * se),
                    "ci_high": math.exp(beta + 1.96 * se),
                    "pvalue": pvalue,
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    order = out["pvalue"].rank(method="first")
    out["fdr_bh"] = (out["pvalue"] * len(out) / order).clip(upper=1.0)
    return out.sort_values(["fdr_bh", "pvalue"]).reset_index(drop=True)


def immune_proxy_correlations(
    scores: pd.DataFrame,
    *,
    group_col: str = "cancer_type",
    proxy_cols: tuple[str, ...] = (
        "BENCHMARK_TCELL_INFLAMED_GEP_score",
        "BENCHMARK_CYTOLYTIC_SCORE_score",
        "BENCHMARK_IFNG_ANTIGEN_PRESENTATION_score",
    ),
) -> pd.DataFrame:
    state_cols = [col for col in scores.columns if col.startswith("CIRCUIT_STATE_") and col.endswith("_score")]
    rows = []
    for cancer_type, group in scores.groupby(group_col, dropna=False):
        for state_col in state_cols:
            for proxy_col in proxy_cols:
                if proxy_col not in group.columns:
                    continue
                work = group[[state_col, proxy_col]].apply(pd.to_numeric, errors="coerce").dropna()
                if work.shape[0] < 10:
                    continue
                rho, pvalue = spearmanr(work[state_col], work[proxy_col])
                rows.append(
                    {
                        "cancer_type": cancer_type,
                        "state": state_col,
                        "immune_proxy": proxy_col,
                        "n": int(work.shape[0]),
                        "spearman_rho": float(rho),
                        "pvalue": float(pvalue),
                    }
                )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    order = out["pvalue"].rank(method="first")
    out["fdr_bh"] = (out["pvalue"] * len(out) / order).clip(upper=1.0)
    return out.sort_values(["fdr_bh", "pvalue"]).reset_index(drop=True)
