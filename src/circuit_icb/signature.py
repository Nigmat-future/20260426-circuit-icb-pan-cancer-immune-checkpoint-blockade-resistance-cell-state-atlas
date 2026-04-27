from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


@dataclass(frozen=True)
class StateSignature:
    name: str
    up: tuple[str, ...]
    down: tuple[str, ...] = ()
    description: str = ""

    @property
    def total_genes(self) -> int:
        return len(set(self.up) | set(self.down))


def _normalize_gene(gene: str) -> str:
    return str(gene).strip().upper()


def load_signatures(path: str | Path) -> dict[str, StateSignature]:
    """Load state signatures from a YAML file."""

    with Path(path).open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}
    states = payload.get("states", payload)
    signatures: dict[str, StateSignature] = {}
    for name, spec in states.items():
        up = tuple(_normalize_gene(g) for g in spec.get("up", []) if str(g).strip())
        down = tuple(_normalize_gene(g) for g in spec.get("down", []) if str(g).strip())
        signatures[str(name)] = StateSignature(
            name=str(name),
            up=up,
            down=down,
            description=str(spec.get("description", "")),
        )
    return signatures


def prepare_expression(expression: pd.DataFrame, sample_column: str = "sample") -> pd.DataFrame:
    """Return a sample-by-gene matrix with normalized gene symbols."""

    expr = expression.copy()
    if sample_column in expr.columns:
        expr = expr.set_index(sample_column)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr.columns = [_normalize_gene(col) for col in expr.columns]
    if expr.columns.duplicated().any():
        expr = expr.T.groupby(level=0).mean().T
    return expr


def rank_state_scores(
    expression: pd.DataFrame,
    signatures: dict[str, StateSignature],
    *,
    min_genes_present: int = 3,
    weights: dict[str, dict[str, float]] | None = None,
) -> pd.DataFrame:
    """Compute rank-based state scores for a sample-by-gene expression matrix.

    The score is the mean within-sample percentile rank of up genes minus the
    mean within-sample percentile rank of down genes. If no down genes are
    supplied, the score is simply the mean percentile rank of up genes.

    If *weights* is provided (state_name -> {gene: weight}), each gene's
    percentile rank is multiplied by its weight before averaging.
    """

    expr = prepare_expression(expression)
    needed_genes = sorted(
        {
            gene
            for signature in signatures.values()
            for gene in (*signature.up, *signature.down)
            if gene in expr.columns
        }
    )
    rank_lookup = _percentile_ranks_for_genes(expr, needed_genes)
    records = pd.DataFrame(index=expr.index)

    for state, signature in signatures.items():
        up_present = [gene for gene in signature.up if gene in rank_lookup.columns]
        down_present = [gene for gene in signature.down if gene in rank_lookup.columns]
        total = signature.total_genes
        present = len(set(up_present) | set(down_present))

        w = weights.get(state) if weights else None
        up_score = _weighted_rank_mean(rank_lookup, up_present, w, expr.index)
        if signature.down:
            down_score = _weighted_rank_mean(rank_lookup, down_present, w, expr.index)
            if down_score is None:
                down_score = pd.Series(0.0, index=expr.index)
            score = up_score - down_score if up_score is not None else pd.Series(float("nan"), index=expr.index)
        else:
            score = up_score if up_score is not None else pd.Series(float("nan"), index=expr.index)

        if present < min_genes_present:
            score = pd.Series(float("nan"), index=expr.index)

        records[f"{state}_score"] = score
        records[f"{state}_n_genes_present"] = present
        records[f"{state}_n_genes_total"] = total
        records[f"{state}_missing_fraction"] = 1.0 - (present / total if total else 0.0)

    records.insert(0, "sample", records.index)
    return records.reset_index(drop=True)


def _weighted_rank_mean(
    rank_lookup: pd.DataFrame,
    genes: list[str],
    weights: dict[str, float] | None,
    index: pd.Index,
) -> pd.Series | None:
    if not genes:
        return None
    if weights is None:
        return rank_lookup[genes].mean(axis=1)
    w = np.array([weights.get(g, 1.0) for g in genes], dtype=float)
    w_sum = w.sum()
    if w_sum == 0:
        return rank_lookup[genes].mean(axis=1)
    return (rank_lookup[genes] * w).sum(axis=1) / w_sum


def ssgsea_state_scores(
    expression: pd.DataFrame,
    signatures: dict[str, StateSignature],
    *,
    min_genes_present: int = 3,
) -> pd.DataFrame:
    """Compute ssGSEA-style enrichment scores per sample and state.

    For each sample, genes are ranked by expression. The enrichment score is
    computed as the difference between a weighted running sum hitting the
    gene set and a uniform penalty for genes outside the set.
    """

    expr = prepare_expression(expression)
    records = pd.DataFrame(index=expr.index)

    for state, signature in signatures.items():
        up_present = set(g for g in signature.up if g in expr.columns)
        down_present = set(g for g in signature.down if g in expr.columns)
        gene_set = up_present | down_present
        total = signature.total_genes
        present = len(gene_set)

        if present < min_genes_present:
            records[f"{state}_score"] = float("nan")
        else:
            scores = []
            for _, row in expr.iterrows():
                vals = row.dropna().sort_values(ascending=False)
                n = len(vals)
                n_hit = sum(1 for g in vals.index if g in up_present)
                n_miss = n - n_hit
                if n_hit == 0 or n_miss == 0:
                    scores.append(0.0)
                    continue
                es = 0.0
                running = 0.0
                for gene in vals.index:
                    if gene in up_present:
                        running += 1.0 / n_hit
                    else:
                        running -= 1.0 / n_miss
                    if abs(running) > abs(es):
                        es = running
                if signature.down and down_present:
                    es_down = 0.0
                    running = 0.0
                    n_hit_d = sum(1 for g in vals.index if g in down_present)
                    n_miss_d = n - n_hit_d
                    if n_hit_d > 0 and n_miss_d > 0:
                        for gene in vals.index:
                            if gene in down_present:
                                running += 1.0 / n_hit_d
                            else:
                                running -= 1.0 / n_miss_d
                            if abs(running) > abs(es_down):
                                es_down = running
                        es = es - es_down
                scores.append(es)
            records[f"{state}_score"] = scores

        records[f"{state}_n_genes_present"] = present
        records[f"{state}_n_genes_total"] = total
        records[f"{state}_missing_fraction"] = 1.0 - (present / total if total else 0.0)

    records.insert(0, "sample", records.index)
    return records.reset_index(drop=True)


def _percentile_ranks_for_genes(expr: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """Compute within-sample percentile ranks only for requested genes."""

    if not genes:
        return pd.DataFrame(index=expr.index)

    matrix = expr.to_numpy(dtype=float, copy=True)
    matrix[~np.isfinite(matrix)] = np.nan
    col_index = {gene: idx for idx, gene in enumerate(expr.columns)}
    requested_indices = [col_index[gene] for gene in genes]
    requested_values = matrix[:, requested_indices]

    row_min = np.nanmin(matrix, axis=1)
    fill_values = np.where(np.isfinite(row_min), row_min - 1.0, 0.0)
    matrix = np.where(np.isfinite(matrix), matrix, fill_values[:, None])
    requested_values = np.where(np.isfinite(requested_values), requested_values, fill_values[:, None])

    sorted_values = np.sort(matrix, axis=1)
    n_samples, n_genes = matrix.shape
    ranks = np.empty((n_samples, len(genes)), dtype=float)

    for row in range(n_samples):
        left = np.searchsorted(sorted_values[row], requested_values[row], side="left")
        right = np.searchsorted(sorted_values[row], requested_values[row], side="right")
        average_rank_1based = (left + 1 + right) / 2.0
        ranks[row] = average_rank_1based / n_genes

    return pd.DataFrame(ranks, index=expr.index, columns=genes)
