from __future__ import annotations

import math

import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

try:
    from sklearn.linear_model import SGDClassifier
except ImportError:
    SGDClassifier = None  # type: ignore[assignment,misc]

from circuit_icb.signature import prepare_expression
from circuit_icb.state_validation import harmonize_response, state_columns

MIN_GROUP_SIZE_DEFAULT = 10


def bootstrap_metric_ci(
    y: pd.Series,
    scores: pd.Series,
    *,
    metric: str = "auroc",
    n_bootstrap: int = 1000,
    seed: int = 7,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    y_values = np.asarray(y, dtype=int)
    score_values = np.asarray(scores, dtype=float)
    metrics = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, len(y_values), len(y_values))
        if len(np.unique(y_values[idx])) < 2:
            continue
        if metric == "auprc":
            metrics.append(average_precision_score(y_values[idx], score_values[idx]))
        else:
            metrics.append(roc_auc_score(y_values[idx], score_values[idx]))
    if not metrics:
        return float("nan"), float("nan")
    return float(np.quantile(metrics, 0.025)), float(np.quantile(metrics, 0.975))


def evaluate_states_with_ci(
    scores: pd.DataFrame,
    *,
    outcome_col: str = "response_original",
    cohort_col: str = "cohort_id",
    n_bootstrap: int = 1000,
    seed: int = 7,
) -> pd.DataFrame:
    work = scores.copy()
    work["response_binary"] = work[outcome_col].map(harmonize_response)
    work = work.loc[work["response_binary"].isin(["responder", "non_responder"])].copy()
    work["y"] = (work["response_binary"] == "responder").astype(int)
    rows = []
    groups: list[tuple[str, pd.DataFrame]] = [("pooled", work)]
    if cohort_col in work.columns:
        groups.extend((str(name), group) for name, group in work.groupby(cohort_col))
    for cohort, group in groups:
        if group["y"].nunique() < 2:
            continue
        for score_col in state_columns(group):
            values = pd.to_numeric(group[score_col], errors="coerce")
            mask = values.notna()
            if mask.sum() < 4 or group.loc[mask, "y"].nunique() < 2:
                continue
            y = group.loc[mask, "y"]
            x = values.loc[mask]
            auc = roc_auc_score(y, x)
            auprc = average_precision_score(y, x)
            auc_low_raw, auc_high_raw = bootstrap_metric_ci(y, x, metric="auroc", n_bootstrap=n_bootstrap, seed=seed)
            pr_low, pr_high = bootstrap_metric_ci(y, x, metric="auprc", n_bootstrap=n_bootstrap, seed=seed + 1)
            if auc >= 0.5:
                best_auc = auc
                direction = "higher_in_responders"
                auc_low = auc_low_raw
                auc_high = auc_high_raw
            else:
                best_auc = 1.0 - auc
                direction = "higher_in_non_responders"
                auc_low = 1.0 - auc_high_raw if np.isfinite(auc_high_raw) else float("nan")
                auc_high = 1.0 - auc_low_raw if np.isfinite(auc_low_raw) else float("nan")
            rows.append(
                {
                    "cohort": cohort,
                    "state": score_col,
                    "n": int(mask.sum()),
                    "n_responders": int(y.sum()),
                    "n_non_responders": int(mask.sum() - y.sum()),
                    "auroc_high_responder": auc,
                    "auroc_best_direction": best_auc,
                    "direction": direction,
                    "auroc_ci_low": auc_low,
                    "auroc_ci_high": auc_high,
                    "auprc_high_responder": auprc,
                    "auprc_ci_low": pr_low,
                    "auprc_ci_high": pr_high,
                    "signature_class": "benchmark" if score_col.startswith("BENCHMARK_") else "circuit_state",
                }
            )
    return pd.DataFrame(rows)


def fixed_or_random_effects_meta(per_cohort: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for state, group in per_cohort.loc[per_cohort["cohort"].ne("pooled")].groupby("state"):
        g = group.dropna(subset=["auroc_high_responder", "n"])
        if g.shape[0] < 2:
            continue
        effects = g["auroc_high_responder"].astype(float).to_numpy() - 0.5
        variances = np.maximum(1.0 / (12.0 * g["n"].astype(float).to_numpy()), 1e-6)
        weights = 1.0 / variances
        fixed = float(np.sum(weights * effects) / np.sum(weights))
        q = float(np.sum(weights * (effects - fixed) ** 2))
        df = max(len(effects) - 1, 1)
        c = float(np.sum(weights) - (np.sum(weights**2) / np.sum(weights)))
        tau2 = max((q - df) / c, 0.0) if c > 0 else 0.0
        random_weights = 1.0 / (variances + tau2)
        pooled = float(np.sum(random_weights * effects) / np.sum(random_weights))
        se = math.sqrt(1.0 / np.sum(random_weights))
        z = pooled / se if se else float("nan")
        pvalue = math.erfc(abs(z) / math.sqrt(2.0)) if np.isfinite(z) else float("nan")
        rows.append(
            {
                "state": state,
                "n_cohorts": int(g.shape[0]),
                "random_effect_auc_delta": pooled,
                "random_effect_auroc": pooled + 0.5,
                "se": se,
                "ci_low": pooled + 0.5 - 1.96 * se,
                "ci_high": pooled + 0.5 + 1.96 * se,
                "tau2": tau2,
                "q": q,
                "pvalue": pvalue,
                "directions_higher_responder": int((g["direction"] == "higher_in_responders").sum()),
                "directions_higher_non_responder": int((g["direction"] == "higher_in_non_responders").sum()),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    order = out["pvalue"].rank(method="first").astype(float)
    out["fdr_bh"] = (out["pvalue"] * len(out) / order).clip(upper=1.0)
    out = out.sort_values(["fdr_bh", "random_effect_auroc"], ascending=[True, False])
    return out


def leave_one_group_logistic(
    df: pd.DataFrame,
    *,
    group_col: str,
    outcome_col: str = "response_original",
    feature_prefixes: tuple[str, ...] = ("CIRCUIT_STATE_",),
    min_test_size: int = MIN_GROUP_SIZE_DEFAULT,
    use_elastic_net: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    work = df.copy()
    work["response_binary"] = work[outcome_col].map(harmonize_response)
    work = work.loc[work["response_binary"].isin(["responder", "non_responder"])].copy()
    work["y"] = (work["response_binary"] == "responder").astype(int)
    features = [col for col in state_columns(work) if col.startswith(feature_prefixes)]
    rows = []
    pred_rows = []
    if group_col not in work.columns or not features:
        return pd.DataFrame(rows), pd.DataFrame(pred_rows)
    for group_name in sorted(work[group_col].dropna().astype(str).unique()):
        train = work.loc[work[group_col].astype(str) != group_name].copy()
        test = work.loc[work[group_col].astype(str) == group_name].copy()
        if train["y"].nunique() < 2 or test["y"].nunique() < 2 or len(test) < min_test_size:
            rows.append({"heldout_group": group_name, "group_col": group_col, "status": "skipped", "reason": f"insufficient samples (n={len(test)}, min={min_test_size})" if len(test) < min_test_size else "insufficient class balance", "n_train": len(train), "n_test": len(test)})
            continue
        x_train = train[features].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        x_test = test[features].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        if use_elastic_net and SGDClassifier is not None:
            model = make_pipeline(
                StandardScaler(),
                SGDClassifier(
                    loss="log_loss",
                    penalty="elasticnet",
                    l1_ratio=0.5,
                    alpha=1e-3,
                    class_weight="balanced",
                    random_state=7,
                    max_iter=2000,
                ),
            )
        else:
            model = make_pipeline(StandardScaler(), LogisticRegression(solver="liblinear", class_weight="balanced", random_state=7))
        model.fit(x_train, train["y"])
        if hasattr(model[-1], "predict_proba"):
            prob = model.predict_proba(x_test)[:, 1]
        else:
            prob = model.decision_function(x_test)
            prob = 1.0 / (1.0 + np.exp(-prob))
        rows.append(
            {
                "heldout_group": group_name,
                "group_col": group_col,
                "status": "evaluated",
                "reason": "",
                "n_train": len(train),
                "n_test": len(test),
                "n_test_responders": int(test["y"].sum()),
                "n_test_non_responders": int(len(test) - test["y"].sum()),
                "auroc": roc_auc_score(test["y"], prob),
                "auprc": average_precision_score(test["y"], prob),
                "brier": brier_score_loss(test["y"], prob),
                "model_type": "elastic_net" if use_elastic_net else "logistic",
            }
        )
        for sample, y_value, p_value in zip(test["sample"], test["y"], prob):
            pred_rows.append({"sample": sample, "heldout_group": group_name, "group_col": group_col, "y": int(y_value), "predicted_probability": float(p_value)})
    return pd.DataFrame(rows), pd.DataFrame(pred_rows)


def calibration_table(predictions: pd.DataFrame, *, n_bins: int = 5) -> pd.DataFrame:
    if predictions.empty or predictions["y"].nunique() < 2:
        return pd.DataFrame(columns=["bin_index", "mean_predicted", "fraction_positive", "n"])
    frac_pos, mean_pred = calibration_curve(predictions["y"], predictions["predicted_probability"], n_bins=n_bins, strategy="quantile")
    bins = pd.qcut(predictions["predicted_probability"], q=min(n_bins, predictions.shape[0]), duplicates="drop")
    counts = bins.value_counts().sort_index().to_numpy()
    return pd.DataFrame(
        {
            "bin_index": range(1, len(mean_pred) + 1),
            "mean_predicted": mean_pred,
            "fraction_positive": frac_pos,
            "n": counts[: len(mean_pred)],
        }
    )


def random_signature_null(
    expression: pd.DataFrame,
    clinical: pd.DataFrame,
    state_score_col: str,
    *,
    outcome_col: str = "response_original",
    signature_size: int,
    excluded_genes: set[str],
    n_permutations: int = 1000,
    seed: int = 7,
) -> pd.DataFrame:
    expr = prepare_expression(expression)
    merged = clinical[["sample", outcome_col]].copy()
    merged["response_binary"] = merged[outcome_col].map(harmonize_response)
    valid_samples = merged.loc[merged["response_binary"].isin(["responder", "non_responder"]), "sample"].astype(str)
    expr = expr.loc[expr.index.astype(str).isin(set(valid_samples))].copy()
    merged = merged.loc[merged["sample"].astype(str).isin(set(expr.index.astype(str)))].copy()
    if expr.shape[0] < 4 or merged["response_binary"].nunique() < 2:
        return pd.DataFrame()
    y = (merged.set_index("sample").loc[expr.index, "response_binary"] == "responder").astype(int)
    available = [gene for gene in expr.columns if gene not in excluded_genes]
    if len(available) < signature_size:
        return pd.DataFrame()
    rng = np.random.default_rng(seed)
    rank_matrix = expr.rank(axis=1, pct=True, method="average")
    rows = []
    for idx in range(n_permutations):
        genes = rng.choice(available, size=signature_size, replace=False).tolist()
        values = rank_matrix[genes].mean(axis=1, skipna=True)
        mask = values.notna()
        if mask.sum() < 4 or y.loc[mask].nunique() < 2:
            continue
        auc = roc_auc_score(y.loc[mask], values.loc[mask])
        rows.append({"state": state_score_col, "random_index": idx, "random_auroc_best_direction": max(auc, 1.0 - auc)})
    return pd.DataFrame(rows)


def direction_stability(per_cohort: pd.DataFrame) -> pd.DataFrame:
    """Flag states whose direction flips across cohorts.

    Returns one row per state with the fraction of cohorts agreeing on
    the majority direction and a stability flag.
    """
    cohort_only = per_cohort.loc[per_cohort["cohort"].ne("pooled")].copy()
    rows = []
    for state, group in cohort_only.groupby("state"):
        n = len(group)
        if n == 0:
            continue
        n_higher_resp = int((group["direction"] == "higher_in_responders").sum())
        n_higher_non = n - n_higher_resp
        majority = "higher_in_responders" if n_higher_resp >= n_higher_non else "higher_in_non_responders"
        agreement = max(n_higher_resp, n_higher_non) / n
        rows.append({
            "state": state,
            "n_cohorts": n,
            "majority_direction": majority,
            "n_agree": max(n_higher_resp, n_higher_non),
            "n_disagree": min(n_higher_resp, n_higher_non),
            "direction_agreement": agreement,
            "is_stable": agreement >= 0.8,
        })
    return pd.DataFrame(rows)


def circuit_vs_benchmark_comparison(per_cohort: pd.DataFrame) -> pd.DataFrame:
    """Side-by-side comparison of CIRCUIT_STATE vs BENCHMARK signatures per cohort."""
    cohort_only = per_cohort.loc[per_cohort["cohort"].ne("pooled")].copy()
    if cohort_only.empty:
        return pd.DataFrame()

    rows = []
    for cohort, group in cohort_only.groupby("cohort"):
        circuit = group.loc[group["signature_class"] == "circuit_state"]
        bench = group.loc[group["signature_class"] == "benchmark"]
        circuit_best = circuit["auroc_best_direction"].max() if not circuit.empty else float("nan")
        circuit_mean = circuit["auroc_best_direction"].mean() if not circuit.empty else float("nan")
        bench_best = bench["auroc_best_direction"].max() if not bench.empty else float("nan")
        bench_mean = bench["auroc_best_direction"].mean() if not bench.empty else float("nan")
        rows.append({
            "cohort": cohort,
            "n_samples": int(group["n"].iloc[0]) if "n" in group.columns else 0,
            "circuit_best_auroc": circuit_best,
            "circuit_mean_auroc": circuit_mean,
            "benchmark_best_auroc": bench_best,
            "benchmark_mean_auroc": bench_mean,
            "gap_best": bench_best - circuit_best if np.isfinite(bench_best) and np.isfinite(circuit_best) else float("nan"),
            "gap_mean": bench_mean - circuit_mean if np.isfinite(bench_mean) and np.isfinite(circuit_mean) else float("nan"),
        })
    return pd.DataFrame(rows)


def repeated_stratified_cv(
    df: pd.DataFrame,
    *,
    outcome_col: str = "response_original",
    feature_prefixes: tuple[str, ...] = ("CIRCUIT_STATE_",),
    n_splits: int = 5,
    n_repeats: int = 20,
    seed: int = 7,
    use_elastic_net: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Repeated stratified k-fold CV on pooled samples.

    Returns (fold_stats, per_sample_predictions).
    Each sample gets an out-of-fold predicted probability averaged across
    the repeats where it appeared in the test fold.
    """
    work = df.copy()
    work["response_binary"] = work[outcome_col].map(harmonize_response)
    work = work.loc[work["response_binary"].isin(["responder", "non_responder"])].copy()
    work["y"] = (work["response_binary"] == "responder").astype(int)
    features = [col for col in state_columns(work) if col.startswith(feature_prefixes)]
    if not features or work["y"].nunique() < 2:
        return pd.DataFrame(), pd.DataFrame()

    x = work[features].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy()
    y = work["y"].to_numpy()
    samples = work["sample"].to_numpy()

    rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=seed)

    fold_rows = []
    sample_preds: dict[str, list[float]] = {s: [] for s in samples}

    for fold_idx, (train_idx, test_idx) in enumerate(rskf.split(x, y)):
        if len(np.unique(y[train_idx])) < 2 or len(np.unique(y[test_idx])) < 2:
            continue
        if use_elastic_net and SGDClassifier is not None:
            model = make_pipeline(
                StandardScaler(),
                SGDClassifier(
                    loss="log_loss", penalty="elasticnet",
                    l1_ratio=0.5, alpha=1e-3,
                    class_weight="balanced", random_state=seed,
                    max_iter=2000,
                ),
            )
        else:
            model = make_pipeline(
                StandardScaler(),
                LogisticRegression(solver="liblinear", class_weight="balanced", random_state=seed),
            )
        model.fit(x[train_idx], y[train_idx])
        if hasattr(model[-1], "predict_proba"):
            prob = model.predict_proba(x[test_idx])[:, 1]
        else:
            dec = model.decision_function(x[test_idx])
            prob = 1.0 / (1.0 + np.exp(-dec))

        auc = roc_auc_score(y[test_idx], prob)
        auprc = average_precision_score(y[test_idx], prob)
        brier = brier_score_loss(y[test_idx], prob)
        fold_rows.append({
            "fold": fold_idx,
            "repeat": fold_idx // n_splits,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "auroc": auc,
            "auprc": auprc,
            "brier": brier,
        })
        for i, ti in enumerate(test_idx):
            sample_preds[samples[ti]].append(float(prob[i]))

    fold_stats = pd.DataFrame(fold_rows)

    pred_rows = []
    for sample_name, preds in sample_preds.items():
        if preds:
            mask = samples == sample_name
            idx = np.where(mask)[0][0]
            pred_rows.append({
                "sample": sample_name,
                "y": int(y[idx]),
                "mean_predicted_probability": float(np.mean(preds)),
                "std_predicted_probability": float(np.std(preds)),
                "n_folds": len(preds),
            })
    predictions = pd.DataFrame(pred_rows)

    return fold_stats, predictions


def repeated_cv_summary(fold_stats: pd.DataFrame) -> dict[str, float]:
    """Summarize repeated CV results into mean +/- std across repeats."""
    if fold_stats.empty:
        return {}
    per_repeat = fold_stats.groupby("repeat").agg(
        auroc_mean=("auroc", "mean"),
        auprc_mean=("auprc", "mean"),
        brier_mean=("brier", "mean"),
    )
    return {
        "cv_auroc_mean": float(per_repeat["auroc_mean"].mean()),
        "cv_auroc_std": float(per_repeat["auroc_mean"].std()),
        "cv_auroc_ci_low": float(np.quantile(per_repeat["auroc_mean"], 0.025)),
        "cv_auroc_ci_high": float(np.quantile(per_repeat["auroc_mean"], 0.975)),
        "cv_auprc_mean": float(per_repeat["auprc_mean"].mean()),
        "cv_auprc_std": float(per_repeat["auprc_mean"].std()),
        "cv_brier_mean": float(per_repeat["brier_mean"].mean()),
        "cv_brier_std": float(per_repeat["brier_mean"].std()),
        "n_repeats": int(per_repeat.shape[0]),
        "n_folds_total": int(fold_stats.shape[0]),
    }
