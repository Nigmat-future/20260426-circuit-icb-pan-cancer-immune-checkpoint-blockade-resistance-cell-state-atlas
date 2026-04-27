from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

try:
    from sklearn.linear_model import SGDClassifier
except ImportError:
    SGDClassifier = None  # type: ignore[assignment,misc]

MIN_COHORT_SIZE = 10


RESPONDER_LABELS = {
    "CR",
    "PR",
    "PRCR",
    "COMPLETE RESPONSE",
    "PARTIAL RESPONSE",
    "R",
    "E",
    "DCB",
    "FAVOURABLE",
    "FAVORABLE",
    "RESPONDER",
    "ICB_PR",
}

NON_RESPONDER_LABELS = {
    "PD",
    "SD",
    "PROGRESSIVE DISEASE",
    "STABLE DISEASE",
    "NR",
    "NE",
    "NDB",
    "UNFAVOURABLE",
    "UNFAVORABLE",
    "NON_RESPONDER",
    "NON-RESPONDER",
    "NONRESPONDER",
    "RESIS",
    "RESISTANT",
    "POST-ICI (RESISTANT)",
}


def harmonize_response(value: object) -> str:
    if pd.isna(value):
        return "unknown"
    text = str(value).strip()
    upper = text.upper()
    if upper in RESPONDER_LABELS:
        return "responder"
    if upper in NON_RESPONDER_LABELS:
        return "non_responder"
    if upper in {"ICB_SD"}:
        return "non_responder"
    if upper in {"STABLE_OR_MIXED"}:
        return "stable_or_mixed"
    if upper in {"UT", "UNTREATED", "NOICB", "UNKNOWN", "NAN", ""}:
        return "unknown"
    return "unknown"


def state_columns(df: pd.DataFrame) -> list[str]:
    return [col for col in df.columns if col.endswith("_activity") or col.endswith("_score")]


def evaluate_univariable_states(
    df: pd.DataFrame,
    *,
    outcome_col: str,
    cohort_col: str | None = None,
) -> pd.DataFrame:
    """Compute univariable AUROC/AUPRC for every state score."""

    work = df.copy()
    work["response_binary"] = work[outcome_col].map(harmonize_response)
    work = work.loc[work["response_binary"].isin(["responder", "non_responder"])].copy()
    work["y"] = (work["response_binary"] == "responder").astype(int)

    rows = []
    for score_col in state_columns(work):
        values = pd.to_numeric(work[score_col], errors="coerce")
        mask = values.notna()
        if mask.sum() < 4 or work.loc[mask, "y"].nunique() < 2:
            continue
        auc = roc_auc_score(work.loc[mask, "y"], values.loc[mask])
        auprc = average_precision_score(work.loc[mask, "y"], values.loc[mask])
        rows.append(
            {
                "cohort": "pooled" if cohort_col is None else "pooled",
                "state": score_col,
                "n": int(mask.sum()),
                "n_responders": int(work.loc[mask, "y"].sum()),
                "n_non_responders": int(mask.sum() - work.loc[mask, "y"].sum()),
                "auroc_high_responder": auc,
                "auroc_best_direction": max(auc, 1.0 - auc),
                "direction": "higher_in_responders" if auc >= 0.5 else "higher_in_non_responders",
                "auprc_high_responder": auprc,
            }
        )

    if cohort_col is not None and cohort_col in work.columns:
        for cohort, group in work.groupby(cohort_col):
            if group["y"].nunique() < 2:
                continue
            for score_col in state_columns(group):
                values = pd.to_numeric(group[score_col], errors="coerce")
                mask = values.notna()
                if mask.sum() < 4 or group.loc[mask, "y"].nunique() < 2:
                    continue
                auc = roc_auc_score(group.loc[mask, "y"], values.loc[mask])
                auprc = average_precision_score(group.loc[mask, "y"], values.loc[mask])
                rows.append(
                    {
                        "cohort": cohort,
                        "state": score_col,
                        "n": int(mask.sum()),
                        "n_responders": int(group.loc[mask, "y"].sum()),
                        "n_non_responders": int(mask.sum() - group.loc[mask, "y"].sum()),
                        "auroc_high_responder": auc,
                        "auroc_best_direction": max(auc, 1.0 - auc),
                        "direction": "higher_in_responders" if auc >= 0.5 else "higher_in_non_responders",
                        "auprc_high_responder": auprc,
                    }
                )
    return pd.DataFrame(rows)


def leave_one_cohort_logistic(
    df: pd.DataFrame,
    *,
    outcome_col: str,
    cohort_col: str,
    min_test_size: int = MIN_COHORT_SIZE,
    use_elastic_net: bool = True,
) -> pd.DataFrame:
    """Run leave-one-cohort-out logistic validation using all state columns."""

    work = df.copy()
    work["response_binary"] = work[outcome_col].map(harmonize_response)
    work = work.loc[work["response_binary"].isin(["responder", "non_responder"])].copy()
    work["y"] = (work["response_binary"] == "responder").astype(int)
    scores = state_columns(work)
    rows = []

    for cohort in sorted(work[cohort_col].dropna().unique()):
        train = work.loc[work[cohort_col] != cohort].copy()
        test = work.loc[work[cohort_col] == cohort].copy()
        if train["y"].nunique() < 2 or test["y"].nunique() < 2 or len(test) < min_test_size:
            rows.append(
                {
                    "heldout_cohort": cohort,
                    "status": "skipped",
                    "reason": f"insufficient samples (n={len(test)}, min={min_test_size})" if len(test) < min_test_size else "insufficient class balance",
                    "n_train": len(train),
                    "n_test": len(test),
                }
            )
            continue
        x_train = train[scores].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        x_test = test[scores].apply(pd.to_numeric, errors="coerce").fillna(0.0)
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
            model = make_pipeline(
                StandardScaler(),
                LogisticRegression(solver="liblinear", class_weight="balanced", random_state=7),
            )
        model.fit(x_train, train["y"])
        if hasattr(model[-1], "predict_proba"):
            prob = model.predict_proba(x_test)[:, 1]
        else:
            prob = model.decision_function(x_test)
            prob = 1.0 / (1.0 + np.exp(-prob))
        rows.append(
            {
                "heldout_cohort": cohort,
                "status": "evaluated",
                "reason": "",
                "n_train": len(train),
                "n_test": len(test),
                "n_test_responders": int(test["y"].sum()),
                "n_test_non_responders": int(len(test) - test["y"].sum()),
                "auroc": roc_auc_score(test["y"], prob),
                "auprc": average_precision_score(test["y"], prob),
                "mean_predicted_probability": float(np.mean(prob)),
                "model_type": "elastic_net" if use_elastic_net else "logistic",
            }
        )
    return pd.DataFrame(rows)
