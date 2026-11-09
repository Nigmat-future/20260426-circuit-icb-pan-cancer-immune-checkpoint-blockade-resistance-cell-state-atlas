"""Generate the manuscript figures for the CIRCUIT-ICB melanoma paper.

Every plotted value is read from results/tables/melanoma_only/. The only
quantity computed here is the Wilson 95% interval for each calibration bin,
derived from the bin size and observed response fraction in that table.

Sizes follow Bioinformatics column widths (86 mm single, 178 mm double).
Each figure is written as PDF and SVG with editable text, plus a 600 dpi PNG.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "tables" / "melanoma_only"
VAL = RES / "validation"
FIG = ROOT / "manuscript" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

MM = 1 / 25.4
SINGLE = 86 * MM
DOUBLE = 178 * MM

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.linewidth": 0.6,
    "axes.edgecolor": "#272727",
    "axes.labelcolor": "#272727",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.color": "#272727",
    "ytick.color": "#272727",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.major.pad": 2,
    "ytick.major.pad": 2,
    "axes.labelpad": 3,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.dpi": 600,
    "figure.dpi": 150,
})

C = {
    "circuit": "#0F4D92",
    "circuit_mid": "#3775BA",
    "circuit_pale": "#9DB8DA",
    "bench": "#767676",
    "bench_light": "#CFCECE",
    "accent": "#B64342",
    "accent_pale": "#E9A6A1",
    "ink": "#272727",
    "muted": "#4D4D4D",
}

# Biological themes as named in manuscript Table 2.
STATE_THEMES = {
    "CIRCUIT_STATE_01": "Oxidative stress / Wnt",
    "CIRCUIT_STATE_02": "Interferon / ECM",
    "CIRCUIT_STATE_03": "Dedifferentiation",
    "CIRCUIT_STATE_04": "Protein folding",
    "CIRCUIT_STATE_05": "Ribosome biogenesis",
    "CIRCUIT_STATE_06": "Signaling / secretory",
}
BENCH_NAMES = {
    "BENCHMARK_IFNG_ANTIGEN_PRESENTATION_score": "IFNG / antigen presentation",
    "BENCHMARK_TCELL_INFLAMED_GEP_score": "T-cell-inflamed GEP",
    "BENCHMARK_CYTOLYTIC_SCORE_score": "Cytolytic activity",
    "BENCHMARK_CELL_CYCLE_score": "Cell cycle",
    "BENCHMARK_HYPOXIA_score": "Hypoxia",
    "BENCHMARK_EMT_TGFB_score": "EMT / TGF-β",
    "BENCHMARK_HOUSEKEEPING_CONTROL_score": "Housekeeping (control)",
}
COHORT_N = {"GSE78220": 27, "GSE91061": 49}


def state_name(state: str) -> str:
    return state.replace("CIRCUIT_", "").replace("_score", "")


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def add_panel_labels(fig, pairs, dx: float = 0.0, dy: float = 0.01) -> None:
    """Place bold panel letters at the top-left of each axes' full extent."""
    fig.draw_without_rendering()
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    for ax, letter in pairs:
        bb = ax.get_tightbbox(renderer).transformed(inv)
        fig.text(bb.x0 + dx, bb.y1 + dy, letter, fontsize=8, fontweight="bold",
                 ha="left", va="bottom", color=C["ink"])


def save(fig, stem: str) -> None:
    for ext in ("pdf", "svg", "png"):
        fig.savefig(FIG / f"{stem}.{ext}", bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print(f"{stem} saved")


def wilson(k: np.ndarray, n: np.ndarray, z: float = 1.959964) -> tuple[np.ndarray, np.ndarray]:
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return centre - half, centre + half


# ---------------------------------------------------------------------------
# Figure 1. NMF component selection
# ---------------------------------------------------------------------------
def fig1_k_scan() -> None:
    df = read(RES / "melanoma_k_scan.tsv")
    k = df["n_components"].to_numpy()
    stab = df["stability"].to_numpy()
    err = df["reconstruction_error"].to_numpy()
    err_sd = df["reconstruction_error_std"].to_numpy()
    sel = int(df.loc[df["stability"] >= 0.95, "n_components"].max())
    is_sel = k == sel

    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, figsize=(SINGLE, 66 * MM), sharex=True, layout="constrained",
        gridspec_kw={"height_ratios": [1.2, 1]},
    )

    ax_a.axhline(0.95, color=C["bench"], lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax_a.text(k.max(), 0.94, "0.95 threshold", ha="right", va="top", fontsize=5.5, color=C["bench"])
    ax_a.plot(k, stab, color=C["circuit"], lw=1.1, zorder=2)
    ax_a.plot(k, stab, "o", ms=3.6, mfc=C["circuit"], mec="white", mew=0.5, zorder=3)
    ax_a.plot(k[is_sel], stab[is_sel], "o", ms=7.5, mfc="none", mec=C["accent"], mew=1.0, zorder=4)
    ax_a.text(sel - 0.25, 0.80, f"k = {sel}\nselected", ha="right", va="center",
              fontsize=6, color=C["accent"], linespacing=1.1)
    ax_a.set_ylim(0.6, 1.05)
    ax_a.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0])
    ax_a.set_ylabel("Assignment stability")

    ax_b.errorbar(k, err, yerr=err_sd, color=C["muted"], lw=1.0, ms=3.2, marker="o",
                  mfc=C["muted"], mec="white", mew=0.5, elinewidth=0.6, capsize=0, zorder=2)
    ax_b.plot(k[is_sel], err[is_sel], "o", ms=7.0, mfc="none", mec=C["accent"], mew=1.0, zorder=3)
    ax_b.set_ylabel("Reconstruction error")
    ax_b.set_xlabel("Number of components (k)")
    ax_b.set_xticks(k)
    ax_b.set_xlim(k.min() - 0.4, k.max() + 0.4)

    fig.align_ylabels([ax_a, ax_b])
    add_panel_labels(fig, [(ax_a, "A"), (ax_b, "B")])
    save(fig, "fig1_k_scan")


# ---------------------------------------------------------------------------
# Figure 2. Marker specificity per state
# ---------------------------------------------------------------------------
def fig2_marker_heatmap() -> None:
    markers = read(RES / "melanoma_top15_markers.tsv")
    states = sorted(markers["state"].unique())
    top = (markers.sort_values(["state", "rank"]).groupby("state").head(10))
    xmax = np.ceil(top["specificity"].max()) + 0.3

    fig, axes = plt.subplots(2, 3, figsize=(DOUBLE, 96 * MM), sharex=True, layout="constrained")
    for i, (ax, state) in enumerate(zip(axes.flat, states)):
        sub = top[top["state"] == state].sort_values("specificity", ascending=True)
        y = np.arange(len(sub))
        ax.barh(y, sub["specificity"], height=0.66, color=C["circuit_mid"], lw=0, zorder=2)
        ax.axvline(2.0, color=C["accent"], lw=0.6, ls=(0, (3, 2)), zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["gene"], fontstyle="italic")
        ax.tick_params(axis="y", length=0, pad=2)
        ax.set_ylim(-0.6, len(sub) - 0.4 + (0.9 if i == 0 else 0))
        ax.set_xlim(0, xmax)
        ax.set_xticks(np.arange(0, xmax, 2))
        ax.set_title(rf"$\bf{{{state_name(state).replace('_', chr(92) + '_')}}}$   {STATE_THEMES[state]}",
                     loc="left", fontsize=6.5, pad=3)
        if i == 0:
            ax.text(2.1, len(sub) - 0.05, "2.0× filter", ha="left", va="center",
                    fontsize=5.5, color=C["accent"])
    for ax in axes[-1]:
        ax.set_xlabel("Specificity ratio")
    save(fig, "fig2_markers")


# ---------------------------------------------------------------------------
# Figure 3. Pooled AUROC forest plot and permutation null
# ---------------------------------------------------------------------------
def fig3_auroc_and_permutation() -> None:
    uni = read(VAL / "melanoma_univariable_ci.tsv")
    pvals = read(VAL / "melanoma_random_signature_pvalues.tsv").set_index("state")
    null = read(VAL / "melanoma_random_signature_null.tsv")

    pooled = uni[uni["cohort"] == "pooled"].drop_duplicates("state")
    circuit = pooled[pooled["signature_class"] == "circuit_state"].sort_values(
        "auroc_best_direction", ascending=False)
    bench = pooled[pooled["signature_class"] == "benchmark"].copy()
    bench["is_control"] = bench["state"].str.contains("HOUSEKEEPING")
    bench = bench.sort_values(["is_control", "auroc_best_direction"], ascending=[True, False])

    rows = []  # (y, kind, record)
    y = 0.0
    rows.append((y, "header", "CIRCUIT states"))
    for _, r in circuit.iterrows():
        y += 1
        rows.append((y, "circuit", r))
    y += 1.6
    rows.append((y, "header", "Benchmark signatures"))
    for _, r in bench.iterrows():
        y += 1
        rows.append((y, "bench", r))
    y_max = y

    fig = plt.figure(figsize=(DOUBLE, 78 * MM), layout="constrained")
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.25, 1.2], wspace=0.02)
    ax_a = fig.add_subplot(gs[0])
    ax_t = fig.add_subplot(gs[1], sharey=ax_a)
    ax_b = fig.add_subplot(gs[2])

    ticks, labels = [], []
    col_auc, col_dir, col_p = 0.02, 0.62, 0.88
    for yy, kind, rec in rows:
        if kind == "header":
            ticks.append(yy)
            labels.append(rec)
            continue
        auc, lo, hi = rec["auroc_best_direction"], rec["auroc_ci_low"], rec["auroc_ci_high"]
        nr = rec["direction"] == "higher_in_non_responders"
        if kind == "circuit":
            p = float(pvals.loc[rec["state"], "permutation_p"])
            colour = C["circuit"]
            face = colour if p < 0.05 else "white"
            label = state_name(rec["state"])
            p_txt = f"{p:.3f}"
        else:
            p = None
            colour = C["bench"]
            face = "white" if rec["is_control"] else colour
            label = BENCH_NAMES[rec["state"]]
            p_txt = "–"
        ax_a.plot([lo, hi], [yy, yy], color=colour, lw=1.1, solid_capstyle="round", zorder=2)
        ax_a.plot(auc, yy, "o", ms=4.2, mfc=face, mec=colour, mew=0.9, zorder=3)
        ticks.append(yy)
        labels.append(label)
        ax_t.text(col_auc, yy, f"{auc:.3f} ({lo:.3f}–{hi:.3f})", va="center", ha="left", fontsize=6)
        ax_t.text(col_dir, yy, "NR" if nr else "R", va="center", ha="center", fontsize=6,
                  color=C["accent"] if nr else C["circuit_mid"], fontweight="bold")
        ax_t.text(col_p, yy, p_txt, va="center", ha="center", fontsize=6,
                  fontweight="bold" if (p is not None and p < 0.05) else "normal")

    for x, txt in ((col_auc, "AUROC (95% CI)"), (col_dir, "Higher in"), (col_p, "P (perm.)")):
        ax_t.text(x, 0, txt, va="center", ha="left" if x == col_auc else "center",
                  fontsize=6, color=C["muted"], fontweight="bold")

    ax_a.axvline(0.5, color=C["bench"], lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax_a.set_yticks(ticks)
    ax_a.set_yticklabels(labels)
    for tick, (_, kind, _) in zip(ax_a.get_yticklabels(), rows):
        if kind == "header":
            tick.set_fontweight("bold")
            tick.set_fontsize(6.5)
    ax_a.tick_params(axis="y", length=0)
    ax_a.set_ylim(y_max + 0.7, -0.7)
    ax_a.set_xlim(0.35, 0.88)
    ax_a.set_xticks([0.4, 0.5, 0.6, 0.7, 0.8])
    ax_a.spines["left"].set_visible(False)
    ax_a.set_xlabel("Pooled AUROC (best direction)")

    ax_t.set_xlim(0, 1)
    ax_t.axis("off")

    # Panel B: null distribution of size-matched random signatures for STATE_01.
    target = "CIRCUIT_STATE_01_score"
    vals = null.loc[null["state"] == target, "random_auroc_best_direction"].to_numpy()
    obs = float(pvals.loc[target, "auroc_best_direction"])
    p_obs = float(pvals.loc[target, "permutation_p"])
    n_ge = int((vals >= obs).sum())
    bins = np.arange(0.50, vals.max() + 0.0125, 0.0125)
    counts, edges = np.histogram(vals, bins=bins)
    centres = (edges[:-1] + edges[1:]) / 2
    tail = edges[:-1] >= obs - 1e-12
    ax_b.bar(centres, counts, width=np.diff(edges) * 0.92,
             color=np.where(tail, C["accent_pale"], C["bench_light"]), lw=0, zorder=2)
    ax_b.axvline(obs, color=C["circuit"], lw=1.2, zorder=3)
    ax_b.text(obs + 0.006, counts.max() * 0.98,
              f"Observed\n{state_name(target)}\nAUROC = {obs:.3f}\nP = {p_obs:.3f}",
              ha="left", va="top", fontsize=6, color=C["circuit"], linespacing=1.15)
    ax_b.text(obs + 0.006, counts.max() * 0.50, f"{n_ge} of {len(vals):,}\nrandom ≥ observed",
              ha="left", va="top", fontsize=5.5, color=C["accent"], linespacing=1.15)
    ax_b.set_xlim(0.49, max(0.86, vals.max() + 0.01))
    ax_b.set_xlabel("AUROC of size-matched random signature")
    ax_b.set_ylabel("Number of random signatures")

    add_panel_labels(fig, [(ax_a, "A"), (ax_b, "B")])
    save(fig, "fig3_auroc_permutation")


# ---------------------------------------------------------------------------
# Figure 4. Cross-validated comparison and calibration
# ---------------------------------------------------------------------------
def fig4_cv_and_calibration() -> None:
    models = [
        ("CIRCUIT states", "6 features", read(VAL / "melanoma_cv_summary.tsv").iloc[0], C["circuit"]),
        ("Benchmark signatures", "7 features", read(VAL / "melanoma_cv_benchmark_summary.tsv").iloc[0], C["bench"]),
        ("Combined", "13 features", read(VAL / "melanoma_cv_combined_summary.tsv").iloc[0], C["circuit_pale"]),
    ]
    cal = read(VAL / "melanoma_calibration.tsv")

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(DOUBLE, 58 * MM), layout="constrained",
        gridspec_kw={"width_ratios": [1.8, 1]},
    )

    for i, (name, feat, s, colour) in enumerate(models):
        m, lo, hi = s["cv_auroc_mean"], s["cv_auroc_ci_low"], s["cv_auroc_ci_high"]
        ax_a.plot([lo, hi], [i, i], color=colour, lw=2.2, solid_capstyle="round", zorder=2)
        ax_a.plot(m, i, "o", ms=6, mfc=colour, mec="white", mew=0.8, zorder=3)
        ax_a.text(0.80, i, f"{m:.3f}  ({lo:.3f}–{hi:.3f})", va="center", ha="left", fontsize=6)
    ax_a.text(0.80, -0.75, "Mean  (95% interval)", va="center", ha="left", fontsize=6,
              color=C["muted"], fontweight="bold")
    ax_a.axvline(0.5, color=C["bench"], lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax_a.text(0.505, 2.45, "chance", fontsize=5.5, color=C["bench"], va="center", ha="left")
    ax_a.set_yticks(range(len(models)))
    ax_a.set_yticklabels([f"{n}\n{f}" for n, f, _, _ in models], linespacing=1.2)
    for tick in ax_a.get_yticklabels():
        tick.set_multialignment("right")
    ax_a.tick_params(axis="y", length=0)
    ax_a.set_ylim(2.6, -1.05)
    ax_a.set_xlim(0.40, 1.0)
    ax_a.set_xticks([0.4, 0.5, 0.6, 0.7, 0.8])
    ax_a.spines["bottom"].set_bounds(0.40, 0.80)
    ax_a.spines["left"].set_visible(False)
    ax_a.set_xlabel("Cross-validated AUROC (50 × 5-fold)")

    n = cal["n"].to_numpy(dtype=float)
    frac = cal["fraction_positive"].to_numpy(dtype=float)
    pred = cal["mean_predicted"].to_numpy(dtype=float)
    k = np.round(frac * n)
    lo, hi = wilson(k, n)
    ax_b.plot([0, 1], [0, 1], color=C["bench"], lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax_b.text(0.88, 0.895, "perfect calibration", rotation=45, rotation_mode="anchor",
              ha="center", va="bottom", fontsize=5.5, color=C["bench"])
    ax_b.plot(pred, frac, color=C["circuit"], lw=0.8, alpha=0.6, zorder=2)
    ax_b.errorbar(pred, frac, yerr=[frac - lo, hi - frac], fmt="o", color=C["circuit"],
                  ms=4.2, mec="white", mew=0.6, elinewidth=0.8, capsize=0, zorder=3)
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1)
    ax_b.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax_b.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax_b.set_aspect("equal")
    ax_b.set_xlabel("Mean predicted probability")
    ax_b.set_ylabel("Observed response rate")

    add_panel_labels(fig, [(ax_a, "A"), (ax_b, "B")])
    save(fig, "fig4_cv_calibration")


# ---------------------------------------------------------------------------
# Figure 5. Direction of association per cohort
# ---------------------------------------------------------------------------
def fig5_direction_stability() -> None:
    uni = read(VAL / "melanoma_univariable_ci.tsv").drop_duplicates(["cohort", "state"])
    stab = read(VAL / "melanoma_direction_stability.tsv").set_index("state")
    per = uni[uni["cohort"].isin(COHORT_N) & (uni["signature_class"] == "circuit_state")].copy()

    # Orient every AUROC and its bootstrap interval toward responders (exact: the
    # interval was computed in a fixed direction and flipped for non-responder states).
    nr = per["direction"].eq("higher_in_non_responders")
    per["auc_r"] = per["auroc_high_responder"]
    per["lo_r"] = np.where(nr, 1 - per["auroc_ci_high"], per["auroc_ci_low"])
    per["hi_r"] = np.where(nr, 1 - per["auroc_ci_low"], per["auroc_ci_high"])

    states = sorted(per["state"].unique())
    offsets = {"GSE78220": -0.17, "GSE91061": 0.17}
    shapes = {"GSE78220": "o", "GSE91061": "s"}

    fig, ax = plt.subplots(figsize=(SINGLE, 72 * MM), layout="constrained")
    ax.axvspan(0, 0.5, color="#F7F7F7", lw=0, zorder=0)
    ax.axvline(0.5, color=C["bench"], lw=0.6, ls=(0, (3, 2)), zorder=1)
    for i, state in enumerate(states):
        stable = bool(stab.loc[state, "is_stable"])
        colour = C["circuit"] if stable else C["bench"]
        sub = per[per["state"] == state].set_index("cohort")
        pts = []
        for cohort, off in offsets.items():
            r = sub.loc[cohort]
            yy = i + off
            ax.plot([r["lo_r"], r["hi_r"]], [yy, yy], color=colour, lw=0.8, alpha=0.55,
                    solid_capstyle="round", zorder=2)
            ax.plot(r["auc_r"], yy, shapes[cohort], ms=4.2, mfc=colour, mec="white", mew=0.5, zorder=4)
            pts.append((r["auc_r"], yy))
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=colour, lw=0.6, zorder=3)
        agree = int(stab.loc[state, "n_agree"])
        total = int(stab.loc[state, "n_cohorts"])
        ax.text(1.02, i, f"{agree}/{total}", transform=ax.get_yaxis_transform(), va="center",
                ha="left", fontsize=6, color=colour, fontweight="bold" if stable else "normal")
    ax.text(1.02, -0.85, "Agree", transform=ax.get_yaxis_transform(), va="center", ha="left",
            fontsize=6, color=C["muted"], fontweight="bold")
    ax.text(0.03, -0.85, "← higher in non-responders", va="center", ha="left", fontsize=5.5,
            color=C["accent"])
    ax.text(0.97, -0.85, "higher in responders →", va="center", ha="right", fontsize=5.5,
            color=C["circuit_mid"])
    ax.set_yticks(range(len(states)))
    ax.set_yticklabels([state_name(s) for s in states])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(states) - 0.45, -1.3)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Per-cohort AUROC (oriented to responders)")

    handles = [Line2D([], [], ls="", marker=shapes[c], ms=4.2, mfc=C["muted"], mec="white", mew=0.5,
                      label=f"{c} (n = {COHORT_N[c]})") for c in offsets]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, handletextpad=0.3,
               columnspacing=1.2, fontsize=6)
    save(fig, "fig5_direction_stability")


if __name__ == "__main__":
    fig1_k_scan()
    fig2_marker_heatmap()
    fig3_auroc_and_permutation()
    fig4_cv_and_calibration()
    fig5_direction_stability()
    print(f"\nAll figures saved to {FIG}")
