"""Generate all manuscript figures for CIRCUIT-ICB melanoma paper."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

RES = ROOT / "results" / "tables" / "melanoma_only"
VAL = RES / "validation"
FIG = ROOT / "manuscript" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# Global style
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "legend.frameon": False,
    "legend.fontsize": 7,
    "figure.dpi": 300,
})
PALETTE_CIRCUIT = "#2166AC"
PALETTE_BENCH = "#B2182B"
PALETTE_COMBINED = "#4DAF4A"
STATE_COLORS = sns.color_palette("Set2", 6)


def fig1_k_scan():
    """Figure 1: NMF stability and reconstruction error vs k."""
    df = pd.read_csv(RES / "melanoma_k_scan.tsv", sep="\t")

    fig, ax1 = plt.subplots(figsize=(3.5, 2.5))
    color_err = "#333333"
    color_stab = PALETTE_CIRCUIT

    ax1.set_xlabel("Number of components (k)")
    ax1.set_ylabel("Reconstruction error", color=color_err)
    bars = ax1.bar(df["n_components"], df["reconstruction_error"],
                   color="#CCCCCC", edgecolor=color_err, linewidth=0.5, width=0.6, zorder=1)
    ax1.tick_params(axis="y", labelcolor=color_err)
    ax1.set_ylim(200, 310)

    ax2 = ax1.twinx()
    ax2.set_ylabel("Assignment stability", color=color_stab)
    ax2.plot(df["n_components"], df["stability"], "o-", color=color_stab,
             linewidth=1.5, markersize=5, zorder=3)
    ax2.axhline(0.95, color=color_stab, linestyle="--", linewidth=0.7, alpha=0.5)
    ax2.set_ylim(0.5, 1.05)
    ax2.tick_params(axis="y", labelcolor=color_stab)

    # Mark selected k
    sel = df.loc[df["stability"] >= 0.95, "n_components"].max()
    ax1.axvline(sel + 0.5, color="#666666", linestyle=":", linewidth=0.8)
    ax1.text(sel + 0.7, 208, f"k={int(sel)} selected", fontsize=6.5, color="#444444")

    fig.tight_layout()
    fig.savefig(FIG / "fig1_k_scan.pdf", bbox_inches="tight")
    fig.savefig(FIG / "fig1_k_scan.png", bbox_inches="tight")
    plt.close()
    print("Figure 1 saved")


def fig2_marker_heatmap():
    """Figure 2: Top marker genes heatmap per state with specificity."""
    markers = pd.read_csv(RES / "melanoma_top15_markers.tsv", sep="\t")
    states = sorted(markers["state"].unique())

    fig, axes = plt.subplots(1, len(states), figsize=(9, 3.2), sharey=False)
    for i, state in enumerate(states):
        ax = axes[i]
        sub = markers[markers["state"] == state].head(10)
        genes = sub["gene"].tolist()[::-1]
        specs = sub["specificity"].tolist()[::-1]
        colors = [sns.light_palette(STATE_COLORS[i], reverse=True)[-1] if s >= 3
                  else STATE_COLORS[i] for s in specs]
        ax.barh(range(len(genes)), specs, color=colors, edgecolor="none", height=0.7)
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes, fontsize=5.5, fontfamily="monospace")
        ax.set_xlim(0, 8)
        ax.set_title(state.replace("CIRCUIT_STATE_", "S"), fontsize=7, fontweight="bold")
        if i == 0:
            ax.set_xlabel("Specificity ratio")
        else:
            ax.set_xlabel("")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.text(0.5, -0.02, "Specificity ratio (loading in state / mean loading across states)", ha="center", fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_markers.pdf", bbox_inches="tight")
    fig.savefig(FIG / "fig2_markers.png", bbox_inches="tight")
    plt.close()
    print("Figure 2 saved")


def fig3_auroc_and_permutation():
    """Figure 3: Forest plot (A) and permutation null (B)."""
    uni = pd.read_csv(VAL / "melanoma_univariable_ci.tsv", sep="\t")
    pvals = pd.read_csv(VAL / "melanoma_random_signature_pvalues.tsv", sep="\t")

    pooled = uni[uni["cohort"] == "pooled"].copy()
    circuit = pooled[pooled["signature_class"] == "circuit_state"].copy()
    bench = pooled[pooled["signature_class"] == "benchmark"].copy()

    circuit = circuit.sort_values("auroc_best_direction", ascending=False)
    bench = bench.sort_values("auroc_best_direction", ascending=False)

    # --- Panel A: Forest plot ---
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.5, 3.5), gridspec_kw={"width_ratios": [1.4, 1]})

    all_sigs = pd.concat([circuit, bench])
    labels = []
    for _, r in all_sigs.iterrows():
        name = r["state"].replace("CIRCUIT_STATE_", "CIRCUIT ").replace("_score", "")
        name = name.replace("BENCHMARK_", "BENCH ").replace("_score", "")
        if len(name) > 22:
            name = name[:22]
        labels.append(name)

    y_pos = range(len(labels))
    colors = [PALETTE_CIRCUIT if r["signature_class"] == "circuit_state" else PALETTE_BENCH
              for _, r in all_sigs.iterrows()]

    for i, (_, r) in enumerate(all_sigs.iterrows()):
        auc = r["auroc_best_direction"]
        ci_lo = r["auroc_ci_low"]
        ci_hi = r["auroc_ci_high"]
        ax_a.plot([ci_lo, ci_hi], [i, i], color=colors[i], linewidth=1.5, solid_capstyle="round")
        ax_a.plot(auc, i, "o", color=colors[i], markersize=5)

    ax_a.axvline(0.5, color="#999999", linestyle="--", linewidth=0.5)
    ax_a.set_yticks(list(y_pos))
    ax_a.set_yticklabels(labels[::-1] if False else labels, fontsize=6)
    ax_a.set_xlabel("AUROC (best direction)")
    ax_a.set_xlim(0.4, 0.85)
    ax_a.invert_yaxis()
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=PALETTE_CIRCUIT, label="CIRCUIT states"),
                       Patch(facecolor=PALETTE_BENCH, label="Benchmark signatures")]
    ax_a.legend(handles=legend_elements, loc="lower right", fontsize=6)

    # --- Panel B: Permutation null distribution for STATE_01 ---
    null = pd.read_csv(VAL / "melanoma_random_signature_null.tsv", sep="\t")
    state01_null = null[null["state"] == "CIRCUIT_STATE_01_score"]["random_auroc_best_direction"]

    ax_b.hist(state01_null, bins=40, color="#DDDDDD", edgecolor="#AAAAAA", linewidth=0.3, density=True)
    obs = pvals.loc[pvals["state"] == "CIRCUIT_STATE_01_score", "auroc_best_direction"].values[0]
    ax_b.axvline(obs, color=PALETTE_CIRCUIT, linewidth=1.5, linestyle="-")
    ax_b.text(obs + 0.01, ax_b.get_ylim()[1] * 0.85,
              f"Observed\nAUROC={obs:.3f}\nP=0.049", fontsize=6, color=PALETTE_CIRCUIT)
    ax_b.set_xlabel("Null AUROC (random gene sets)")
    ax_b.set_ylabel("Density")
    ax_b.set_title("STATE_01 vs 1000 random signatures", fontsize=7)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)

    fig.text(0.01, 0.97, "A", fontsize=10, fontweight="bold")
    fig.text(0.52, 0.97, "B", fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_auroc_permutation.pdf", bbox_inches="tight")
    fig.savefig(FIG / "fig3_auroc_permutation.png", bbox_inches="tight")
    plt.close()
    print("Figure 3 saved")


def fig4_cv_and_calibration():
    """Figure 4: CV comparison (A) and calibration curve (B)."""
    cv_circuit = pd.read_csv(VAL / "melanoma_cv_summary.tsv", sep="\t").iloc[0]
    cv_bench = pd.read_csv(VAL / "melanoma_cv_benchmark_summary.tsv", sep="\t").iloc[0]
    cv_comb = pd.read_csv(VAL / "melanoma_cv_combined_summary.tsv", sep="\t").iloc[0]
    cal = pd.read_csv(VAL / "melanoma_calibration.tsv", sep="\t")
    folds = pd.read_csv(VAL / "melanoma_cv_folds.tsv", sep="\t")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6, 2.8), gridspec_kw={"width_ratios": [1.1, 1]})

    # --- Panel A: Repeated CV AUROC comparison ---
    models = ["CIRCUIT\nstates", "Benchmark\nsignatures", "Combined"]
    means = [cv_circuit["cv_auroc_mean"], cv_bench["cv_auroc_mean"], cv_comb["cv_auroc_mean"]]
    stds = [cv_circuit["cv_auroc_std"], cv_bench["cv_auroc_std"], cv_comb["cv_auroc_std"]]
    ci_lo = [cv_circuit["cv_auroc_ci_low"], cv_bench["cv_auroc_ci_low"], cv_comb["cv_auroc_ci_low"]]
    ci_hi = [cv_circuit["cv_auroc_ci_high"], cv_bench["cv_auroc_ci_high"], cv_comb["cv_auroc_ci_high"]]
    bar_colors = [PALETTE_CIRCUIT, PALETTE_BENCH, PALETTE_COMBINED]

    x = np.arange(len(models))
    bars = ax_a.bar(x, means, width=0.55, color=bar_colors, edgecolor="none", zorder=2)
    for i in range(len(models)):
        ax_a.plot([x[i], x[i]], [ci_lo[i], ci_hi[i]], color="black", linewidth=1, zorder=3)

    ax_a.axhline(0.5, color="#999999", linestyle="--", linewidth=0.5, zorder=1)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(models, fontsize=7)
    ax_a.set_ylabel("CV AUROC")
    ax_a.set_ylim(0.35, 0.85)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)

    for i, m in enumerate(means):
        ax_a.text(x[i], m + stds[i] + 0.02, f"{m:.3f}", ha="center", fontsize=6.5, fontweight="bold")

    # --- Panel B: Calibration curve ---
    ax_b.plot([0, 1], [0, 1], "--", color="#999999", linewidth=0.5, label="Perfect calibration")
    ax_b.plot(cal["mean_predicted"], cal["fraction_positive"], "o-",
              color=PALETTE_CIRCUIT, linewidth=1.5, markersize=5, label="CIRCUIT model")
    ax_b.fill_between(cal["mean_predicted"],
                      cal["fraction_positive"] - 0.1,
                      cal["fraction_positive"] + 0.1,
                      alpha=0.15, color=PALETTE_CIRCUIT)
    ax_b.set_xlabel("Mean predicted probability")
    ax_b.set_ylabel("Observed fraction positive")
    ax_b.set_xlim(0, 1)
    ax_a.set_ylim(0.35, 0.85)
    ax_b.set_ylim(0, 0.8)
    ax_b.legend(fontsize=6, loc="upper left")
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)

    fig.text(0.01, 0.97, "A", fontsize=10, fontweight="bold")
    fig.text(0.52, 0.97, "B", fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "fig4_cv_calibration.pdf", bbox_inches="tight")
    fig.savefig(FIG / "fig4_cv_calibration.png", bbox_inches="tight")
    plt.close()
    print("Figure 4 saved")


def fig5_direction_stability():
    """Figure 5: Direction stability across melanoma validation cohorts."""
    stab = pd.read_csv(VAL / "melanoma_direction_stability.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(3.6, 2.8))

    circuit_stab = stab[stab["state"].str.startswith("CIRCUIT_STATE")].copy()
    circuit_stab = circuit_stab.sort_values("direction_agreement", ascending=True)
    labels = [s.replace("CIRCUIT_STATE_", "S").replace("_score", "") for s in circuit_stab["state"]]
    agreements = circuit_stab["direction_agreement"].values
    stable = circuit_stab["is_stable"].values

    bar_colors = [PALETTE_CIRCUIT if s else "#CCCCCC" for s in stable]
    ax.barh(range(len(labels)), agreements, color=bar_colors, edgecolor="none", height=0.6)
    ax.axvline(0.8, color="#999999", linestyle="--", linewidth=0.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Direction agreement across cohorts")
    ax.set_xlim(0, 1.1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for i, v in enumerate(agreements):
        ax.text(v + 0.02, i, f"{v:.0%}", va="center", fontsize=6)
    ax.text(0.82, -0.7, "80% threshold", fontsize=5.5, color="#999999")

    fig.tight_layout()
    fig.savefig(FIG / "fig5_direction_stability.pdf", bbox_inches="tight")
    fig.savefig(FIG / "fig5_direction_stability.png", bbox_inches="tight")
    plt.close()
    print("Figure 5 saved")


if __name__ == "__main__":
    fig1_k_scan()
    fig2_marker_heatmap()
    fig3_auroc_and_permutation()
    fig4_cv_and_calibration()
    fig5_direction_stability()
    print(f"\nAll figures saved to {FIG}")
