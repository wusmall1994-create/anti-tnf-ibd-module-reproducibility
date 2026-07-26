from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "results" / "bulk_validation" / "bulk_validation_enhanced_stats.csv"
FIG_DIR = ROOT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
    }
)

MODULE_ORDER = [
    "myeloid_inflammation",
    "stromal_fibroinflammatory",
    "epithelial_ifn_damage",
    "epithelial_barrier_maturity",
]
COLORS = {
    "myeloid_inflammation": "#4C78A8",
    "stromal_fibroinflammatory": "#B279A2",
    "epithelial_ifn_damage": "#E45756",
    "epithelial_barrier_maturity": "#54A24B",
}


def save_pub_py(fig, filename: Path, dpi: int = 600) -> None:
    fig.savefig(str(filename.with_suffix(".svg")), bbox_inches="tight")
    fig.savefig(str(filename.with_suffix(".png")), dpi=300, bbox_inches="tight")


def main() -> None:
    df = pd.read_csv(IN)
    df["module"] = pd.Categorical(df["module"], MODULE_ORDER, ordered=True)
    df["label"] = (
        df["cohort"]
        + " | "
        + df["comparison_group"]
        + " (R="
        + df["n_responder"].astype(str)
        + ", NR="
        + df["n_non_responder"].astype(str)
        + ")"
    )
    df = df.sort_values(["module", "cohort", "comparison_group"]).reset_index(drop=True)
    df["y"] = range(len(df), 0, -1)

    fig, ax = plt.subplots(figsize=(6.7, 4.2))
    for module, sub in df.groupby("module", observed=True):
        ax.scatter(
            sub["cliffs_delta_responder_vs_non"],
            sub["y"],
            s=sub["fdr_q_lt_0_10"].map({True: 52, False: 30}),
            color=COLORS[str(module)],
            edgecolor="0.2",
            linewidth=0.4,
            alpha=0.9,
            label=str(module),
        )
    ax.axvline(0, color="0.45", lw=0.8, ls="--")
    ax.set_yticks(df["y"])
    ax.set_yticklabels(df["label"])
    ax.set_xlabel("Cliff's delta: responder vs non-responder")
    ax.set_title(
        "Pretreatment module activity is lower in anti-TNF responders for inflammatory/damage programs",
        loc="left",
        fontweight="bold",
    )
    ax.text(
        -0.98,
        0.2,
        "Lower in responders",
        ha="left",
        va="bottom",
        fontsize=7,
        color="0.25",
    )
    ax.text(
        0.98,
        0.2,
        "Higher in responders",
        ha="right",
        va="bottom",
        fontsize=7,
        color="0.25",
    )
    for _, r in df.iterrows():
        ax.text(
            1.05,
            r["y"],
            f"q={r['q_value_fdr']:.3f}",
            va="center",
            ha="left",
            fontsize=6,
            color="0.25" if r["fdr_q_lt_0_10"] else "0.55",
        )
    ax.set_xlim(-1.1, 1.35)
    ax.set_ylim(0, len(df) + 1)
    ax.legend(title="", loc="lower left", bbox_to_anchor=(0, -0.25), ncol=2)
    fig.tight_layout()
    out = FIG_DIR / "fig_bulk_validation_forest"
    save_pub_py(fig, out)
    print(out)


if __name__ == "__main__":
    main()
