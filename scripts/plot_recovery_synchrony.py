from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from pandas.errors import EmptyDataError
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
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


PALETTE = {"Remission": "#3B8BC2", "Non_Remission": "#C95C4A"}


def save_pub_py(fig, filename: Path, dpi: int = 600) -> None:
    fig.savefig(str(filename.with_suffix(".svg")), bbox_inches="tight")
    fig.savefig(str(filename.with_suffix(".png")), dpi=300, bbox_inches="tight")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--recovery-csv", required=True)
    parser.add_argument("--synchrony-csv", required=True)
    parser.add_argument("--label", default="first_pass")
    args = parser.parse_args()

    rec = pd.read_csv(args.recovery_csv)
    try:
        sync = pd.read_csv(args.synchrony_csv)
    except EmptyDataError:
        sync = pd.DataFrame()
    if rec.empty:
        raise SystemExit("Recovery CSV is empty; cannot plot.")

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.6), gridspec_kw={"width_ratios": [1.4, 1]})

    ax = axes[0]
    sns.boxplot(
        data=rec,
        x="module",
        y="recovery",
        hue="remission",
        palette=PALETTE,
        fliersize=0,
        linewidth=0.7,
        ax=ax,
    )
    sns.stripplot(
        data=rec,
        x="module",
        y="recovery",
        hue="remission",
        palette=PALETTE,
        dodge=True,
        size=2.5,
        alpha=0.65,
        linewidth=0.2,
        edgecolor="white",
        ax=ax,
        legend=False,
    )
    ax.axhline(0, color="0.45", lw=0.7, ls="--")
    ax.set_title("Compartment/module recovery toward healthy", loc="left", fontweight="bold")
    ax.set_ylabel("Recovery score")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=35)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], title="", loc="best")

    ax = axes[1]
    if sync.empty:
        ax.text(0.5, 0.5, "Synchrony not available", ha="center", va="center")
        ax.axis("off")
    else:
        sns.boxplot(
            data=sync,
            x="module",
            y="synchrony",
            hue="remission",
            palette=PALETTE,
            fliersize=0,
            linewidth=0.7,
            ax=ax,
        )
        sns.stripplot(
            data=sync,
            x="module",
            y="synchrony",
            hue="remission",
            palette=PALETTE,
            dodge=True,
            size=2.5,
            alpha=0.65,
            linewidth=0.2,
            edgecolor="white",
            ax=ax,
            legend=False,
        )
        ax.set_title("Cross-lineage synchrony", loc="left", fontweight="bold")
        ax.set_ylabel("Synchrony score")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=35)
        for label in ax.get_xticklabels():
            label.set_ha("right")
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:2], labels[:2], title="", loc="best")

    fig.suptitle(
        "Anti-TNF outcome as coordinated recovery toward a healthy reference",
        x=0.01,
        y=1.05,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.tight_layout()
    out = FIG_DIR / f"fig1_{args.label}_recovery_synchrony"
    save_pub_py(fig, out)
    print(out)


if __name__ == "__main__":
    main()
