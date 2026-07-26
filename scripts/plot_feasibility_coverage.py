from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "results" / "feasibility" / "GSE282122_umap_sample_coverage.csv"
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


ORDER = ["CD_Remission", "CD_Non_Remission", "UC_Remission", "UC_Non_Remission", "No_outcome"]
PALETTE = {
    "CD_Remission": "#6BAED6",
    "CD_Non_Remission": "#9ECAE1",
    "UC_Remission": "#F4A582",
    "UC_Non_Remission": "#D6604D",
    "No_outcome": "#BDBDBD",
}


def save_pub_py(fig, filename: Path, dpi: int = 600) -> None:
    fig.savefig(str(filename.with_suffix(".svg")), bbox_inches="tight")
    fig.savefig(str(filename.with_suffix(".png")), dpi=300, bbox_inches="tight")


def main() -> None:
    df = pd.read_csv(IN)
    df["category"] = pd.Categorical(df["category"], ORDER, ordered=True)
    summary = (
        df.groupby("category", observed=True)
        .agg(
            cells=("n_cells_in_umap", "sum"),
            samples=("sample_id", "nunique"),
            patients=("patient", "nunique"),
        )
        .reset_index()
    )
    summary["cells_k"] = summary["cells"] / 1000

    fig, axes = plt.subplots(1, 2, figsize=(6.9, 2.4), gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    sns.barplot(
        data=summary,
        x="category",
        y="cells_k",
        hue="category",
        palette=PALETTE,
        ax=ax,
        legend=False,
        edgecolor="0.25",
        linewidth=0.4,
    )
    ax.set_ylabel("Cells in UMAP object (×10³)")
    ax.set_xlabel("")
    ax.set_title("Cell coverage by outcome group", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=35)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    for p, val in zip(ax.patches, summary["cells_k"]):
        ax.text(
            p.get_x() + p.get_width() / 2,
            p.get_height() + summary["cells_k"].max() * 0.015,
            f"{val:.0f}",
            ha="center",
            va="bottom",
            fontsize=6,
        )

    ax = axes[1]
    long = summary.melt(
        id_vars="category",
        value_vars=["samples", "patients"],
        var_name="metric",
        value_name="count",
    )
    sns.barplot(
        data=long,
        x="category",
        y="count",
        hue="metric",
        palette={"samples": "#9E9E9E", "patients": "#424242"},
        ax=ax,
        edgecolor="0.25",
        linewidth=0.4,
    )
    ax.set_ylabel("Count")
    ax.set_xlabel("")
    ax.set_title("Patient and sample coverage", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=35)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    ax.legend(title="", loc="upper right")

    fig.suptitle(
        "GSE282122 supports expression-level feasibility testing, with limited UC-remission patients",
        x=0.01,
        y=1.05,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.tight_layout()
    out = FIG_DIR / "figS1_feasibility_coverage"
    save_pub_py(fig, out)
    print(out)


if __name__ == "__main__":
    main()
