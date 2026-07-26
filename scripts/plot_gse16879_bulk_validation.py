from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "results" / "bulk_validation" / "GSE16879_bulk_module_scores_clean.csv"
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

MODULES = [
    "myeloid_inflammation",
    "stromal_fibroinflammatory",
    "epithelial_ifn_damage",
    "epithelial_barrier_maturity",
]
PALETTE = {"Responder": "#3B8BC2", "Non_responder": "#C95C4A"}


def save_pub_py(fig, filename: Path, dpi: int = 600) -> None:
    fig.savefig(str(filename.with_suffix(".svg")), bbox_inches="tight")
    fig.savefig(str(filename.with_suffix(".png")), dpi=300, bbox_inches="tight")


def main() -> None:
    df = pd.read_csv(IN)
    df = df[
        df["timepoint"].eq("Before")
        & df["response"].isin(["Responder", "Non_responder"])
        & df["tissue"].eq("Colon")
        & df["disease"].isin(["CD", "UC"])
    ].copy()
    long = df.melt(
        id_vars=["geo_accession", "response", "disease", "tissue"],
        value_vars=MODULES,
        var_name="module",
        value_name="score",
    )
    long["module"] = pd.Categorical(long["module"], MODULES, ordered=True)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=False)
    for ax, disease in zip(axes, ["CD", "UC"]):
        sub = long[long["disease"].eq(disease)]
        sns.boxplot(
            data=sub,
            x="module",
            y="score",
            hue="response",
            palette=PALETTE,
            fliersize=0,
            linewidth=0.7,
            ax=ax,
        )
        sns.stripplot(
            data=sub,
            x="module",
            y="score",
            hue="response",
            palette=PALETTE,
            dodge=True,
            size=2.3,
            alpha=0.7,
            linewidth=0.2,
            edgecolor="white",
            ax=ax,
            legend=False,
        )
        ax.set_title(f"Colon {disease}", loc="left", fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Bulk module score" if disease == "CD" else "")
        ax.tick_params(axis="x", rotation=35)
        for label in ax.get_xticklabels():
            label.set_ha("right")
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:2], labels[:2], title="", loc="best")

    fig.suptitle(
        "GSE16879 pretreatment modules associate with infliximab response",
        x=0.01,
        y=1.05,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.tight_layout()
    out = FIG_DIR / "fig_bulk_GSE16879_pretreatment_modules"
    save_pub_py(fig, out)
    print(out)


if __name__ == "__main__":
    main()
