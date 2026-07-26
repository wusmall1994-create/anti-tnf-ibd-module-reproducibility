from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "results" / "bulk_validation" / "GSE23597_bulk_module_scores_clean.csv"
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
    "epithelial_ifn_damage",
    "stromal_fibroinflammatory",
    "epithelial_barrier_maturity",
]
PALETTE = {"Responder": "#3B8BC2", "Non_responder": "#C95C4A"}


def save_pub_py(fig, filename: Path, dpi: int = 600) -> None:
    fig.savefig(str(filename.with_suffix(".svg")), bbox_inches="tight")
    fig.savefig(str(filename.with_suffix(".png")), dpi=300, bbox_inches="tight")


def main() -> None:
    df = pd.read_csv(IN)
    df = df[df["time"].eq("W0") & df["treatment_group"].eq("IFX")].copy()
    long = df.melt(
        id_vars=["geo_accession", "response_wk8"],
        value_vars=MODULES,
        var_name="module",
        value_name="score",
    )
    long["module"] = pd.Categorical(long["module"], MODULES, ordered=True)
    fig, ax = plt.subplots(figsize=(5.6, 2.7))
    sns.boxplot(
        data=long,
        x="module",
        y="score",
        hue="response_wk8",
        palette=PALETTE,
        fliersize=0,
        linewidth=0.7,
        ax=ax,
    )
    sns.stripplot(
        data=long,
        x="module",
        y="score",
        hue="response_wk8",
        palette=PALETTE,
        dodge=True,
        size=2.4,
        alpha=0.7,
        linewidth=0.2,
        edgecolor="white",
        ax=ax,
        legend=False,
    )
    ax.set_title("GSE23597 baseline IFX-treated UC biopsies", loc="left", fontweight="bold")
    ax.set_ylabel("Bulk module score")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=35)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], title="Week-8 response", loc="best")
    fig.suptitle(
        "Pretreatment inflammatory modules are lower in subsequent infliximab responders",
        x=0.01,
        y=1.05,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.tight_layout()
    out = FIG_DIR / "fig_bulk_GSE23597_baseline_modules"
    save_pub_py(fig, out)
    print(out)


if __name__ == "__main__":
    main()
