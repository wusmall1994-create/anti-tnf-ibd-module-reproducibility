"""Create a publication-style workflow schematic for the anti-TNF IBD project."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "figures"


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
    }
)


COLORS = {
    "navy": "#0B2545",
    "blue": "#2E74B5",
    "light_blue": "#E8F1FA",
    "teal": "#3A8F8F",
    "light_teal": "#E6F4F1",
    "orange": "#C9772A",
    "light_orange": "#FFF1E2",
    "gray": "#555555",
    "light_gray": "#F4F6F9",
    "border": "#B8C2CC",
    "green": "#3C7D3F",
    "red": "#9B1C1C",
}


def rounded_box(ax, xy, w, h, title, body, face, edge=None, title_color=None, body_fs=6.2, title_fs=7.3):
    x, y = xy
    edge = edge or COLORS["border"]
    title_color = title_color or COLORS["navy"]
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=0.9,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(
        x + 0.025,
        y + h - 0.06,
        title,
        ha="left",
        va="top",
        weight="bold",
        color=title_color,
        fontsize=title_fs,
    )
    ax.text(
        x + 0.025,
        y + h - 0.115,
        body,
        ha="left",
        va="top",
        color=COLORS["gray"],
        fontsize=body_fs,
        linespacing=1.25,
    )
    return patch


def arrow(ax, start, end, color="#6B7280", rad=0.0):
    arr = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=12,
        linewidth=1.1,
        color=color,
        connectionstyle=f"arc3,rad={rad}",
    )
    ax.add_patch(arr)


def small_tag(ax, xy, text, face, color):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        0.16,
        0.043,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=0,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(x + 0.08, y + 0.022, text, ha="center", va="center", fontsize=5.7, color=color, weight="bold")


def main():
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax.text(
        0.02,
        0.965,
        "Public-data workflow for interpretable anti-TNF response modules in IBD",
        ha="left",
        va="top",
        fontsize=10.5,
        weight="bold",
        color=COLORS["navy"],
    )
    ax.text(
        0.02,
        0.915,
        "A feasibility gate prevents overinterpretation of single-cell matrices; independent bulk cohorts test baseline mucosal modules.",
        ha="left",
        va="top",
        fontsize=6.8,
        color=COLORS["gray"],
    )

    # Main horizontal workflow
    boxes = [
        (
            (0.035, 0.67),
            0.19,
            0.18,
            "1. Public resources",
            "TAURUS / GSE282122\nGSE16879\nGSE23597\nGSE14580 overlap audit",
            COLORS["light_gray"],
            COLORS["border"],
        ),
        (
            (0.285, 0.67),
            0.19,
            0.18,
            "2. TAURUS gate",
            "216 biopsy libraries scored\n~1M cells audited\nNo lineage H5AD in GEO",
            COLORS["light_blue"],
            COLORS["blue"],
        ),
        (
            (0.535, 0.67),
            0.19,
            0.18,
            "3. Bulk validation",
            "Pretreatment modules\nGSE16879: CD + UC\nGSE23597: UC IFX baseline",
            COLORS["light_teal"],
            COLORS["teal"],
        ),
        (
            (0.775, 0.67),
            0.19,
            0.18,
            "4. Evidence synthesis",
            "Effect sizes + FDR\nOverlap audit\nBounded claims",
            COLORS["light_orange"],
            COLORS["orange"],
        ),
    ]

    for spec in boxes:
        rounded_box(ax, spec[0], spec[1], spec[2], spec[3], spec[4], spec[5], edge=spec[6])

    for x0, x1 in [(0.225, 0.285), (0.475, 0.535), (0.725, 0.775)]:
        arrow(ax, (x0, 0.76), (x1, 0.76))

    # Lower branch: negative gate and positive validation signal.
    rounded_box(
        ax,
        (0.135, 0.39),
        0.34,
        0.17,
        "Whole-biopsy recovery gate",
        "Recovery toward healthy was computable,\nbut did not robustly separate\nremission from non-remission.",
        "#FFFFFF",
        edge=COLORS["blue"],
        body_fs=6.0,
        title_fs=7.4,
    )
    small_tag(ax, (0.29, 0.505), "Neutral gate", "#FDECEC", COLORS["red"])

    rounded_box(
        ax,
        (0.535, 0.39),
        0.34,
        0.17,
        "Pretreatment module signal",
        "Responders showed lower inflammatory,\nepithelial-damage and fibroinflammatory\nactivity in public bulk cohorts.",
        "#FFFFFF",
        edge=COLORS["teal"],
        body_fs=6.0,
        title_fs=7.4,
    )
    small_tag(ax, (0.695, 0.505), "Primary evidence", "#E6F4EA", COLORS["green"])

    arrow(ax, (0.38, 0.67), (0.31, 0.56), color=COLORS["blue"], rad=0.08)
    arrow(ax, (0.63, 0.67), (0.705, 0.56), color=COLORS["teal"], rad=-0.08)

    # Bottom conclusion band.
    band = FancyBboxPatch(
        (0.08, 0.13),
        0.84,
        0.14,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=1.0,
        edgecolor=COLORS["navy"],
        facecolor="#F8FAFC",
    )
    ax.add_patch(band)
    ax.text(
        0.105,
        0.235,
        "Defensible conclusion",
        ha="left",
        va="top",
        fontsize=8,
        color=COLORS["navy"],
        weight="bold",
    )
    ax.text(
        0.105,
        0.19,
        "Pretreatment mucosal-state modules, especially lower myeloid inflammation, are associated with later anti-TNF response.",
        ha="left",
        va="top",
        fontsize=6.8,
        color=COLORS["gray"],
    )
    ax.text(
        0.105,
        0.155,
        "Boundary: current GEO matrices do not support lineage-specific single-cell localization without annotated H5AD data.",
        ha="left",
        va="top",
        fontsize=6.8,
        color=COLORS["gray"],
    )

    arrow(ax, (0.305, 0.39), (0.34, 0.27), color=COLORS["blue"], rad=0.08)
    arrow(ax, (0.705, 0.39), (0.66, 0.27), color=COLORS["teal"], rad=-0.08)

    ax.text(
        0.02,
        0.035,
        "IFX, infliximab; FDR, false-discovery rate; H5AD, annotated AnnData object.",
        ha="left",
        va="bottom",
        fontsize=6,
        color="#6B7280",
    )

    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / "fig1_public_data_workflow_schematic"
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(base)


if __name__ == "__main__":
    main()
