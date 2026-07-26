from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "bulk_validation"
MODULES = [
    "myeloid_inflammation",
    "epithelial_ifn_damage",
    "epithelial_barrier_maturity",
    "stromal_fibroinflammatory",
]


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Cliff's delta for x versus y; positive means x tends to be larger."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    gt = 0
    lt = 0
    for val in x:
        gt += np.sum(val > y)
        lt += np.sum(val < y)
    return float((gt - lt) / (len(x) * len(y)))


def get_values_gse16879() -> pd.DataFrame:
    df = pd.read_csv(OUT / "GSE16879_bulk_module_scores_clean.csv")
    df = df[
        df["timepoint"].eq("Before")
        & df["response"].isin(["Responder", "Non_responder"])
        & df["tissue"].eq("Colon")
        & df["disease"].isin(["CD", "UC"])
    ].copy()
    rows = []
    for disease, sub in df.groupby("disease"):
        for module in MODULES:
            rows.append(
                {
                    "cohort": "GSE16879",
                    "comparison_group": f"Colon {disease}",
                    "independent_status": "independent",
                    "module": module,
                    "responder_values": sub[sub["response"].eq("Responder")][module].to_numpy(),
                    "non_responder_values": sub[sub["response"].eq("Non_responder")][module].to_numpy(),
                }
            )
    return pd.DataFrame(rows)


def get_values_gse23597() -> pd.DataFrame:
    df = pd.read_csv(OUT / "GSE23597_bulk_module_scores_clean.csv")
    df = df[
        df["time"].eq("W0")
        & df["treatment_group"].eq("IFX")
        & df["response_wk8"].isin(["Responder", "Non_responder"])
    ].copy()
    rows = []
    for module in MODULES:
        rows.append(
            {
                "cohort": "GSE23597",
                "comparison_group": "UC IFX baseline",
                "independent_status": "independent",
                "module": module,
                "responder_values": df[df["response_wk8"].eq("Responder")][module].to_numpy(),
                "non_responder_values": df[df["response_wk8"].eq("Non_responder")][module].to_numpy(),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    values = pd.concat([get_values_gse16879(), get_values_gse23597()], ignore_index=True)
    rows = []
    for _, r in values.iterrows():
        x = np.asarray(r["responder_values"], dtype=float)
        y = np.asarray(r["non_responder_values"], dtype=float)
        x = x[~np.isnan(x)]
        y = y[~np.isnan(y)]
        if len(x) >= 2 and len(y) >= 2:
            u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        else:
            u, p = np.nan, np.nan
        med_x = float(np.median(x)) if len(x) else np.nan
        med_y = float(np.median(y)) if len(y) else np.nan
        rows.append(
            {
                "cohort": r["cohort"],
                "comparison_group": r["comparison_group"],
                "independent_status": r["independent_status"],
                "module": r["module"],
                "n_responder": len(x),
                "n_non_responder": len(y),
                "median_responder": med_x,
                "median_non_responder": med_y,
                "delta_median_responder_minus_non": med_x - med_y if len(x) and len(y) else np.nan,
                "cliffs_delta_responder_vs_non": cliffs_delta(x, y),
                "mannwhitney_u": u,
                "p_value": p,
            }
        )

    res = pd.DataFrame(rows)
    mask = res["p_value"].notna()
    res["q_value_fdr"] = np.nan
    if mask.any():
        res.loc[mask, "q_value_fdr"] = multipletests(res.loc[mask, "p_value"], method="fdr_bh")[1]
    res["nominal_p_lt_0_05"] = res["p_value"] < 0.05
    res["fdr_q_lt_0_10"] = res["q_value_fdr"] < 0.10
    res["direction"] = np.where(
        res["delta_median_responder_minus_non"] > 0,
        "higher_in_responder",
        np.where(res["delta_median_responder_minus_non"] < 0, "lower_in_responder", "tie"),
    )
    out_path = OUT / "bulk_validation_enhanced_stats.csv"
    res.to_csv(out_path, index=False)

    # Module-level directional summary across independent comparisons.
    summary = (
        res.groupby("module", as_index=False)
        .agg(
            n_comparisons=("cohort", "count"),
            n_lower_in_responder=("direction", lambda s: int((s == "lower_in_responder").sum())),
            n_higher_in_responder=("direction", lambda s: int((s == "higher_in_responder").sum())),
            n_nominal=("nominal_p_lt_0_05", "sum"),
            n_fdr_q_lt_0_10=("fdr_q_lt_0_10", "sum"),
            median_cliffs_delta=("cliffs_delta_responder_vs_non", "median"),
            median_delta=("delta_median_responder_minus_non", "median"),
        )
    )
    summary_path = OUT / "bulk_validation_enhanced_module_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(res.to_string(index=False))
    print("\nMODULE SUMMARY")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
