"""TAURUS baseline cell-composition analysis: Remission vs Non-Remission.

Two layers:
1) Within-lineage cell-state proportions (state cells / lineage cells per sample).
2) Cross-lineage infiltration proxies from per-sample cell counts:
   myeloid/epithelial ratio and myeloid share of (myeloid+fibperi+epithelial).

Statistical unit = patient (sample-level proportions averaged per patient).
Baseline filter per author-team guidance: Treatment == 'Pre', Inflammation_score > 6.5,
CD/UC treated patients with outcome labels only.
"""

import json
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "taurus_zenodo"
OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "taurus_lineage"
OUT_DIR.mkdir(parents=True, exist_ok=True)

LINEAGES = {
    "myeloid": "myeloid_final.h5ad",
    "fibperi": "fibperi_final.h5ad",
    "epicolonic": "epicolonic_final.h5ad",
}
CUTOFF = 6.5


def cliffs_delta(x, y):
    if len(x) == 0 or len(y) == 0:
        return np.nan
    gt = sum((xi > y).sum() for xi in x)
    lt = sum((xi < y).sum() for xi in x)
    return (gt - lt) / (len(x) * len(y))


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(ranked) / (np.arange(len(ranked)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    return np.clip(q[np.argsort(order)], 0, 1)


def baseline_obs(h5ad_file):
    a = ad.read_h5ad(DATA_DIR / h5ad_file, backed="r")
    obs = a.obs[["sample_id", "Patient", "Disease", "Treatment",
                 "Inflammation_score", "Remission_status", "final_analysis"]].copy()
    a.file.close()
    obs = obs[obs["Treatment"] == "Pre"]
    obs = obs[obs["Inflammation_score"] > CUTOFF]
    obs = obs[obs["Remission_status"].isin(["Remission", "Non_Remission"])]
    obs = obs[obs["Disease"].isin(["CD", "UC"])]
    return obs


def compare_patient_level(df, value_col, keys=("Disease",)):
    """df: one row per patient x stratum with value_col; returns stats rows."""
    rows = []
    group_cols = list(keys) + ["stratum"]
    for key_vals, grp in df.groupby(group_cols, observed=True):
        if not isinstance(key_vals, tuple):
            key_vals = (key_vals,)
        rem = grp.loc[grp["Remission_status"] == "Remission", value_col].to_numpy()
        non = grp.loc[grp["Remission_status"] == "Non_Remission", value_col].to_numpy()
        if len(rem) < 3 or len(non) < 3:
            continue
        p = stats.mannwhitneyu(rem, non, alternative="two-sided").pvalue
        row = dict(zip(group_cols, key_vals))
        row.update({
            "n_remission": len(rem), "n_non_remission": len(non),
            "median_remission": np.median(rem), "median_non_remission": np.median(non),
            "median_diff": np.median(rem) - np.median(non),
            "mannwhitney_p": p, "cliffs_delta": cliffs_delta(rem, non),
        })
        rows.append(row)
    res = pd.DataFrame(rows)
    if len(res):
        res["fdr_q"] = bh_fdr(res["mannwhitney_p"].to_numpy())
        res = res.sort_values("mannwhitney_p")
    return res


def main():
    # ---------- layer 1: within-lineage state proportions ----------
    all_prop_rows = []
    sample_counts = {}  # lineage -> per-sample cell counts (for layer 2)
    for lineage, h5 in LINEAGES.items():
        obs = baseline_obs(h5)
        print(f"{lineage}: baseline treated cells used = {len(obs)}")
        # per-sample proportions of each state within this lineage
        per_sample = (obs.groupby(["Disease", "sample_id", "Patient", "Remission_status", "final_analysis"], observed=True)
                      .size().reset_index(name="n"))
        totals = per_sample.groupby("sample_id", observed=True)["n"].transform("sum")
        per_sample["proportion"] = per_sample["n"] / totals
        # keep states with >=2% mean proportion to avoid ultra-sparse strata
        keep = per_sample.groupby("final_analysis", observed=True)["proportion"].mean()
        keep = keep[keep >= 0.02].index
        per_sample = per_sample[per_sample["final_analysis"].isin(keep)]
        # patient level: mean proportion across the patient's samples
        pat = (per_sample.groupby(["Disease", "Patient", "Remission_status", "final_analysis"], observed=True)
               .agg(proportion=("proportion", "mean")).reset_index())
        pat = pat.rename(columns={"final_analysis": "stratum"})
        res = compare_patient_level(pat, "proportion")
        if len(res):
            res.insert(0, "lineage", lineage)
            res["analysis"] = "within_lineage_proportion"
            all_prop_rows.append(res)
        sample_counts[lineage] = (obs.groupby(["Disease", "sample_id", "Patient", "Remission_status"], observed=True)
                                  .size().reset_index(name="n_cells"))

    prop_res = pd.concat(all_prop_rows, ignore_index=True) if all_prop_rows else pd.DataFrame()
    if len(prop_res):
        prop_res["fdr_q_global"] = bh_fdr(prop_res["mannwhitney_p"].to_numpy())
        prop_res = prop_res.sort_values("mannwhitney_p")
    prop_res.to_csv(OUT_DIR / "composition_within_lineage_proportions.csv", index=False)

    # ---------- layer 2: cross-lineage infiltration proxies ----------
    counts = None
    for lineage, sc in sample_counts.items():
        sc = sc.rename(columns={"n_cells": f"n_{lineage}"})
        counts = sc if counts is None else counts.merge(
            sc, on=["Disease", "sample_id", "Patient", "Remission_status"], how="outer")
    counts = counts.fillna(0)
    counts["myeloid_over_epithelial"] = counts["n_myeloid"] / counts["n_epicolonic"].clip(lower=1)
    total3 = counts["n_myeloid"] + counts["n_fibperi"] + counts["n_epicolonic"]
    counts["myeloid_share_3lineages"] = counts["n_myeloid"] / total3.clip(lower=1)
    counts["stromal_share_3lineages"] = counts["n_fibperi"] / total3.clip(lower=1)

    proxy_rows = []
    for proxy in ["myeloid_over_epithelial", "myeloid_share_3lineages", "stromal_share_3lineages"]:
        pat = (counts.groupby(["Disease", "Patient", "Remission_status"], observed=True)
               .agg(**{proxy: (proxy, "mean")}).reset_index())
        pat["stratum"] = proxy
        res = compare_patient_level(pat, proxy)
        if len(res):
            res.insert(0, "lineage", "cross_lineage")
            res["analysis"] = "infiltration_proxy"
            proxy_rows.append(res)
    proxy_res = pd.concat(proxy_rows, ignore_index=True) if proxy_rows else pd.DataFrame()
    proxy_res.to_csv(OUT_DIR / "composition_infiltration_proxies.csv", index=False)

    print("\n=== infiltration proxies (patient level) ===")
    print(proxy_res.to_string(index=False) if len(proxy_res) else "none")
    print("\n=== within-lineage proportions: top 20 ===")
    print(prop_res.head(20).to_string(index=False) if len(prop_res) else "none")

    meta = {
        "baseline_filter": f"Treatment == 'Pre' and Inflammation_score > {CUTOFF}",
        "statistical_unit": "patient",
        "min_group_size": 3,
        "min_mean_proportion_for_state": 0.02,
        "note": "cross-lineage proxies use the three downloaded lineage objects only; "
                "myeloid share is relative to myeloid+fibperi+epithelial, not all cells.",
    }
    (OUT_DIR / "composition_run_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"\noutputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()
