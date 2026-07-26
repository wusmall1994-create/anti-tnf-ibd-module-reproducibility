"""TAURUS lineage-level module localization (myeloid lineage).

Scores the prespecified myeloid-inflammation module within the annotated
TAURUS myeloid object and compares Remission vs Non-Remission at the
sample-by-cell-state level (baseline samples only).

Scoring convention matches the manuscript whole-biopsy scoring:
counts normalized to 10,000 per cell, log1p, mean over available genes,
then averaged across cells within each sample x cell-state stratum.

Author-team guidance (Zenodo record note): baseline analyses use samples with
Inflammation_score > 6.5. We apply this filter and also record the
unfiltered counts for the audit trail.
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

MODULES = {
    "myeloid_inflammation": [
        "IL1B", "TNF", "CXCL8", "SPP1", "OSM", "IL1RN", "CCL3", "CCL4", "NFKBIA",
    ],
}

H5AD = "myeloid_final.h5ad"
INFLAMMATION_SCORE_CUTOFF = 6.5


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Cliff's delta for x vs y (>0 means x higher)."""
    n_x, n_y = len(x), len(y)
    if n_x == 0 or n_y == 0:
        return np.nan
    gt = sum((xi > y).sum() for xi in x)
    lt = sum((xi < y).sum() for xi in x)
    return (gt - lt) / (n_x * n_y)


def main():
    adata = ad.read_h5ad(DATA_DIR / H5AD, backed="r")
    genes = [g for g in MODULES["myeloid_inflammation"] if g in adata.var_names]
    missing = [g for g in MODULES["myeloid_inflammation"] if g not in adata.var_names]
    print(f"module genes found: {len(genes)}/9; missing: {missing}")

    sub = adata[:, genes].to_memory()
    obs = sub.obs.copy()

    # Normalize to 10k counts per cell, log1p, mean over module genes.
    counts = sub.X
    lib = np.asarray(counts.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    normed = counts.multiply(10000.0 / lib[:, None])
    normed.data = np.log1p(normed.data)
    obs["module_score"] = np.asarray(normed.mean(axis=1)).ravel()

    print("\ncell-level Treatment x Remission_status (all cells):")
    print(pd.crosstab(obs["Treatment"], obs["Remission_status"]))

    # Baseline filter per author-team guidance.
    pre = obs[obs["Treatment"] == "Pre"].copy()
    print(f"\nbaseline cells: {len(pre)}; "
          f"baseline cells with Inflammation_score > {INFLAMMATION_SCORE_CUTOFF}: "
          f"{(pre['Inflammation_score'] > INFLAMMATION_SCORE_CUTOFF).sum()}")
    pre = pre[pre["Inflammation_score"] > INFLAMMATION_SCORE_CUTOFF]

    # Keep treated patients with an outcome label.
    pre = pre[pre["Remission_status"].isin(["Remission", "Non_Remission"])]
    pre = pre[pre["Disease"].isin(["CD", "UC"])]
    print(f"baseline high-inflammation treated cells used: {len(pre)}")

    # Aggregate in two steps to avoid pseudoreplication:
    #   1) mean module score per sample x cell state (with cell support);
    #   2) patient x cell state = cell-weighted mean across that patient's samples.
    # The statistical unit for group comparisons is the PATIENT, not the sample.
    sample_agg = (
        pre.groupby(["Disease", "sample_id", "Patient", "Remission_status", "final_analysis"], observed=True)
        .agg(module_score=("module_score", "mean"), n_cells=("module_score", "size"))
        .reset_index()
    )
    sample_agg = sample_agg[sample_agg["n_cells"] >= 10]  # minimal cell support per stratum

    def _wmean(g):
        return np.average(g["module_score"], weights=g["n_cells"])

    agg = (
        sample_agg.groupby(["Disease", "Patient", "Remission_status", "final_analysis"], observed=True)
        .apply(lambda g: pd.Series({
            "module_score": _wmean(g),
            "n_cells": g["n_cells"].sum(),
            "n_samples": g["sample_id"].nunique(),
        }), include_groups=False)
        .reset_index()
    )

    # Patient-level group sizes per disease.
    samples = (
        agg[["Disease", "Patient", "Remission_status"]]
        .drop_duplicates()
        .groupby(["Disease", "Remission_status"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    print("\nPATIENT counts per disease stratum (statistical unit):")
    print(samples)

    # Per disease x cell-state comparison: Remission vs Non_Remission.
    rows = []
    for (dz, state), grp in agg.groupby(["Disease", "final_analysis"], observed=True):
        rem = grp.loc[grp["Remission_status"] == "Remission", "module_score"].to_numpy()
        non = grp.loc[grp["Remission_status"] == "Non_Remission", "module_score"].to_numpy()
        if len(rem) < 3 or len(non) < 3:
            continue
        p = stats.mannwhitneyu(rem, non, alternative="two-sided").pvalue
        rows.append({
            "disease": dz,
            "cell_state": state,
            "n_remission": len(rem),
            "n_non_remission": len(non),
            "median_remission": np.median(rem),
            "median_non_remission": np.median(non),
            "median_diff_rem_minus_non": np.median(rem) - np.median(non),
            "mannwhitney_p": p,
            "cliffs_delta": cliffs_delta(rem, non),
        })
    res = pd.DataFrame(rows)
    if len(res):
        res = res.sort_values("mannwhitney_p")
        # BH-FDR within this lineage localization family.
        pvals = res["mannwhitney_p"].to_numpy()
        order = np.argsort(pvals)
        ranked = pvals[order]
        q = ranked * len(ranked) / (np.arange(len(ranked)) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        res["fdr_q"] = np.clip(q[np.argsort(order)], 0, 1)
    res.to_csv(OUT_DIR / "myeloid_lineage_module_by_state.csv", index=False)
    agg.to_csv(OUT_DIR / "myeloid_sample_state_module_scores.csv", index=False)

    meta = {
        "h5ad": H5AD,
        "module_genes_used": genes,
        "module_genes_missing": missing,
        "baseline_filter": f"Treatment == 'Pre' and Inflammation_score > {INFLAMMATION_SCORE_CUTOFF}",
        "min_cells_per_sample_state": 10,
        "n_result_rows": int(len(res)),
    }
    (OUT_DIR / "myeloid_lineage_run_meta.json").write_text(json.dumps(meta, indent=2))

    print("\n=== top results (sorted by p) ===")
    if len(res):
        print(res.head(20).to_string(index=False))
    else:
        print("no strata met the n>=3 per-group threshold")
    print(f"\noutputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()
