from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "recovery"
OUT_DIR.mkdir(parents=True, exist_ok=True)


MODULE_COLUMNS = [
    "myeloid_inflammation",
    "epithelial_ifn_damage",
    "epithelial_barrier_maturity",
    "stromal_fibroinflammatory",
]


def zscore_against_healthy(df: pd.DataFrame, modules: list[str]) -> pd.DataFrame:
    out = df.copy()
    healthy = out[out["disease"].eq("Healthy")]
    for m in modules:
        mu = healthy[m].mean(skipna=True)
        sd = healthy[m].std(skipna=True)
        if pd.isna(sd) or sd == 0:
            sd = out[m].std(skipna=True)
        if pd.isna(sd) or sd == 0:
            out[m + "_z"] = np.nan
        else:
            out[m + "_z"] = (out[m] - mu) / sd
    return out


def winsorize_series(s: pd.Series, lower: float = 0.05, upper: float = 0.95) -> pd.Series:
    if s.dropna().empty:
        return s
    lo = s.quantile(lower)
    hi = s.quantile(upper)
    return s.clip(lo, hi)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--score-csv",
        action="append",
        required=True,
        help="One or more *_sample_module_scores.csv files.",
    )
    parser.add_argument("--label", default="combined")
    parser.add_argument(
        "--pair-by-site",
        action="store_true",
        help="Compute pre/post recovery within patient and biopsy site before testing.",
    )
    args = parser.parse_args()

    frames = []
    for p in args.score_csv:
        path = Path(p)
        label = path.name.replace("_sample_module_scores.csv", "")
        df = pd.read_csv(path)
        df["lineage"] = label
        frames.append(df)
    all_scores = pd.concat(frames, ignore_index=True, sort=False)

    modules = [m for m in MODULE_COLUMNS if m in all_scores.columns]
    all_scores = zscore_against_healthy(all_scores, modules)

    rows = []
    group_cols = ["patient", "lineage", "Category"]
    if args.pair_by_site and "site" in all_scores.columns:
        group_cols.append("site")

    for keys, g in all_scores.groupby(group_cols, dropna=False):
        key_map = dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,)))
        patient = key_map.get("patient")
        lineage = key_map.get("lineage")
        category = key_map.get("Category")
        site = key_map.get("site", "")
        if not isinstance(patient, str) or not patient or patient == "nan":
            continue
        if str(category) == "nan" or category == "No_outcome":
            continue
        pre = g[g["treatment"].eq("Pre")]
        post = g[g["treatment"].eq("Post")]
        if pre.empty or post.empty:
            continue

        # If multiple sites exist, use average as first-pass patient-level summary.
        pre_mean = pre[[m + "_z" for m in modules]].mean()
        post_mean = post[[m + "_z" for m in modules]].mean()
        for m in modules:
            pre_dist = abs(pre_mean[m + "_z"])
            post_dist = abs(post_mean[m + "_z"])
            recovery = np.nan if pre_dist == 0 or pd.isna(pre_dist) else 1 - post_dist / pre_dist
            delta_abs_distance = pre_dist - post_dist
            rows.append(
                {
                    "patient": patient,
                    "lineage": lineage,
                    "category": category,
                    "site": site,
                    "disease": str(category).split("_", 1)[0],
                    "remission": "Non_Remission"
                    if "Non_Remission" in str(category)
                    else "Remission",
                    "module": m,
                    "pre_z": pre_mean[m + "_z"],
                    "post_z": post_mean[m + "_z"],
                    "pre_abs_distance_to_healthy": pre_dist,
                    "post_abs_distance_to_healthy": post_dist,
                    "recovery": recovery,
                    "delta_abs_distance_to_healthy": delta_abs_distance,
                    "n_pre_samples": pre.shape[0],
                    "n_post_samples": post.shape[0],
                }
            )

    recovery_df = pd.DataFrame(rows)
    if not recovery_df.empty:
        recovery_df["recovery_winsorized"] = recovery_df.groupby(
            ["disease", "module"], group_keys=False
        )["recovery"].apply(winsorize_series)
        # Testing should remain patient-level. If site-aware pairing produces multiple
        # rows per patient/module/lineage, average them before group comparisons.
        test_recovery_df = (
            recovery_df.groupby(
                ["patient", "lineage", "category", "disease", "remission", "module"],
                dropna=False,
                as_index=False,
            )
            .agg(
                pre_z=("pre_z", "mean"),
                post_z=("post_z", "mean"),
                pre_abs_distance_to_healthy=("pre_abs_distance_to_healthy", "mean"),
                post_abs_distance_to_healthy=("post_abs_distance_to_healthy", "mean"),
                recovery=("recovery", "mean"),
                recovery_winsorized=("recovery_winsorized", "mean"),
                delta_abs_distance_to_healthy=("delta_abs_distance_to_healthy", "mean"),
                n_pre_samples=("n_pre_samples", "sum"),
                n_post_samples=("n_post_samples", "sum"),
            )
        )
    else:
        test_recovery_df = recovery_df
    rec_path = OUT_DIR / f"{args.label}_patient_lineage_module_recovery.csv"
    recovery_df.to_csv(rec_path, index=False)
    test_rec_path = OUT_DIR / f"{args.label}_patient_lineage_module_recovery_for_tests.csv"
    test_recovery_df.to_csv(test_rec_path, index=False)

    sync_rows = []
    if not test_recovery_df.empty:
        for (patient, category, module), g in test_recovery_df.groupby(["patient", "category", "module"]):
            value_for_sync = "recovery_winsorized" if "recovery_winsorized" in g.columns else "recovery"
            vals = g[value_for_sync].dropna().to_numpy()
            if vals.size < 2:
                continue
            sync_rows.append(
                {
                    "patient": patient,
                    "category": category,
                    "disease": str(category).split("_", 1)[0],
                    "remission": "Non_Remission"
                    if "Non_Remission" in str(category)
                    else "Remission",
                    "module": module,
                    "mean_recovery": float(np.mean(vals)),
                    "sd_recovery_across_lineages": float(np.std(vals, ddof=1))
                    if vals.size > 1
                    else np.nan,
                    "synchrony": -float(np.std(vals, ddof=1)) if vals.size > 1 else np.nan,
                    "n_lineages": int(vals.size),
                }
            )
    sync_df = pd.DataFrame(sync_rows)
    sync_path = OUT_DIR / f"{args.label}_patient_synchrony.csv"
    sync_df.to_csv(sync_path, index=False)

    tests = []
    for df_name, df, value_col in [
        ("recovery", test_recovery_df, "recovery"),
        ("recovery_winsorized", test_recovery_df, "recovery_winsorized"),
        ("delta_abs_distance", test_recovery_df, "delta_abs_distance_to_healthy"),
        ("synchrony", sync_df, "synchrony"),
    ]:
        if df.empty:
            continue
        for keys, g in df.groupby(["disease", "module"]):
            remission = g[g["remission"].eq("Remission")][value_col].dropna()
            non = g[g["remission"].eq("Non_Remission")][value_col].dropna()
            if remission.size < 2 or non.size < 2:
                p = np.nan
                stat = np.nan
            else:
                stat, p = stats.mannwhitneyu(remission, non, alternative="two-sided")
            tests.append(
                {
                    "analysis": df_name,
                    "disease": keys[0],
                    "module": keys[1],
                    "n_remission": int(remission.size),
                    "n_non_remission": int(non.size),
                    "median_remission": float(np.median(remission)) if remission.size else np.nan,
                    "median_non_remission": float(np.median(non)) if non.size else np.nan,
                    "mannwhitney_u": stat,
                    "p_value": p,
                }
            )

    tests_df = pd.DataFrame(tests)
    tests_path = OUT_DIR / f"{args.label}_recovery_synchrony_tests.csv"
    tests_df.to_csv(tests_path, index=False)

    audit = {
        "input_score_csvs": args.score_csv,
        "modules": modules,
        "outputs": {
            "recovery": str(rec_path),
            "recovery_for_tests": str(test_rec_path),
            "synchrony": str(sync_path),
            "tests": str(tests_path),
        },
        "n_recovery_rows": int(recovery_df.shape[0]),
        "n_test_recovery_rows": int(test_recovery_df.shape[0]),
        "n_synchrony_rows": int(sync_df.shape[0]),
        "pair_by_site": args.pair_by_site,
    }
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
