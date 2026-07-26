from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "bulk_validation"


def load_tests(path: Path, cohort: str, independent: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["cohort"] = cohort
    df["independent_status"] = independent
    df["direction_supports_response"] = df["delta_responder_minus_non"].apply(
        lambda x: "higher_in_responder" if x > 0 else ("lower_in_responder" if x < 0 else "tie")
    )
    return df


def main() -> None:
    frames = []
    p = OUT / "GSE16879_bulk_module_tests.csv"
    if p.exists():
        frames.append(load_tests(p, "GSE16879", "primary_bulk_validation"))
    p = OUT / "GSE14580_bulk_module_tests.csv"
    if p.exists():
        frames.append(load_tests(p, "GSE14580", "overlaps_GSE16879_UC_subset"))
    p = OUT / "GSE23597_bulk_module_tests.csv"
    if p.exists():
        frames.append(load_tests(p, "GSE23597", "independent_bulk_validation"))
    if not frames:
        raise SystemExit("No bulk test files found.")
    all_tests = pd.concat(frames, ignore_index=True, sort=False)
    out_all = OUT / "bulk_validation_all_module_tests.csv"
    all_tests.to_csv(out_all, index=False)

    informative = all_tests[
        all_tests["n_responder"].fillna(0).ge(2) & all_tests["n_non_responder"].fillna(0).ge(2)
    ].copy()
    summary = (
        informative.groupby(["cohort", "independent_status", "module"], dropna=False)
        .agg(
            n_responder=("n_responder", "max"),
            n_non_responder=("n_non_responder", "max"),
            median_delta=("delta_responder_minus_non", "median"),
            min_p=("p_value", "min"),
        )
        .reset_index()
    )
    summary["nominal_p_lt_0_05"] = summary["min_p"] < 0.05
    out_summary = OUT / "bulk_validation_summary_by_module.csv"
    summary.to_csv(out_summary, index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
