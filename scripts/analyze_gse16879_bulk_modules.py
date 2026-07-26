from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "bulk_validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MODULES = [
    "myeloid_inflammation",
    "epithelial_ifn_damage",
    "epithelial_barrier_maturity",
    "stromal_fibroinflammatory",
]


def get_field(row: pd.Series, prefix: str) -> str:
    pat = re.compile(rf"^{re.escape(prefix)}\s*:\s*(.*)$", re.I)
    vals = []
    for c in row.index:
        if not str(c).startswith("characteristics_ch1"):
            continue
        m = pat.match(str(row[c]).strip())
        if m:
            vals.append(m.group(1).strip())
    return vals[0] if vals else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(OUT_DIR / "GSE16879_bulk_module_scores.csv"))
    parser.add_argument("--label", default="GSE16879")
    args = parser.parse_args()

    in_path = Path(args.input)
    df = pd.read_csv(in_path)
    for field in [
        "tissue",
        "disease",
        "response to infliximab",
        "before or after first infliximab treatment",
    ]:
        df[field] = df.apply(lambda r, f=field: get_field(r, f), axis=1)

    df["response"] = df["response to infliximab"].replace(
        {"Yes": "Responder", "No": "Non_responder", "Not applicable": "Control"}
    )
    df["timepoint"] = df["before or after first infliximab treatment"].replace(
        {
            "Before first infliximab treatment": "Before",
            "After first infliximab treatment": "After",
            "Not applicable": "Control",
        }
    )

    clean_path = OUT_DIR / f"{args.label}_bulk_module_scores_clean.csv"
    df.to_csv(clean_path, index=False)

    tests = []
    for tissue in sorted(df["tissue"].dropna().unique()):
        for disease in sorted(df["disease"].dropna().unique()):
            sub = df[(df["tissue"].eq(tissue)) & (df["disease"].eq(disease))]
            if sub.empty:
                continue
            before = sub[sub["timepoint"].eq("Before")]
            for m in MODULES:
                r = before[before["response"].eq("Responder")][m].dropna()
                nr = before[before["response"].eq("Non_responder")][m].dropna()
                if len(r) >= 2 and len(nr) >= 2:
                    stat, p = stats.mannwhitneyu(r, nr, alternative="two-sided")
                else:
                    stat, p = np.nan, np.nan
                tests.append(
                    {
                        "comparison": "pretreatment_responder_vs_non_responder",
                        "tissue": tissue,
                        "disease": disease,
                        "module": m,
                        "n_responder": len(r),
                        "n_non_responder": len(nr),
                        "median_responder": float(np.median(r)) if len(r) else np.nan,
                        "median_non_responder": float(np.median(nr)) if len(nr) else np.nan,
                        "delta_responder_minus_non": float(np.median(r) - np.median(nr))
                        if len(r) and len(nr)
                        else np.nan,
                        "p_value": p,
                    }
                )

    tests_df = pd.DataFrame(tests)
    tests_path = OUT_DIR / f"{args.label}_bulk_module_tests.csv"
    tests_df.to_csv(tests_path, index=False)

    summary = {
        "input": str(in_path),
        "clean_scores": str(clean_path),
        "tests": str(tests_path),
        "n_samples": int(df.shape[0]),
        "response_counts": df["response"].value_counts(dropna=False).to_dict(),
        "timepoint_counts": df["timepoint"].value_counts(dropna=False).to_dict(),
        "disease_counts": df["disease"].value_counts(dropna=False).to_dict(),
        "tissue_counts": df["tissue"].value_counts(dropna=False).to_dict(),
    }
    summary_path = OUT_DIR / f"{args.label}_bulk_module_analysis_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(summary_path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
