from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "results" / "bulk_validation" / "GSE23597_bulk_module_scores.csv"
OUT_DIR = ROOT / "results" / "bulk_validation"

MODULES = [
    "myeloid_inflammation",
    "epithelial_ifn_damage",
    "epithelial_barrier_maturity",
    "stromal_fibroinflammatory",
]


def get_field(row: pd.Series, prefix: str) -> str:
    pat = re.compile(rf"^{re.escape(prefix)}\s*:\s*(.*)$", re.I)
    for c in row.index:
        if not str(c).startswith("characteristics_ch1"):
            continue
        m = pat.match(str(row[c]).strip())
        if m:
            return m.group(1).strip()
    return ""


def main() -> None:
    df = pd.read_csv(IN)
    for f in ["dosage_time_resp", "dose", "time", "subject", "wk8 response", "wk30 response", "tissue"]:
        df[f] = df.apply(lambda r, field=f: get_field(r, field), axis=1)
    df["treatment_group"] = df["dose"].replace({"5mg/kg": "IFX", "10mg/kg": "IFX", "placebo": "Placebo"})
    df["response_wk8"] = df["wk8 response"].replace({"Yes": "Responder", "No": "Non_responder"})

    clean = OUT_DIR / "GSE23597_bulk_module_scores_clean.csv"
    df.to_csv(clean, index=False)

    tests = []
    baseline = df[df["time"].eq("W0") & df["treatment_group"].eq("IFX")]
    for m in MODULES:
        r = baseline[baseline["response_wk8"].eq("Responder")][m].dropna()
        nr = baseline[baseline["response_wk8"].eq("Non_responder")][m].dropna()
        if len(r) >= 2 and len(nr) >= 2:
            stat, p = stats.mannwhitneyu(r, nr, alternative="two-sided")
        else:
            stat, p = np.nan, np.nan
        tests.append(
            {
                "comparison": "IFX_W0_wk8_responder_vs_non_responder",
                "module": m,
                "n_responder": len(r),
                "n_non_responder": len(nr),
                "median_responder": float(np.median(r)) if len(r) else np.nan,
                "median_non_responder": float(np.median(nr)) if len(nr) else np.nan,
                "delta_responder_minus_non": float(np.median(r) - np.median(nr)) if len(r) and len(nr) else np.nan,
                "p_value": p,
            }
        )

    tests_df = pd.DataFrame(tests)
    tests_path = OUT_DIR / "GSE23597_bulk_module_tests.csv"
    tests_df.to_csv(tests_path, index=False)
    summary = {
        "input": str(IN),
        "clean_scores": str(clean),
        "tests": str(tests_path),
        "n_samples": int(df.shape[0]),
        "dose_counts": df["dose"].value_counts(dropna=False).to_dict(),
        "time_counts": df["time"].value_counts(dropna=False).to_dict(),
        "wk8_response_counts": df["response_wk8"].value_counts(dropna=False).to_dict(),
        "baseline_ifx_wk8_response_counts": baseline["response_wk8"].value_counts(dropna=False).to_dict(),
    }
    summary_path = OUT_DIR / "GSE23597_bulk_module_analysis_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(summary_path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
