from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "recovery"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def direction(row: pd.Series) -> str:
    mr = row.get("median_remission")
    mn = row.get("median_non_remission")
    if pd.isna(mr) or pd.isna(mn):
        return "insufficient"
    if mr > mn:
        return "remission_higher"
    if mr < mn:
        return "non_remission_higher"
    return "tie"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", required=True)
    parser.add_argument("--label", required=True)
    args = parser.parse_args()

    tests = pd.read_csv(args.tests)
    if tests.empty:
        raise SystemExit("Empty tests table")
    tests["direction"] = tests.apply(direction, axis=1)
    tests["effect_delta_median"] = tests["median_remission"] - tests["median_non_remission"]
    tests["nominal"] = tests["p_value"].lt(0.05)

    rows = []
    for analysis, g in tests.groupby("analysis"):
        informative = g[g["direction"].ne("insufficient")]
        rows.append(
            {
                "analysis": analysis,
                "n_comparisons": int(g.shape[0]),
                "n_informative": int(informative.shape[0]),
                "n_remission_higher": int((informative["direction"] == "remission_higher").sum()),
                "n_non_remission_higher": int(
                    (informative["direction"] == "non_remission_higher").sum()
                ),
                "n_nominal_p_lt_0_05": int(g["nominal"].sum()),
                "median_effect_delta": float(informative["effect_delta_median"].median())
                if not informative.empty
                else np.nan,
            }
        )
    summary = pd.DataFrame(rows)

    # Conservative text decision: this is a gate, not a manuscript conclusion.
    recovery = summary[summary["analysis"].isin(["recovery_winsorized", "delta_abs_distance"])]
    if recovery.empty:
        decision = "NO_GO_NO_RECOVERY_ANALYSIS"
        rationale = "No recovery or delta-distance analysis rows were available."
    else:
        rem_higher = int(recovery["n_remission_higher"].sum())
        non_higher = int(recovery["n_non_remission_higher"].sum())
        nominal = int(recovery["n_nominal_p_lt_0_05"].sum())
        if rem_higher >= non_higher and (rem_higher >= 4 or nominal >= 1):
            decision = "PROVISIONAL_GO_TO_ANNOTATED_LINEAGE_GATE"
            rationale = (
                "Recovery-related analyses show at least as many remission-favouring "
                "directions as non-remission-favouring directions. Treat as provisional "
                "until full sample and sensitivity analyses are complete."
            )
        elif rem_higher > 0:
            decision = "WEAK_SIGNAL_NEEDS_FULL_AND_SENSITIVITY_CHECKS"
            rationale = (
                "Some remission-favouring directions exist, but the pattern is not strong "
                "enough for a manuscript claim."
            )
        else:
            decision = "PIVOT_OR_REDEFINE"
            rationale = (
                "The predefined whole-biopsy recovery scores do not currently show a "
                "remission-favouring pattern."
            )

    out_summary = OUT_DIR / f"{args.label}_gate_summary.csv"
    out_json = OUT_DIR / f"{args.label}_gate_decision.json"
    summary.to_csv(out_summary, index=False)
    out_json.write_text(
        json.dumps(
            {
                "label": args.label,
                "tests": args.tests,
                "decision": decision,
                "rationale": rationale,
                "summary_csv": str(out_summary),
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    print(out_json.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
