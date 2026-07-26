from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from score_bulk_series_matrix_modules import (
    MODULES,
    collapse_to_gene,
    load_platform_mapping,
    parse_series_matrix,
)
from analyze_gse16879_bulk_modules import get_field as get_field_gse16879
from analyze_gse23597_bulk_modules import get_field as get_field_gse23597
from enhance_bulk_validation_stats import cliffs_delta


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "bulk_geo"
OUT = ROOT / "results" / "bulk_validation"


def zscore_gene_expression(gene_expr: pd.DataFrame) -> pd.DataFrame:
    """Z-score each gene across samples within one GEO series."""
    means = gene_expr.mean(axis=1)
    stds = gene_expr.std(axis=1, ddof=0).replace(0, np.nan)
    return gene_expr.sub(means, axis=0).div(stds, axis=0)


def score_modules_from_gene_expr(gene_expr: pd.DataFrame, drop_gene: str | None = None) -> pd.DataFrame:
    rows = pd.DataFrame({"geo_accession": gene_expr.columns})
    for module, genes in MODULES.items():
        use_genes = [g for g in genes if g != drop_gene and g in gene_expr.index]
        rows[module] = gene_expr.loc[use_genes].mean(axis=0).values if use_genes else np.nan
    return rows


def load_zscore_scores(gse: str) -> pd.DataFrame:
    matrix = DATA / f"{gse}_series_matrix.txt.gz"
    expr, meta, platforms = parse_series_matrix(matrix)
    platform = platforms[0] if len(platforms) == 1 else "GPL570"
    mapping = load_platform_mapping(platform, download=False)
    gene_expr = collapse_to_gene(expr, mapping)
    z_gene_expr = zscore_gene_expression(gene_expr)
    scores = score_modules_from_gene_expr(z_gene_expr)
    return scores.merge(meta, on="geo_accession", how="left"), z_gene_expr


def load_zscore_bundle(gse: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return metadata, z-scored gene expression and default module scores."""
    matrix = DATA / f"{gse}_series_matrix.txt.gz"
    expr, meta, platforms = parse_series_matrix(matrix)
    platform = platforms[0] if len(platforms) == 1 else "GPL570"
    mapping = load_platform_mapping(platform, download=False)
    gene_expr = collapse_to_gene(expr, mapping)
    z_gene_expr = zscore_gene_expression(gene_expr)
    scores = score_modules_from_gene_expr(z_gene_expr)
    return meta, z_gene_expr, scores


def prepare_gse16879(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for field in [
        "tissue",
        "disease",
        "response to infliximab",
        "before or after first infliximab treatment",
    ]:
        df[field] = df.apply(lambda r, f=field: get_field_gse16879(r, f), axis=1)
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
    return df


def prepare_gse23597(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for f in ["dosage_time_resp", "dose", "time", "subject", "wk8 response", "wk30 response", "tissue"]:
        df[f] = df.apply(lambda r, field=f: get_field_gse23597(r, field), axis=1)
    df["treatment_group"] = df["dose"].replace({"5mg/kg": "IFX", "10mg/kg": "IFX", "placebo": "Placebo"})
    df["response_wk8"] = df["wk8 response"].replace({"Yes": "Responder", "No": "Non_responder"})
    return df


def test_values(x: pd.Series, y: pd.Series) -> dict[str, float | int | bool]:
    x = pd.to_numeric(x, errors="coerce").dropna().to_numpy(dtype=float)
    y = pd.to_numeric(y, errors="coerce").dropna().to_numpy(dtype=float)
    if len(x) >= 2 and len(y) >= 2:
        _, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    else:
        p = np.nan
    med_x = float(np.median(x)) if len(x) else np.nan
    med_y = float(np.median(y)) if len(y) else np.nan
    delta = med_x - med_y if len(x) and len(y) else np.nan
    return {
        "n_responder": int(len(x)),
        "n_non_responder": int(len(y)),
        "median_responder": med_x,
        "median_non_responder": med_y,
        "delta_median_responder_minus_non": float(delta),
        "cliffs_delta_responder_vs_non": cliffs_delta(x, y),
        "p_value": float(p) if not np.isnan(p) else np.nan,
        "direction": "higher_in_responder" if delta > 0 else ("lower_in_responder" if delta < 0 else "tie"),
    }


def zscore_tests() -> pd.DataFrame:
    meta_16879, _, scores_16879 = load_zscore_bundle("GSE16879")
    gse16879 = prepare_gse16879(scores_16879.merge(meta_16879, on="geo_accession", how="left"))
    meta_23597, _, scores_23597 = load_zscore_bundle("GSE23597")
    gse23597 = prepare_gse23597(scores_23597.merge(meta_23597, on="geo_accession", how="left"))

    rows = []
    base = gse16879[
        gse16879["timepoint"].eq("Before")
        & gse16879["response"].isin(["Responder", "Non_responder"])
        & gse16879["tissue"].eq("Colon")
        & gse16879["disease"].isin(["CD", "UC"])
    ].copy()
    for disease, sub in base.groupby("disease"):
        for module in MODULES:
            r = sub[sub["response"].eq("Responder")][module]
            nr = sub[sub["response"].eq("Non_responder")][module]
            row = {
                "cohort": "GSE16879",
                "comparison_group": f"Colon {disease}",
                "module": module,
                "scoring": "gene_zscore_mean",
            }
            row.update(test_values(r, nr))
            rows.append(row)

    base = gse23597[
        gse23597["time"].eq("W0")
        & gse23597["treatment_group"].eq("IFX")
        & gse23597["response_wk8"].isin(["Responder", "Non_responder"])
    ].copy()
    for module in MODULES:
        r = base[base["response_wk8"].eq("Responder")][module]
        nr = base[base["response_wk8"].eq("Non_responder")][module]
        row = {
            "cohort": "GSE23597",
            "comparison_group": "UC IFX baseline",
            "module": module,
            "scoring": "gene_zscore_mean",
        }
        row.update(test_values(r, nr))
        rows.append(row)

    out = pd.DataFrame(rows)
    mask = out["p_value"].notna()
    out["q_value_fdr"] = np.nan
    out.loc[mask, "q_value_fdr"] = multipletests(out.loc[mask, "p_value"], method="fdr_bh")[1]
    return out


def leave_one_gene_out_direction_audit() -> pd.DataFrame:
    """Check whether responder-lower/higher direction depends on one gene."""
    rows = []
    primary = pd.read_csv(OUT / "bulk_validation_enhanced_stats.csv")
    primary_key = primary.set_index(["cohort", "comparison_group", "module"])["direction"].to_dict()
    bundles = {
        "GSE16879": (*load_zscore_bundle("GSE16879")[:2], prepare_gse16879),
        "GSE23597": (*load_zscore_bundle("GSE23597")[:2], prepare_gse23597),
    }
    for gse, (meta, z_expr, prepare) in bundles.items():
        for module, genes in MODULES.items():
            present = [g for g in genes if g in z_expr.index]
            for drop_gene in present:
                scores = score_modules_from_gene_expr(z_expr, drop_gene=drop_gene)
                df = prepare(scores.merge(meta, on="geo_accession", how="left"))
                if gse == "GSE16879":
                    base = df[
                        df["timepoint"].eq("Before")
                        & df["response"].isin(["Responder", "Non_responder"])
                        & df["tissue"].eq("Colon")
                        & df["disease"].isin(["CD", "UC"])
                    ].copy()
                    groups = [(f"Colon {d}", sub, "response") for d, sub in base.groupby("disease")]
                else:
                    base = df[
                        df["time"].eq("W0")
                        & df["treatment_group"].eq("IFX")
                        & df["response_wk8"].isin(["Responder", "Non_responder"])
                    ].copy()
                    groups = [("UC IFX baseline", base, "response_wk8")]
                for comparison_group, sub, response_col in groups:
                    r = sub[sub[response_col].eq("Responder")][module]
                    nr = sub[sub[response_col].eq("Non_responder")][module]
                    tested = test_values(r, nr)
                    primary_direction = primary_key.get((gse, comparison_group, module), "")
                    rows.append(
                        {
                            "cohort": gse,
                            "comparison_group": comparison_group,
                            "module": module,
                            "dropped_gene": drop_gene,
                            "primary_direction": primary_direction,
                            "loo_direction": tested["direction"],
                            "direction_matches_primary": tested["direction"] == primary_direction,
                            "cliffs_delta_responder_vs_non": tested["cliffs_delta_responder_vs_non"],
                            "p_value": tested["p_value"],
                        }
                    )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    z = zscore_tests()
    z_path = OUT / "bulk_module_zscore_sensitivity_stats.csv"
    z.to_csv(z_path, index=False)

    loo = leave_one_gene_out_direction_audit()
    loo_path = OUT / "bulk_module_leave_one_gene_out_direction_audit.csv"
    loo.to_csv(loo_path, index=False)

    module_rows = []
    for module, sub in loo.groupby("module"):
        module_rows.append(
            {
                "module": module,
                "n_leave_one_gene_tests": int(sub.shape[0]),
                "n_direction_mismatches": int((~sub["direction_matches_primary"]).sum()),
                "direction_match_rate": float(sub["direction_matches_primary"].mean()),
            }
        )
    summary = {
        "zscore_stats": str(z_path.relative_to(ROOT)),
        "leave_one_gene_out_audit": str(loo_path.relative_to(ROOT)),
        "zscore_direction_summary": z.groupby("module")["direction"].value_counts().unstack(fill_value=0).to_dict(),
        "leave_one_gene_out_summary": module_rows,
    }
    summary_path = OUT / "bulk_module_score_sensitivity_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
