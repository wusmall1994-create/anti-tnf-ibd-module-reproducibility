from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, hypergeom, norm, t
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_bulk_series_matrix_modules import (  # noqa: E402
    MODULES,
    collapse_to_gene,
    load_platform_mapping,
    parse_series_matrix,
)


OUT = ROOT / "results" / "enhancement_gate_20260805"
OUT.mkdir(parents=True, exist_ok=True)

PUBLISHED_SIGNATURES: dict[str, list[str]] = {
    # Arijs et al., Gut 2009, PMID 19700435: top five genes in the combined cohorts.
    "Arijs_2009_top5": ["TNFRSF11B", "STC1", "PTGS2", "IL13RA2", "IL11"],
    # West et al., Nat Med 2017, Fig. 2a: OSM plus 21 consistently associated cytokines/chemokines.
    "West_2017_OSM22": [
        "OSM", "IL1B", "IL1A", "IL6", "IL11", "CSF2", "CSF3", "IFNG", "IL17A", "IL22",
        "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CCL2", "CCL3", "CCL4",
        "CXCL9", "CXCL10", "CXCL11",
    ],
    # West et al. also reported the two-gene mucosal OSM/OSMR axis across anti-TNF cohorts.
    "West_2017_OSM_OSMR": ["OSM", "OSMR"],
}

ALL_SIGNATURES = {**MODULES, **PUBLISHED_SIGNATURES}


def strip_value(value: object, prefix: str) -> str:
    text = "" if pd.isna(value) else str(value)
    return text.split(":", 1)[1].strip() if text.lower().startswith(prefix.lower() + ":") else text.strip()


def eligible_metadata(gse: str, scores: pd.DataFrame) -> pd.DataFrame:
    x = scores.copy()
    if gse == "GSE16879":
        x["tissue"] = x["characteristics_ch1"].map(lambda z: strip_value(z, "tissue"))
        x["disease"] = x["characteristics_ch1_2"].map(lambda z: strip_value(z, "disease"))
        x["response"] = x["characteristics_ch1_3"].map(lambda z: strip_value(z, "response to infliximab"))
        x["visit"] = x["characteristics_ch1_4"].map(
            lambda z: strip_value(z, "before or after first infliximab treatment")
        )
        x = x[(x.tissue == "Colon") & (x.visit == "Before first infliximab treatment")]
        x = x[x.response.isin(["Yes", "No"])]
        x["stratum"] = "GSE16879_" + x["disease"].astype(str) + "_colon"
        x["activity"] = np.nan
        x["age"] = np.nan
    elif gse == "GSE23597":
        x["dose"] = x["characteristics_ch1_2"].map(lambda z: strip_value(z, "dose"))
        x["visit"] = x["characteristics_ch1_3"].map(lambda z: strip_value(z, "time"))
        x["response"] = x["characteristics_ch1_5"].map(lambda z: strip_value(z, "wk8 response"))
        x = x[(x.visit == "W0") & (x.dose != "placebo") & x.response.isin(["Yes", "No"])]
        x["disease"] = "UC"
        x["tissue"] = "Colon"
        x["stratum"] = "GSE23597_UC_colon"
        x["activity"] = np.nan
        x["age"] = np.nan
    elif gse == "GSE92415":
        x["disease"] = x["characteristics_ch1_3"].map(lambda z: strip_value(z, "disease"))
        x["age"] = pd.to_numeric(x["characteristics_ch1_4"].map(lambda z: strip_value(z, "age")), errors="coerce")
        x["treatment"] = x["characteristics_ch1_5"].map(lambda z: strip_value(z, "treatment"))
        x["visit"] = x["characteristics_ch1_6"].map(lambda z: strip_value(z, "visit"))
        x["response"] = x["characteristics_ch1_7"].map(lambda z: strip_value(z, "wk6response"))
        x["activity"] = pd.to_numeric(
            x["characteristics_ch1_8"].map(lambda z: strip_value(z, "mayo score")), errors="coerce"
        )
        x = x[(x.treatment == "golimumab") & (x.visit == "Week 0") & x.response.isin(["Yes", "No"])]
        x["tissue"] = "Colon"
        x["stratum"] = "GSE92415_UC_colon"
    else:
        raise ValueError(gse)
    x["response_binary"] = (x["response"] == "Yes").astype(int)
    return x.reset_index(drop=True)


def load_emtab7604() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load Expression Atlas raw counts and long-form SDRF for E-MTAB-7604."""
    base = ROOT / "data" / "bulk_geo" / "E-MTAB-7604"
    counts = pd.read_csv(base / "E-MTAB-7604-raw-counts.tsv", sep="\t")
    counts = counts.dropna(subset=["Gene Name"]).copy()
    sample_cols = [c for c in counts.columns if c.startswith("ERR")]
    numeric = counts[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    numeric.index = counts["Gene Name"].astype(str).str.upper()
    numeric = numeric.groupby(level=0).sum()
    libsize = numeric.sum(axis=0)
    logcpm = np.log2(numeric.div(libsize, axis=1) * 1_000_000 + 0.5)

    sdrf = pd.read_csv(
        base / "E-MTAB-7604.condensed-sdrf.tsv", sep="\t", header=None,
        names=["study", "blank", "sample", "kind", "field", "value", "ontology"],
    )
    meta = sdrf[sdrf.kind.eq("characteristic")].pivot_table(
        index="sample", columns="field", values="value", aggfunc="first"
    ).reset_index()
    meta["response"] = meta["clinical history"].map(
        {"responder": "Yes", "non-responder": "No"}
    )
    meta["response_binary"] = (meta.response == "Yes").astype(int)
    meta["tissue"] = meta["organism part"].str.capitalize()
    meta["disease"] = "IBD"
    meta["stratum"] = "E-MTAB-7604_IBD_" + meta["organism part"].str.lower()
    meta["activity"] = np.nan
    meta["age"] = np.nan
    meta = meta.rename(columns={"sample": "geo_accession"})
    return logcpm, meta


def gene_z_scores(gene_expr: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    sub = gene_expr.loc[:, sample_ids].astype(float)
    means = sub.mean(axis=1)
    sds = sub.std(axis=1, ddof=1).replace(0, np.nan)
    return sub.sub(means, axis=0).div(sds, axis=0)


def cliffs_delta(response: np.ndarray, nonresponse: np.ndarray) -> float:
    return float((np.greater.outer(response, nonresponse).sum() - np.less.outer(response, nonresponse).sum()) /
                 (len(response) * len(nonresponse)))


def auc_nonresponse_high(response: np.ndarray, nonresponse: np.ndarray) -> float:
    # Probability that a randomly selected non-responder has a higher score than a responder.
    greater = np.greater.outer(nonresponse, response).sum()
    ties = np.equal.outer(nonresponse, response).sum()
    return float((greater + 0.5 * ties) / (len(response) * len(nonresponse)))


def hedges_g(response: np.ndarray, nonresponse: np.ndarray) -> tuple[float, float]:
    n1, n0 = len(response), len(nonresponse)
    df = n1 + n0 - 2
    pooled = math.sqrt(((n1 - 1) * np.var(response, ddof=1) + (n0 - 1) * np.var(nonresponse, ddof=1)) / df)
    d = (np.mean(response) - np.mean(nonresponse)) / pooled
    j = 1 - 3 / (4 * df - 1)
    g = j * d
    var_g = (n1 + n0) / (n1 * n0) + (g * g) / (2 * df)
    return float(g), float(var_g)


def reml_meta(effects: np.ndarray, variances: np.ndarray) -> dict[str, float]:
    k = len(effects)
    w_fixed = 1 / variances
    mu_fixed = np.sum(w_fixed * effects) / np.sum(w_fixed)
    q = float(np.sum(w_fixed * (effects - mu_fixed) ** 2))
    q_p = float(chi2.sf(q, k - 1))
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if q > 0 else 0.0

    def nll(tau2: float) -> float:
        w = 1 / (variances + tau2)
        mu = np.sum(w * effects) / np.sum(w)
        return 0.5 * (
            np.sum(np.log(variances + tau2))
            + np.log(np.sum(w))
            + np.sum(w * (effects - mu) ** 2)
        )

    upper = max(10.0, float(np.var(effects, ddof=1) * 10 + np.max(variances)))
    opt = minimize_scalar(nll, bounds=(0.0, upper), method="bounded")
    tau2 = max(0.0, float(opt.x))
    w = 1 / (variances + tau2)
    mu = float(np.sum(w * effects) / np.sum(w))
    se_naive = math.sqrt(1 / np.sum(w))
    # Hartung-Knapp variance estimator and t interval for small k.
    q_re = float(np.sum(w * (effects - mu) ** 2))
    hk_scale = max(q_re / (k - 1), 1e-12)
    se_hk = math.sqrt(hk_scale / np.sum(w))
    crit = float(t.ppf(0.975, k - 1))
    return {
        "k": k,
        "pooled_g_reml": mu,
        "se_naive": se_naive,
        "ci95_low_hk": mu - crit * se_hk,
        "ci95_high_hk": mu + crit * se_hk,
        "p_hk": float(2 * t.sf(abs(mu / se_hk), k - 1)) if se_hk > 0 else np.nan,
        "tau2_reml": tau2,
        "Q": q,
        "Q_p": q_p,
        "I2_percent": i2,
    }


def logistic_rows(gse92415: pd.DataFrame, signature_cols: list[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for signature in signature_cols:
        for adjusted, covars in [(False, []), (True, ["activity", "age"])]:
            cols = [signature, "response_binary", *covars]
            d = gse92415[cols].dropna().copy()
            for col in [signature, *covars]:
                sd = d[col].std(ddof=1)
                d[col] = (d[col] - d[col].mean()) / sd
            X = sm.add_constant(d[[signature, *covars]], has_constant="add")
            try:
                fit = sm.GLM(d.response_binary, X, family=sm.families.Binomial()).fit(cov_type="HC3")
                beta = float(fit.params[signature]); se = float(fit.bse[signature])
                rows.append({
                    "signature": signature,
                    "model": "adjusted_Mayo_age" if adjusted else "unadjusted",
                    "n": len(d),
                    "n_responder": int(d.response_binary.sum()),
                    "n_non_responder": int((1 - d.response_binary).sum()),
                    "OR_response_per_SD": math.exp(beta),
                    "CI95_low": math.exp(beta - 1.96 * se),
                    "CI95_high": math.exp(beta + 1.96 * se),
                    "p_value": float(fit.pvalues[signature]),
                })
            except Exception as exc:
                rows.append({"signature": signature, "model": "adjusted_Mayo_age" if adjusted else "unadjusted",
                             "n": len(d), "error": str(exc)})
    return rows


def emtab_site_adjusted_rows(d0: pd.DataFrame, signature_cols: list[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for signature in signature_cols:
        d = d0[[signature, "response_binary", "tissue", "treatment"]].dropna().copy()
        d[signature] = (d[signature] - d[signature].mean()) / d[signature].std(ddof=1)
        d["ileum"] = (d.tissue == "Ileum").astype(int)
        d["infliximab"] = (d.treatment == "infliximab").astype(int)
        X = sm.add_constant(d[[signature, "ileum", "infliximab"]], has_constant="add")
        try:
            fit = sm.GLM(d.response_binary, X, family=sm.families.Binomial()).fit(cov_type="HC3")
            beta = float(fit.params[signature]); se = float(fit.bse[signature])
            rows.append({
                "signature": signature, "model": "adjusted_tissue_treatment", "n": len(d),
                "n_responder": int(d.response_binary.sum()),
                "n_non_responder": int((1 - d.response_binary).sum()),
                "OR_response_per_SD": math.exp(beta),
                "CI95_low": math.exp(beta - 1.96 * se),
                "CI95_high": math.exp(beta + 1.96 * se),
                "p_value": float(fit.pvalues[signature]),
            })
        except Exception as exc:
            rows.append({"signature": signature, "model": "adjusted_tissue_treatment",
                         "n": len(d), "error": str(exc)})
    return rows


def paired_auc_bootstrap(d: pd.DataFrame, n_boot: int = 10000) -> pd.DataFrame:
    orientation = {
        "myeloid_inflammation": 1,
        "epithelial_ifn_damage": 1,
        "epithelial_barrier_maturity": -1,
        "stromal_fibroinflammatory": 1,
        **{name: 1 for name in PUBLISHED_SIGNATURES},
    }
    rng = np.random.default_rng(20260805)
    responders = np.flatnonzero(d.response_binary.to_numpy() == 1)
    nonresponders = np.flatnonzero(d.response_binary.to_numpy() == 0)
    rows = []
    for current in MODULES:
        current_values = orientation[current] * d[current].to_numpy()
        current_auc = auc_nonresponse_high(current_values[responders], current_values[nonresponders])
        for benchmark in PUBLISHED_SIGNATURES:
            benchmark_values = orientation[benchmark] * d[benchmark].to_numpy()
            benchmark_auc = auc_nonresponse_high(
                benchmark_values[responders], benchmark_values[nonresponders]
            )
            diffs = np.empty(n_boot)
            for i in range(n_boot):
                ri = rng.choice(responders, size=len(responders), replace=True)
                ni = rng.choice(nonresponders, size=len(nonresponders), replace=True)
                a1 = auc_nonresponse_high(current_values[ri], current_values[ni])
                a0 = auc_nonresponse_high(benchmark_values[ri], benchmark_values[ni])
                diffs[i] = a1 - a0
            rows.append({
                "current_module": current,
                "published_signature": benchmark,
                "n_responder": len(responders),
                "n_non_responder": len(nonresponders),
                "oriented_AUC_current": current_auc,
                "oriented_AUC_published": benchmark_auc,
                "delta_AUC_current_minus_published": current_auc - benchmark_auc,
                "bootstrap_CI95_low": float(np.quantile(diffs, 0.025)),
                "bootstrap_CI95_high": float(np.quantile(diffs, 0.975)),
                "bootstrap_two_sided_p": float(2 * min(np.mean(diffs <= 0), np.mean(diffs >= 0))),
                "n_bootstrap": n_boot,
            })
    return pd.DataFrame(rows)


def main() -> None:
    cohort_specs = {
        "GSE16879": (ROOT / "data/bulk_geo/GSE16879_series_matrix.txt.gz", "GPL570"),
        "GSE23597": (ROOT / "data/bulk_geo/GSE23597_series_matrix.txt.gz", "GPL570"),
        "GSE92415": (ROOT / "data/bulk_geo/GSE92415_series_matrix.txt.gz", "GPL13158"),
    }
    effect_rows: list[dict[str, object]] = []
    availability_rows: list[dict[str, object]] = []
    eligible_frames: dict[str, pd.DataFrame] = {}
    measured_universe: set[str] = set()

    for gse, (path, platform) in cohort_specs.items():
        expr, meta, _ = parse_series_matrix(path)
        mapping = load_platform_mapping(platform, download=False)
        gene_expr = collapse_to_gene(expr, mapping)
        measured_universe.update(gene_expr.index)
        blank_scores = pd.DataFrame({"geo_accession": gene_expr.columns}).merge(meta, on="geo_accession", how="left")
        eligible = eligible_metadata(gse, blank_scores)
        z = gene_z_scores(gene_expr, eligible.geo_accession.tolist())
        for name, genes in ALL_SIGNATURES.items():
            present = [gene for gene in genes if gene in z.index]
            eligible[name] = z.loc[present, eligible.geo_accession].mean(axis=0).to_numpy() if present else np.nan
            availability_rows.append({
                "cohort": gse,
                "signature": name,
                "n_requested": len(genes),
                "n_present": len(present),
                "coverage": len(present) / len(genes),
                "present_genes": ";".join(present),
                "missing_genes": ";".join(g for g in genes if g not in present),
            })
        eligible_frames[gse] = eligible

        for stratum, d in eligible.groupby("stratum"):
            for name in ALL_SIGNATURES:
                r = d.loc[d.response_binary == 1, name].dropna().to_numpy()
                nr = d.loc[d.response_binary == 0, name].dropna().to_numpy()
                g, var_g = hedges_g(r, nr)
                effect_rows.append({
                    "cohort": gse,
                    "stratum": stratum,
                    "signature": name,
                    "signature_type": "current_module" if name in MODULES else "published_benchmark",
                    "n_responder": len(r),
                    "n_non_responder": len(nr),
                    "mean_responder": float(np.mean(r)),
                    "mean_non_responder": float(np.mean(nr)),
                    "cliffs_delta_R_vs_NR": cliffs_delta(r, nr),
                    "AUC_nonresponse_high": auc_nonresponse_high(r, nr),
                    "hedges_g_R_minus_NR": g,
                    "variance_g": var_g,
                    "se_g": math.sqrt(var_g),
                    "ci95_low_g": g - 1.96 * math.sqrt(var_g),
                    "ci95_high_g": g + 1.96 * math.sqrt(var_g),
                })

    # Independent RNA-seq cohort from Expression Atlas (19 R, 25 NR), with
    # explicit ileum/colon and anti-TNF agent metadata.
    gene_expr, eligible = load_emtab7604()
    measured_universe.update(gene_expr.index)
    z = gene_z_scores(gene_expr, eligible.geo_accession.tolist())
    for name, genes in ALL_SIGNATURES.items():
        present = [gene for gene in genes if gene in z.index]
        eligible[name] = z.loc[present, eligible.geo_accession].mean(axis=0).to_numpy() if present else np.nan
        availability_rows.append({
            "cohort": "E-MTAB-7604", "signature": name, "n_requested": len(genes),
            "n_present": len(present), "coverage": len(present) / len(genes),
            "present_genes": ";".join(present),
            "missing_genes": ";".join(g for g in genes if g not in present),
        })
    eligible_frames["E-MTAB-7604"] = eligible
    for stratum, d in eligible.groupby("stratum"):
        for name in ALL_SIGNATURES:
            r = d.loc[d.response_binary == 1, name].dropna().to_numpy()
            nr = d.loc[d.response_binary == 0, name].dropna().to_numpy()
            g, var_g = hedges_g(r, nr)
            effect_rows.append({
                "cohort": "E-MTAB-7604", "stratum": stratum, "signature": name,
                "signature_type": "current_module" if name in MODULES else "published_benchmark",
                "n_responder": len(r), "n_non_responder": len(nr),
                "mean_responder": float(np.mean(r)), "mean_non_responder": float(np.mean(nr)),
                "cliffs_delta_R_vs_NR": cliffs_delta(r, nr),
                "AUC_nonresponse_high": auc_nonresponse_high(r, nr),
                "hedges_g_R_minus_NR": g, "variance_g": var_g, "se_g": math.sqrt(var_g),
                "ci95_low_g": g - 1.96 * math.sqrt(var_g),
                "ci95_high_g": g + 1.96 * math.sqrt(var_g),
            })

    effects = pd.DataFrame(effect_rows)
    effects.to_csv(OUT / "cross_cohort_signature_effects.csv", index=False)
    pd.DataFrame(availability_rows).to_csv(OUT / "signature_gene_availability.csv", index=False)

    meta_rows = []
    for signature, d in effects.groupby("signature"):
        # Secondary analysis: all disease/tissue strata.
        result = reml_meta(d.hedges_g_R_minus_NR.to_numpy(), d.variance_g.to_numpy())
        meta_rows.append({"signature": signature, "signature_type": d.signature_type.iloc[0],
                          "meta_scope": "all_strata_secondary", **result})
        # Primary analysis: collapse the two disjoint GSE16879 disease strata into one
        # study-level fixed-effect estimate, yielding four independent cohorts.
        study_effects = []
        study_vars = []
        for cohort, c in d.groupby("cohort"):
            w = 1 / c.variance_g.to_numpy()
            study_effects.append(float(np.sum(w * c.hedges_g_R_minus_NR.to_numpy()) / np.sum(w)))
            study_vars.append(float(1 / np.sum(w)))
        result_primary = reml_meta(np.asarray(study_effects), np.asarray(study_vars))
        meta_rows.append({"signature": signature, "signature_type": d.signature_type.iloc[0],
                          "meta_scope": "four_independent_cohorts_primary", **result_primary})

        # Colon-only sensitivity analysis; prevents the E-MTAB ileal subgroup
        # from changing the anatomical target of the main mucosal comparison.
        colon = d[d.stratum.str.endswith("_colon")]
        study_effects = []
        study_vars = []
        for cohort, c in colon.groupby("cohort"):
            w = 1 / c.variance_g.to_numpy()
            study_effects.append(float(np.sum(w * c.hedges_g_R_minus_NR.to_numpy()) / np.sum(w)))
            study_vars.append(float(1 / np.sum(w)))
        if len(study_effects) >= 3:
            result_colon = reml_meta(np.asarray(study_effects), np.asarray(study_vars))
            meta_rows.append({"signature": signature, "signature_type": d.signature_type.iloc[0],
                              "meta_scope": "four_independent_colon_cohorts_sensitivity",
                              **result_colon})
    pd.DataFrame(meta_rows).to_csv(OUT / "random_effects_meta_analysis.csv", index=False)

    overlap_rows = []
    universe_n = max(len(measured_universe), 1)
    for current_name, current_genes in MODULES.items():
        a = set(current_genes)
        for published_name, published_genes in PUBLISHED_SIGNATURES.items():
            b = set(published_genes)
            overlap = a & b
            union = a | b
            p_hyper = float(hypergeom.sf(len(overlap) - 1, universe_n, len(a), len(b)))
            overlap_rows.append({
                "current_module": current_name,
                "published_signature": published_name,
                "n_current": len(a), "n_published": len(b), "n_overlap": len(overlap),
                "overlap_genes": ";".join(sorted(overlap)),
                "jaccard": len(overlap) / len(union),
                "overlap_coefficient": len(overlap) / min(len(a), len(b)),
                "hypergeometric_p": p_hyper,
                "gene_universe_n": universe_n,
            })
    pd.DataFrame(overlap_rows).to_csv(OUT / "module_vs_published_signature_overlap.csv", index=False)

    adjusted = pd.DataFrame(logistic_rows(eligible_frames["GSE92415"], list(ALL_SIGNATURES)))
    adjusted.to_csv(OUT / "GSE92415_disease_activity_adjusted_models.csv", index=False)
    pd.DataFrame(emtab_site_adjusted_rows(eligible_frames["E-MTAB-7604"], list(ALL_SIGNATURES))).to_csv(
        OUT / "E-MTAB-7604_tissue_treatment_adjusted_models.csv", index=False
    )
    paired_auc_bootstrap(eligible_frames["GSE92415"]).to_csv(
        OUT / "GSE92415_paired_auc_comparisons.csv", index=False
    )
    paired_auc_bootstrap(eligible_frames["E-MTAB-7604"]).to_csv(
        OUT / "E-MTAB-7604_paired_auc_comparisons.csv", index=False
    )

    cohort_summary = []
    for gse, d in eligible_frames.items():
        for stratum, s in d.groupby("stratum"):
            cohort_summary.append({
                "cohort": gse, "stratum": stratum, "n": len(s),
                "n_responder": int(s.response_binary.sum()),
                "n_non_responder": int((1 - s.response_binary).sum()),
                "tissue": ";".join(sorted(s.tissue.dropna().astype(str).unique())),
                "disease": ";".join(sorted(s.disease.dropna().astype(str).unique())),
                "baseline_activity_available": bool(s.activity.notna().all()),
                "activity_variable": "Total Mayo score" if gse == "GSE92415" else "not in public sample metadata",
            })
    pd.DataFrame(cohort_summary).to_csv(OUT / "included_cohort_summary.csv", index=False)

    manifest = {
        "analysis_date": "2026-08-05",
        "effect_measure": "Hedges g, responder minus non-responder, from within-stratum gene-z-score module means",
        "meta_method": "REML random effects; Hartung-Knapp 95% CI; Q and I2 heterogeneity",
        "activity_adjustment": "GSE92415 logistic regression adjusted for baseline total Mayo score and age",
        "site_adjustment": "E-MTAB-7604 logistic regression adjusted for ileum/colon and anti-TNF agent; colon-only meta-analysis sensitivity",
        "published_signatures": PUBLISHED_SIGNATURES,
        "outputs": sorted(p.name for p in OUT.glob("*.csv")),
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
