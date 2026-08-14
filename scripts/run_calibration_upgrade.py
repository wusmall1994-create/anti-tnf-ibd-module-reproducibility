from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests


ROOT = Path(
    os.environ.get("ANTI_TNF_IBD_ROOT", Path(__file__).resolve().parents[1])
).resolve()
OUT = ROOT / "results" / "calibration_upgrade_20260814"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "scripts"))

from run_enhancement_gate import (  # noqa: E402
    eligible_metadata,
    gene_z_scores,
    hedges_g,
    reml_meta,
    load_emtab7604,
)
from score_bulk_series_matrix_modules import (  # noqa: E402
    MODULES,
    collapse_to_gene,
    load_platform_mapping,
    parse_series_matrix,
)


# Exact gene lists transcribed from Gaujoux et al. 2019 Supplementary Table S1,
# plus two independently published pathway signatures with explicit gene lists.
BENCHMARKS = {
    "Gaujoux_UC_A_20": ["IL13RA2", "CD86", "NR3C1", "TMEM158", "TNFRSF11B", "RGS5", "CD69", "SPG20", "IGFBP5", "WNT5A", "DSE", "ACSL4", "STC1", "ZEB2", "CXCL6", "RBMS1", "EDNRB", "MMP10", "ADGRL2", "ADGRE2"],
    "Gaujoux_UC_B_20": ["TNFRSF11B", "STC1", "CSF3R", "IL1B", "MASP1", "PTGS2", "G0S2", "OSM", "AQP9", "TLR2", "APOBEC3A", "FPR1", "KCNJ15", "FCGR1B", "BCL6", "C5AR1", "GK", "HCAR3", "CXCR2", "CEMIP"],
    "Gaujoux_UC_B_knn_19": ["TREM1", "OSM", "G0S2", "PDE4B", "BCL6", "CREB5", "IL1B", "CSF3R", "TLR2", "C5AR1", "FCN1", "LILRA5", "KCNJ15", "FGR", "LILRA1", "GK3P", "HCAR3", "FPR2", "NAMPT"],
    "Gaujoux_IRRAT_29": ["ADAM9", "ADAMTS1", "AKAP12", "CDH6", "CTSS", "EGR1", "EVI2A", "FCGR3A", "FOS", "ITGB3", "ITGB6", "LCN2", "LTF", "MEGF11", "NFKBIZ", "NNMT", "OLFM4", "OSMR", "PI15", "PTPRC", "PTX3", "RARRES1", "S100A8", "SERPINA3", "SLPI", "SOD2", "VCAN", "TMEM252", "MT-ND6"],
    "Gaujoux_UC_AB_53": ["ACSL4", "BCAT1", "CCL4", "CLEC4A", "COL12A1", "CSF2RB", "CXCL5", "CXCL6", "FAM129A", "FCER1G", "FCGR2A", "FGF7", "GLIS3", "IGFBP5", "IGSF6", "IL11", "IL13RA2", "IL6", "IL7R", "INHBA", "KLHL5", "LCP2", "LILRB1", "LILRB2", "LY96", "MCTP1", "MME", "MMP3", "PAPPA", "PDE4B", "PI15", "PRR16", "PTGS2", "RGS18", "RGS5", "S100A9", "SAMSN1", "SELE", "SERPINE1", "SLC2A3", "SNX10", "SOCS3", "SRGN", "STC1", "TAGAP", "TFPI", "TFPI2", "TLR1", "TNC", "TNFRSF11B", "WNT5A", "ADGRE2", "ADGRL2"],
    "Gaujoux_CDc_20": ["MNDA", "S100A12", "IL13RA2", "IL11", "FCGR2B", "PTGES", "NCF2", "LILRB2", "IL6", "TNFAIP6", "TAGAP", "SLC2A3", "FCN1", "S100A8", "PROK2", "FPR1", "BCL2A1", "G0S2", "CD14", "S100A9"],
    "West_OSM22": ["OSM", "IL1B", "IL1A", "IL6", "IL11", "CSF2", "CSF3", "IFNG", "IL17A", "IL22", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CCL2", "CCL3", "CCL4", "CXCL9", "CXCL10", "CXCL11"],
    "Belarif_IL7R_10": ["IL7R", "IL2RG", "JAK1", "PIK3CA", "LCK", "PTK2B", "EP300", "NMI", "CRLF2", "TSLP"],
}

DEVELOPMENT_OVERLAP = {
    "Gaujoux_UC_A_20": {"GSE16879"},  # GSE14580 is the UC subset of GSE16879.
    "Gaujoux_UC_B_20": {"GSE23597"},  # GSE12251 is an exact subset of GSE23597.
    "Gaujoux_UC_B_knn_19": {"GSE23597"},
    "Gaujoux_IRRAT_29": {"GSE23597"},
    "Gaujoux_UC_AB_53": {"GSE16879", "GSE23597"},
    "Gaujoux_CDc_20": {"GSE16879"},
    "West_OSM22": set(),
    "Belarif_IL7R_10": {"GSE16879", "GSE23597"},
}

ALL = {**MODULES, **BENCHMARKS}
ORIENTATION = {name: 1 for name in ALL}  # higher score expected in non-response
ORIENTATION["epithelial_barrier_maturity"] = -1


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def add_scores(gene_expr: pd.DataFrame, meta: pd.DataFrame, cohort: str):
    z = gene_z_scores(gene_expr, meta.geo_accession.tolist())
    availability = []
    for name, genes in ALL.items():
        present = [g for g in genes if g in z.index]
        meta[name] = z.loc[present, meta.geo_accession].mean(axis=0).to_numpy() if present else np.nan
        availability.append({
            "cohort": cohort, "signature": name, "n_requested": len(genes),
            "n_present": len(present), "coverage": len(present) / len(genes),
            "present_genes": ";".join(present),
            "missing_genes": ";".join(g for g in genes if g not in present),
        })
    return meta, availability


def load_antitnf_cohorts():
    specs = {
        "GSE16879": (ROOT / "data/bulk_geo/GSE16879_series_matrix.txt.gz", "GPL570"),
        "GSE23597": (ROOT / "data/bulk_geo/GSE23597_series_matrix.txt.gz", "GPL570"),
        "GSE92415": (ROOT / "data/bulk_geo/GSE92415_series_matrix.txt.gz", "GPL13158"),
    }
    frames, genes, availability = {}, {}, []
    for cohort, (path, platform) in specs.items():
        expr, raw_meta, _ = parse_series_matrix(path)
        gene_expr = collapse_to_gene(expr, load_platform_mapping(platform, download=False))
        blank = pd.DataFrame({"geo_accession": gene_expr.columns}).merge(raw_meta, on="geo_accession", how="left")
        meta = eligible_metadata(cohort, blank)
        if cohort == "GSE16879": meta["drug"] = "infliximab"
        elif cohort == "GSE23597": meta["drug"] = "infliximab"
        else: meta["drug"] = "golimumab"
        meta, a = add_scores(gene_expr, meta, cohort)
        frames[cohort], genes[cohort] = meta, gene_expr
        availability.extend(a)

    gene_expr, meta = load_emtab7604()
    meta["drug"] = meta["treatment"].str.lower()
    meta, a = add_scores(gene_expr, meta, "E-MTAB-7604")
    frames["E-MTAB-7604"], genes["E-MTAB-7604"] = meta, gene_expr
    availability.extend(a)
    return frames, genes, pd.DataFrame(availability)


def effect_row(cohort, stratum, signature, d):
    r = d.loc[d.response_binary == 1, signature].dropna().to_numpy()
    nr = d.loc[d.response_binary == 0, signature].dropna().to_numpy()
    if len(r) < 2 or len(nr) < 2:
        return None
    g, var = hedges_g(r, nr)
    oriented = ORIENTATION[signature] * d[signature].to_numpy()
    auc = roc_auc_score(d.response_binary, -oriented)
    p = mannwhitneyu(r, nr, alternative="two-sided").pvalue
    return {
        "cohort": cohort, "stratum": stratum, "signature": signature,
        "signature_type": "current_module" if signature in MODULES else "published_benchmark",
        "development_overlap": cohort in DEVELOPMENT_OVERLAP.get(signature, set()),
        "n_responder": len(r), "n_non_responder": len(nr),
        "hedges_g_R_minus_NR": g, "variance_g": var,
        "ci95_low_g": g - 1.96 * math.sqrt(var), "ci95_high_g": g + 1.96 * math.sqrt(var),
        "oriented_AUC_response": auc, "mannwhitney_p": p,
    }


def benchmark(frames):
    rows = []
    for cohort, d in frames.items():
        for stratum, s in d.groupby("stratum"):
            for signature in ALL:
                row = effect_row(cohort, stratum, signature, s)
                if row: rows.append(row)
    effects = pd.DataFrame(rows)
    effects.to_csv(OUT / "benchmark_stratum_effects.csv", index=False)

    study_rows = []
    for (signature, cohort), d in effects.groupby(["signature", "cohort"]):
        w = 1 / d.variance_g.to_numpy()
        g = float(np.sum(w * d.hedges_g_R_minus_NR) / np.sum(w))
        var = float(1 / np.sum(w))
        study_rows.append({"signature": signature, "cohort": cohort, "g": g, "variance": var,
                           "development_overlap": bool(d.development_overlap.any())})
    studies = pd.DataFrame(study_rows)
    studies.to_csv(OUT / "benchmark_study_level_effects.csv", index=False)

    meta_rows, loco_rows = [], []
    for signature, d in studies.groupby("signature"):
        for scope, x in [("all_four_cohorts", d), ("external_only", d[~d.development_overlap])]:
            if len(x) >= 2:
                m = reml_meta(x.g.to_numpy(), x.variance.to_numpy())
                meta_rows.append({"signature": signature, "scope": scope, **m})
        if len(d) >= 3:
            for held in d.cohort:
                x = d[d.cohort != held]
                m = reml_meta(x.g.to_numpy(), x.variance.to_numpy())
                loco_rows.append({"signature": signature, "held_out_cohort": held, **m})
    meta = pd.DataFrame(meta_rows)
    meta.to_csv(OUT / "benchmark_meta_analysis.csv", index=False)
    pd.DataFrame(loco_rows).to_csv(OUT / "benchmark_leave_one_cohort_out.csv", index=False)

    ratings = []
    for signature, d in studies.groupby("signature"):
        ext = d[~d.development_overlap]
        dirs = (ORIENTATION[signature] * ext.g < 0).mean() if len(ext) else np.nan
        mm = meta[(meta.signature == signature) & (meta.scope == "external_only")]
        pooled = mm.iloc[0] if len(mm) else None
        decisive_ci = False
        if pooled is not None:
            if ORIENTATION[signature] == 1:
                decisive_ci = pooled.ci95_high_hk < 0
            else:
                decisive_ci = pooled.ci95_low_hk > 0
        if pooled is not None and len(ext) >= 3 and dirs == 1 and decisive_ci:
            grade = "robust_cross_cohort"
        elif pooled is not None and len(ext) >= 3 and dirs == 1 and pooled.pooled_g_reml * ORIENTATION[signature] < 0:
            grade = "directionally_consistent"
        elif pooled is not None and len(ext) >= 2 and dirs >= 0.67 and pooled.pooled_g_reml * ORIENTATION[signature] < 0:
            grade = "limited_external_evidence"
        else:
            grade = "insufficient"
        ratings.append({"signature": signature, "n_external_cohorts": len(ext),
                        "external_direction_consistency": dirs, "evidence_grade": grade})
    ratings = pd.DataFrame(ratings)
    ratings.to_csv(OUT / "benchmark_evidence_grades.csv", index=False)
    return effects, studies, meta, ratings


def subgroup_analysis(frames):
    rows = []
    for cohort, d in frames.items():
        d = d.copy()
        if cohort == "GSE16879": d["disease_group"] = d["disease"].replace({"Ulcerative Colitis (UC)": "UC"})
        elif cohort in {"GSE23597", "GSE92415"}: d["disease_group"] = "UC"
        else: d["disease_group"] = "mixed_IBD"
        d["tissue_group"] = d["tissue"].str.lower().replace({"colon": "colon", "ileum": "ileum"})
        for axis in ["disease_group", "tissue_group", "drug"]:
            for level, s in d.groupby(axis):
                for signature in MODULES:
                    row = effect_row(cohort, f"{axis}:{level}", signature, s)
                    if row:
                        row.update({"axis": axis, "level": level})
                        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "module_disease_tissue_drug_strata.csv", index=False)
    return out


def clinical_increment(gse92415: pd.DataFrame, repeats=100):
    models = {
        "clinical_Mayo_age": ["activity", "age"],
        "clinical_plus_myeloid": ["activity", "age", "myeloid_inflammation"],
        "clinical_plus_barrier": ["activity", "age", "epithelial_barrier_maturity"],
        "clinical_plus_myeloid_barrier": ["activity", "age", "myeloid_inflammation", "epithelial_barrier_maturity"],
    }
    cols = sorted(set(sum(models.values(), [])) | {"response_binary"})
    d = gse92415[cols].dropna().reset_index(drop=True)
    y = d.response_binary.to_numpy()
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=repeats, random_state=20260814)
    pred_sum = {m: np.zeros(len(d)) for m in models}
    pred_n = {m: np.zeros(len(d)) for m in models}
    fold_rows = []
    for split, (tr, te) in enumerate(cv.split(d, y)):
        for name, features in models.items():
            pipe = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, solver="liblinear", max_iter=2000))
            pipe.fit(d.loc[tr, features], y[tr])
            p = pipe.predict_proba(d.loc[te, features])[:, 1]
            pred_sum[name][te] += p; pred_n[name][te] += 1
            fold_rows.append({"split": split, "model": name, "auc": roc_auc_score(y[te], p),
                              "brier": brier_score_loss(y[te], p)})
    fold = pd.DataFrame(fold_rows)
    fold.to_csv(OUT / "GSE92415_incremental_model_cv_folds.csv", index=False)
    pred = pd.DataFrame({"response_binary": y})
    summary = []
    for name in models:
        p = pred_sum[name] / pred_n[name]
        pred[name] = p
        clipped = np.clip(p, 1e-6, 1 - 1e-6)
        lp = np.log(clipped / (1 - clipped)).reshape(-1, 1)
        calibration_model = LogisticRegression(C=np.inf, solver="lbfgs", max_iter=2000)
        calibration_model.fit(lp, y)
        slope = float(calibration_model.coef_[0, 0])
        summary.append({"model": name, "n": len(y), "repeats": repeats,
                        "oof_auc": roc_auc_score(y, p), "oof_brier": brier_score_loss(y, p),
                        "logistic_calibration_slope": slope,
                        "fold_auc_q025": fold.loc[fold.model == name, "auc"].quantile(.025),
                        "fold_auc_q975": fold.loc[fold.model == name, "auc"].quantile(.975)})
    pred.to_csv(OUT / "GSE92415_incremental_model_oof_predictions.csv", index=False)
    summary = pd.DataFrame(summary)
    base_auc = summary.loc[summary.model == "clinical_Mayo_age", "oof_auc"].iloc[0]
    base_brier = summary.loc[summary.model == "clinical_Mayo_age", "oof_brier"].iloc[0]
    summary["delta_oof_auc_vs_clinical"] = summary.oof_auc - base_auc
    summary["delta_oof_brier_vs_clinical"] = summary.oof_brier - base_brier
    summary.to_csv(OUT / "GSE92415_incremental_model_summary.csv", index=False)

    rng = np.random.default_rng(20260814)
    base = pred["clinical_Mayo_age"].to_numpy()
    reclassification = []
    for name in models:
        if name == "clinical_Mayo_age":
            continue
        new = pred[name].to_numpy()
        def metrics(idx):
            yy, bb, nn = y[idx], base[idx], new[idx]
            events, nonevents = yy == 1, yy == 0
            nri = (np.mean(nn[events] > bb[events]) - np.mean(nn[events] < bb[events])
                   + np.mean(nn[nonevents] < bb[nonevents]) - np.mean(nn[nonevents] > bb[nonevents]))
            ds_base = bb[events].mean() - bb[nonevents].mean()
            ds_new = nn[events].mean() - nn[nonevents].mean()
            return nri, ds_new - ds_base
        point = metrics(np.arange(len(y)))
        boots = []
        for _ in range(2000):
            e = rng.choice(np.where(y == 1)[0], np.sum(y == 1), replace=True)
            ne = rng.choice(np.where(y == 0)[0], np.sum(y == 0), replace=True)
            boots.append(metrics(np.r_[e, ne]))
        boots = np.asarray(boots)
        reclassification.append({
            "model": name, "continuous_NRI": point[0], "NRI_ci95_low": np.quantile(boots[:, 0], .025),
            "NRI_ci95_high": np.quantile(boots[:, 0], .975), "IDI": point[1],
            "IDI_ci95_low": np.quantile(boots[:, 1], .025), "IDI_ci95_high": np.quantile(boots[:, 1], .975),
            "bootstrap_resamples": 2000,
        })
    pd.DataFrame(reclassification).to_csv(OUT / "GSE92415_incremental_reclassification.csv", index=False)

    thresholds = np.arange(.10, .81, .05)
    dca = []
    n = len(y)
    for name in models:
        p = pred[name].to_numpy()
        for t in thresholds:
            positive = p >= t
            tp = np.sum(positive & (y == 1)); fp = np.sum(positive & (y == 0))
            nb = tp / n - fp / n * t / (1 - t)
            dca.append({"model": name, "threshold": t, "net_benefit": nb})
    pd.DataFrame(dca).to_csv(OUT / "GSE92415_decision_curve.csv", index=False)


def strip_prefix(s, prefix):
    s = "" if pd.isna(s) else str(s)
    return s.split(":", 1)[1].strip() if s.lower().startswith(prefix.lower() + ":") else s.strip()


def load_vdz_boundary():
    path = ROOT / "data/bulk_geo/GSE73661_series_matrix.txt.gz"
    expr, meta, _ = parse_series_matrix(path)
    mapping = load_platform_mapping("GPL6244", download=True)
    gene_expr = collapse_to_gene(expr, mapping)
    meta["patient"] = meta.characteristics_ch1.map(lambda x: strip_prefix(x, "study individual number"))
    meta["week"] = meta.characteristics_ch1_2.map(lambda x: strip_prefix(x, "week (w)"))
    meta["regimen"] = meta.characteristics_ch1_3.map(lambda x: strip_prefix(x, "induction therapy_maintenance therapy")).str.lower()
    meta["mayo_endoscopy"] = pd.to_numeric(meta.characteristics_ch1_4.map(lambda x: strip_prefix(x, "mayo endoscopic subscore")), errors="coerce")
    meta["is_vdz"] = meta.regimen.str.contains("vdz", na=False)
    base = meta[meta.is_vdz & meta.week.eq("W0")].copy()
    follow = meta[meta.is_vdz & meta.week.isin(["W6", "W12"])].copy()
    follow["follow_order"] = follow.week.map({"W6": 6, "W12": 12})
    follow = follow.sort_values(["patient", "follow_order"]).drop_duplicates("patient")
    outcome = follow[["patient", "mayo_endoscopy", "week"]].rename(
        columns={"mayo_endoscopy": "outcome_mayo_endoscopy", "week": "outcome_week"}
    )
    base = base.merge(outcome, on="patient", how="inner")
    base["response_binary"] = (base.outcome_mayo_endoscopy <= 1).astype(int)
    base["response"] = base.response_binary.map({1: "Yes", 0: "No"})
    base["tissue"] = "Colon"; base["disease"] = "UC"; base["drug"] = "vedolizumab"
    base["stratum"] = "GSE73661_VDZ_UC_colon"
    base, avail = add_scores(gene_expr, base, "GSE73661_VDZ")
    return base, gene_expr, avail


def load_ust_boundary():
    path = ROOT / "data/bulk_geo/GSE206285_series_matrix.txt.gz"
    if not path.exists(): return None
    try:
        expr, meta, _ = parse_series_matrix(path)
    except (EOFError, OSError, ValueError):
        return None
    mapping = load_platform_mapping("GPL13158", download=False)
    gene_expr = collapse_to_gene(expr, mapping)
    # GEO has stable field order, but identify characteristics by their explicit prefixes.
    char_cols = [c for c in meta if c.startswith("characteristics_ch1")]
    def find(prefix):
        for c in char_cols:
            vals = meta[c].astype(str)
            if vals.str.lower().str.startswith(prefix.lower() + ":").mean() > .5:
                return vals.map(lambda x: strip_prefix(x, prefix))
        return pd.Series([np.nan] * len(meta), index=meta.index)
    meta["treatment_clean"] = find("treatment")
    meta["visit_clean"] = find("visit")
    meta["healing"] = find("mucosal healing at week 8")
    d = meta[meta.treatment_clean.str.contains("ustekinumab", case=False, na=False) & meta.visit_clean.str.contains("week.*0|week i-0", case=False, regex=True, na=False)].copy()
    d = d[d.healing.isin(["Y", "N"])]
    d["response_binary"] = (d.healing == "Y").astype(int)
    d["response"] = d.response_binary.map({1: "Yes", 0: "No"})
    d["tissue"] = "Colon"; d["disease"] = "UC"; d["drug"] = "ustekinumab"
    d["stratum"] = "GSE206285_UST_UC_colon"
    d, avail = add_scores(gene_expr, d, "GSE206285_UST")
    return d, gene_expr, avail


def boundary_analysis(vdz, ust):
    rows, availability = [], []
    for cohort, payload in [("GSE73661_VDZ", vdz), ("GSE206285_UST", ust)]:
        if payload is None: continue
        d, _, a = payload; availability.extend(a)
        for signature in ALL:
            row = effect_row(cohort, d.stratum.iloc[0], signature, d)
            if row: rows.append(row)
    effects = pd.DataFrame(rows)
    if len(effects):
        effects["mannwhitney_fdr_within_cohort"] = effects.groupby("cohort")["mannwhitney_p"].transform(
            lambda p: pd.Series(multipletests(p, method="fdr_bh")[1], index=p.index)
        )
    effects.to_csv(OUT / "cross_treatment_boundary_effects.csv", index=False)
    pd.DataFrame(availability).to_csv(OUT / "cross_treatment_signature_availability.csv", index=False)


def expression_overlap_audit():
    pairs = [("GSE12251", "GSE23597"), ("GSE14580", "GSE16879")]
    rows = []
    mapping = load_platform_mapping("GPL570", download=False)
    for small, large in pairs:
        e1, m1, _ = parse_series_matrix(ROOT / f"data/bulk_geo/{small}_series_matrix.txt.gz")
        e2, m2, _ = parse_series_matrix(ROOT / f"data/bulk_geo/{large}_series_matrix.txt.gz")
        g1, g2 = collapse_to_gene(e1, mapping), collapse_to_gene(e2, mapping)
        genes = g1.index.intersection(g2.index)
        # Stable high-variance genes in the larger accession provide expression fingerprints.
        genes = g2.loc[genes].var(axis=1).nlargest(min(3000, len(genes))).index
        a, b = g1.loc[genes], g2.loc[genes]
        az = (a - a.mean(axis=0)) / a.std(axis=0)
        bz = (b - b.mean(axis=0)) / b.std(axis=0)
        corr = az.T.to_numpy() @ bz.to_numpy() / (len(genes) - 1)
        for i, sample in enumerate(a.columns):
            j = int(np.nanargmax(corr[i]))
            rows.append({"smaller_accession": small, "smaller_sample": sample,
                         "larger_accession": large, "best_match_sample": b.columns[j],
                         "pearson_r": float(corr[i, j]), "n_fingerprint_genes": len(genes)})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "dataset_expression_fingerprint_matches.csv", index=False)
    return out


def independence_audit(fingerprint):
    records = []
    specs = [
        ("GSE16879", "anti-TNF", "pretreatment mucosa", "include", "independent primary cohort"),
        ("GSE23597", "anti-TNF", "pretreatment mucosa", "include", "independent primary cohort"),
        ("GSE92415", "anti-TNF", "pretreatment mucosa", "include", "independent trial cohort"),
        ("E-MTAB-7604", "anti-TNF", "pretreatment mucosa", "include", "independent RNA-seq cohort"),
        ("GSE14580", "anti-TNF", "pretreatment mucosa", "exclude", "UC subset represented in GSE16879"),
        ("GSE12251", "anti-TNF", "pretreatment mucosa", "exclude", "expression-level subset of GSE23597"),
        ("GSE73661_IFX", "anti-TNF", "pretreatment mucosa", "exclude", "independence from earlier Leuven IFX series not established"),
        ("GSE73661_VDZ", "anti-integrin", "pretreatment mucosa", "boundary", "distinct vedolizumab trial population"),
        ("GSE206285_UST", "anti-IL12/23", "pretreatment mucosa", "boundary", "distinct UNIFI trial population"),
    ]
    files = {
        x: ROOT / f"data/bulk_geo/{x.split('_')[0]}_series_matrix.txt.gz"
        for x, *_ in specs if x != "E-MTAB-7604"
    }
    for accession, therapy, compartment, decision, reason in specs:
        path = files.get(accession)
        records.append({"dataset": accession, "therapy_class": therapy, "compartment": compartment,
                        "role": decision, "independence_rationale": reason,
                   "source_file": path.relative_to(ROOT).as_posix() if path else "E-MTAB-7604 raw counts/SDRF",
                        "sha256": sha256(path) if path and path.exists() else "not_local_or_incomplete"})
    pd.DataFrame(records).to_csv(OUT / "dataset_independence_evidence_register.csv", index=False)


def main():
    frames, genes, availability = load_antitnf_cohorts()
    availability.to_csv(OUT / "benchmark_signature_gene_availability.csv", index=False)
    effects, studies, meta, ratings = benchmark(frames)
    subgroup_analysis(frames)
    clinical_increment(frames["GSE92415"])
    fingerprint = expression_overlap_audit()
    independence_audit(fingerprint)
    vdz = load_vdz_boundary()
    ust = load_ust_boundary()
    boundary_analysis(vdz, ust)

    registry = []
    for name, geneset in BENCHMARKS.items():
        registry.append({"signature": name, "n_genes": len(geneset), "genes": ";".join(geneset),
                         "expected_high_group": "non-response",
                         "development_overlap_accessions": ";".join(sorted(DEVELOPMENT_OVERLAP[name]))})
    pd.DataFrame(registry).to_csv(OUT / "benchmark_signature_registry.csv", index=False)
    manifest = {
        "analysis_date": "2026-08-14",
        "anti_tnf_cohorts": list(frames),
        "benchmark_signatures": list(BENCHMARKS),
        "boundary_datasets_completed": ["GSE73661_VDZ"] + (["GSE206285_UST"] if ust is not None else []),
        "GSE206285_status": "complete" if ust is not None else "download_or_parse_incomplete",
        "outputs": sorted(p.name for p in OUT.glob("*")),
    }
    (OUT / "upgrade_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
