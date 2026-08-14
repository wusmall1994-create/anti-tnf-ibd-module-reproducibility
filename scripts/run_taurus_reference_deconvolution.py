from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

ROOT = Path(
    os.environ.get("ANTI_TNF_IBD_ROOT", Path(__file__).resolve().parents[1])
).resolve()
OUT = ROOT / "results" / "calibration_upgrade_20260814"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_calibration_upgrade import ALL, BENCHMARKS, MODULES, hedges_g, load_antitnf_cohorts

H5ADS = {
    "myeloid": ROOT / "data/taurus_zenodo/myeloid_final.h5ad",
    "stroma": ROOT / "data/taurus_zenodo/fibperi_final.h5ad",
    "epithelium": ROOT / "data/taurus_zenodo/epicolonic_final.h5ad",
}
STATES = ["resident_like_macrophage", "inflammatory_monocyte", "other_myeloid", "stroma", "colonic_epithelium"]


def eligible(obs: pd.DataFrame) -> np.ndarray:
    pre = obs["Treatment"].astype(str).str.lower().eq("pre")
    inflamed = obs["Inflammation"].astype(str).str.lower().eq("inflamed")
    score = pd.to_numeric(obs["Inflammation_score"], errors="coerce") > 6.5
    return (pre & inflamed & score).to_numpy()


def aggregate_file(path: Path, compartment: str, block=5000):
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs
    keep = eligible(obs)
    if compartment == "myeloid":
        labels = obs["final_analysis"].astype(str)
        state_masks = {
            "resident_like_macrophage": keep & labels.eq("C1Qhi IL1Blo macro").to_numpy(),
            "inflammatory_monocyte": keep & labels.isin(["S100A8 A9hi mono", "S100A8 A9hi TNFhi IL6pos mono"]).to_numpy(),
            "other_myeloid": keep & ~labels.isin(["C1Qhi IL1Blo macro", "S100A8 A9hi mono", "S100A8 A9hi TNFhi IL6pos mono"]).to_numpy(),
        }
    else:
        state_name = "colonic_epithelium" if compartment == "epithelium" else compartment
        state_masks = {state_name: keep}
    sums = {s: np.zeros(a.n_vars, dtype=np.float64) for s in state_masks}
    counts = {s: int(m.sum()) for s, m in state_masks.items()}
    for start in range(0, a.n_obs, block):
        end = min(a.n_obs, start + block)
        x = a.X[start:end]
        for state, mask in state_masks.items():
            local = mask[start:end]
            if local.any():
                sums[state] += np.asarray(x[local].sum(axis=0)).ravel()
    genes = a.var["gene_symbol"].astype(str).to_numpy()
    a.file.close()
    return genes, sums, counts


def build_reference():
    cache = OUT / "TAURUS_partial_reference_logCPM.csv"
    if cache.exists():
        cached = pd.read_csv(cache, index_col=0)
        if list(cached.columns) == STATES and not cached.isna().any().any():
            return cached, pd.read_csv(OUT / "TAURUS_reference_state_cell_counts.csv")
    pieces, count_rows = [], []
    for compartment, path in H5ADS.items():
        genes, sums, counts = aggregate_file(path, compartment)
        frame = pd.DataFrame(sums, index=genes).groupby(level=0).sum()
        pieces.append(frame)
        count_rows.extend({"state": s, "n_cells": n, "source_object": compartment} for s, n in counts.items())
    raw = pd.concat(pieces, axis=1).fillna(0).reindex(columns=STATES)
    ref = np.log2(raw.div(raw.sum(axis=0), axis=1) * 1e6 + 1)
    ref.to_csv(cache)
    counts = pd.DataFrame(count_rows)
    counts.to_csv(OUT / "TAURUS_reference_state_cell_counts.csv", index=False)
    return ref, counts


def choose_markers(ref: pd.DataFrame, n_each: int):
    excluded = {g.upper() for genes in list(MODULES.values()) + list(BENCHMARKS.values()) for g in genes}
    candidates = ref.loc[[g for g in ref.index if g.upper() not in excluded and not g.upper().startswith(("RPL", "RPS", "MT-"))]]
    marker_rows = []
    for state in STATES:
        other = candidates.drop(columns=state).max(axis=1)
        specificity = candidates[state] - other
        order = specificity[candidates[state] >= 2].sort_values(ascending=False).head(n_each)
        marker_rows.extend({"state": state, "gene": g, "specificity_log2cpm": v, "n_each": n_each} for g, v in order.items())
    return pd.DataFrame(marker_rows)


def rank_columns(frame: pd.DataFrame):
    return pd.DataFrame({c: rankdata(frame[c].to_numpy()) / (len(frame) + 1) for c in frame}, index=frame.index)


def estimate(expr: pd.DataFrame, ref: pd.DataFrame, markers: list[str]):
    genes = sorted(set(markers) & set(expr.index) & set(ref.index))
    r = rank_columns(ref.loc[genes, STATES]).to_numpy()
    out = []
    penalty = 5.0
    design = np.vstack([r, np.ones((1, len(STATES))) * penalty])
    for sample in expr.columns:
        y = rankdata(expr.loc[genes, sample].to_numpy()) / (len(genes) + 1)
        coef, _ = nnls(design, np.r_[y, penalty])
        coef = coef / coef.sum() if coef.sum() else coef
        out.append(coef)
    return pd.DataFrame(out, index=expr.columns, columns=STATES), len(genes)


def synthetic_validation(ref, markers, n=500):
    rng = np.random.default_rng(20260814)
    genes = sorted(set(markers) & set(ref.index))
    r = ref.loc[genes, STATES]
    truth = rng.dirichlet(np.ones(len(STATES)), size=n)
    rows = []
    est = []
    for w in truth:
        pseudo = pd.DataFrame((r.to_numpy() @ w)[:, None], index=genes, columns=["pseudo"])
        e, _ = estimate(pseudo, ref, genes)
        est.append(e.iloc[0].to_numpy())
    est = np.asarray(est)
    for j, state in enumerate(STATES):
        rho = spearmanr(truth[:, j], est[:, j]).statistic
        rmse = np.sqrt(np.mean((truth[:, j] - est[:, j]) ** 2))
        rows.append({"state": state, "spearman_rho": rho, "rmse": rmse, "n_simulations": n})
    return pd.DataFrame(rows)


def associations(frames, estimates):
    rows, atten = [], []
    for cohort, d in frames.items():
        frac = estimates[cohort]
        x = d.merge(frac, left_on="geo_accession", right_index=True, how="inner")
        for stratum, s in x.groupby("stratum"):
            if s.response_binary.nunique() < 2:
                continue
            for state in STATES:
                yes = s.loc[s.response_binary == 1, state].to_numpy()
                no = s.loc[s.response_binary == 0, state].to_numpy()
                if min(len(yes), len(no)) >= 3:
                    g, var = hedges_g(yes, no)
                    rows.append({"cohort": cohort, "stratum": stratum, "state": state, "n_response": len(yes), "n_nonresponse": len(no), "hedges_g_response_minus_nonresponse": g, "variance": var})
            cols = ["myeloid_inflammation", "resident_like_macrophage", "inflammatory_monocyte", "response_binary"]
            z = s[cols].dropna()
            if len(z) >= 20 and z.response_binary.nunique() == 2:
                scaler = StandardScaler()
                m0 = LogisticRegression(C=1.0, solver="liblinear").fit(scaler.fit_transform(z[["myeloid_inflammation"]]), z.response_binary)
                m1 = LogisticRegression(C=1.0, solver="liblinear").fit(scaler.fit_transform(z[["myeloid_inflammation", "resident_like_macrophage", "inflammatory_monocyte"]]), z.response_binary)
                b0, b1 = float(m0.coef_[0, 0]), float(m1.coef_[0, 0])
                atten.append({"cohort": cohort, "stratum": stratum, "n": len(z), "myeloid_beta_unadjusted": b0, "myeloid_beta_composition_adjusted": b1, "absolute_beta_change": abs(b0) - abs(b1), "attenuation_fraction": (abs(b0) - abs(b1)) / abs(b0) if b0 else np.nan})
    return pd.DataFrame(rows), pd.DataFrame(atten)


def main():
    ref, counts = build_reference()
    frames, genes, _ = load_antitnf_cohorts()
    all_estimates, availability = {}, []
    for n_each in (25, 50):
        markers = choose_markers(ref, n_each)
        markers.to_csv(OUT / f"TAURUS_reference_markers_top{n_each}.csv", index=False)
        synthetic_validation(ref, markers.gene.tolist()).to_csv(OUT / f"TAURUS_synthetic_validation_top{n_each}.csv", index=False)
        for cohort, expr in genes.items():
            est, used = estimate(expr, ref, markers.gene.tolist())
            est.index.name = "geo_accession"
            est.reset_index().to_csv(OUT / f"{cohort}_TAURUS_partial_reference_top{n_each}.csv", index=False)
            availability.append({"cohort": cohort, "marker_set": f"top{n_each}", "n_markers_defined": markers.gene.nunique(), "n_markers_used": used})
            if n_each == 50:
                all_estimates[cohort] = est
    assoc, atten = associations(frames, all_estimates)
    assoc.to_csv(OUT / "TAURUS_deconvolution_outcome_associations.csv", index=False)
    atten.to_csv(OUT / "TAURUS_myeloid_module_composition_adjustment.csv", index=False)
    pd.DataFrame(availability).to_csv(OUT / "TAURUS_deconvolution_marker_availability.csv", index=False)
    manifest = {"reference_scope": "TAURUS pretreatment, inflamed, inflammation score >6.5; partial five-state reference", "states": STATES, "target_signatures_excluded_from_marker_selection": list(ALL), "primary_marker_set": "top50_per_state", "sensitivity_marker_set": "top25_per_state", "cell_counts": counts.to_dict(orient="records")}
    (OUT / "TAURUS_deconvolution_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
