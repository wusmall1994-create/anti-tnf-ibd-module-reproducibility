from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "module_scores"
OUT_DIR.mkdir(parents=True, exist_ok=True)


MODULES: dict[str, list[str]] = {
    "myeloid_inflammation": [
        "IL1B",
        "TNF",
        "CXCL8",
        "SPP1",
        "OSM",
        "IL1RN",
        "CCL3",
        "CCL4",
        "NFKBIA",
    ],
    "epithelial_ifn_damage": [
        "ISG15",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "IFI6",
        "MX1",
        "CXCL10",
        "HLA-DRA",
        "HLA-DRB1",
    ],
    "epithelial_barrier_maturity": [
        "MUC2",
        "TFF3",
        "KRT20",
        "CA1",
        "CA2",
        "BEST4",
        "OTOP2",
        "SLC26A3",
    ],
    "stromal_fibroinflammatory": [
        "COL1A1",
        "COL3A1",
        "DCN",
        "LUM",
        "CXCL12",
        "IL11",
        "OSMR",
        "CCL19",
        "CXCL13",
        "WNT2B",
        "RSPO3",
    ],
}


def load_sample_metadata() -> pd.DataFrame:
    meta = pd.read_csv(ROOT / "results" / "feasibility" / "GSE282122_sample_metadata.csv")
    meta["sample_key"] = meta["title"].str.replace("-reup", "", regex=False)
    return meta


def load_outcomes() -> pd.DataFrame:
    paired = pd.read_csv(ROOT / "data" / "metadata" / "paired_sample_list.csv", sep="\t")
    paired["sample_key"] = paired["sample_id"].str.replace("-reup", "", regex=False)
    return paired[["sample_key", "Category"]]


def find_first_existing(columns: Iterable[str], candidates: list[str]) -> str | None:
    cols = set(columns)
    for c in candidates:
        if c in cols:
            return c
    lowered = {c.lower(): c for c in columns}
    for c in candidates:
        if c.lower() in lowered:
            return lowered[c.lower()]
    return None


def normalize_sample_id(value: str) -> str:
    value = str(value)
    return value.replace("-reup", "")


def barcode_to_sample(index: pd.Index) -> pd.Series:
    vals = []
    for cell_id in index.astype(str):
        parts = cell_id.split("-")
        vals.append("-".join(parts[-2:]) if len(parts) >= 4 else "")
    return pd.Series(vals, index=index, dtype="object")


def choose_expression_matrix(adata: ad.AnnData, layer: str | None):
    if layer:
        if layer not in adata.layers:
            raise SystemExit(f"Layer {layer!r} not found. Available layers: {list(adata.layers)}")
        return adata.layers[layer]
    return adata.X


def mean_by_group(matrix, group_codes: np.ndarray, n_groups: int) -> np.ndarray:
    rows = []
    for i in range(n_groups):
        mask = group_codes == i
        if mask.sum() == 0:
            rows.append(np.nan)
            continue
        sub = matrix[mask]
        if sparse.issparse(sub):
            rows.append(np.asarray(sub.mean(axis=0)).ravel())
        else:
            rows.append(np.asarray(sub).mean(axis=0))
    return np.vstack(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, help="Input AnnData/H5AD file")
    parser.add_argument("--label", required=True, help="Short label, e.g. myeloid")
    parser.add_argument("--layer", default=None, help="Optional expression layer")
    parser.add_argument(
        "--sample-column",
        default=None,
        help="Optional obs column containing sample IDs. If omitted, script guesses.",
    )
    parser.add_argument(
        "--celltype-column",
        default=None,
        help="Optional obs column for cell-state annotation to carry into summaries.",
    )
    args = parser.parse_args()

    h5ad = Path(args.h5ad)
    adata = ad.read_h5ad(h5ad)
    var_names = pd.Index(adata.var_names.astype(str).str.upper())
    obs = adata.obs.copy()
    obs["_row_index"] = np.arange(obs.shape[0])

    sample_col = args.sample_column or find_first_existing(
        obs.columns.tolist(), ["sample_id", "SampleID2", "SampleID3", "sample", "Sample", "orig.ident"]
    )
    if sample_col:
        obs["sample_key"] = obs[sample_col].map(normalize_sample_id)
    else:
        obs["sample_key"] = barcode_to_sample(obs.index)

    celltype_col = args.celltype_column or find_first_existing(
        obs.columns.tolist(), ["final_analysis", "minor", "major", "sub_bucket", "bucket", "cell_type"]
    )

    sample_meta = load_sample_metadata()
    outcomes = load_outcomes()
    obs = obs.merge(sample_meta, on="sample_key", how="left", suffixes=("", "_geo"))
    obs = obs.merge(outcomes, on="sample_key", how="left")

    X = choose_expression_matrix(adata, args.layer)
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    else:
        X = X.tocsr()

    present_modules = {}
    for module, genes in MODULES.items():
        idx = [i for i, g in enumerate(var_names) if g in set(genes)]
        present_modules[module] = {
            "requested_genes": genes,
            "present_genes": [var_names[i] for i in idx],
            "indices": idx,
        }

    sample_codes, sample_levels = pd.factorize(obs["sample_key"], sort=True)
    sample_summary = pd.DataFrame({"sample_key": sample_levels})
    sample_summary["n_cells"] = np.bincount(sample_codes, minlength=len(sample_levels))

    for col in ["disease", "treatment", "inflammation", "site", "patient", "batch", "Category"]:
        if col in obs.columns:
            sample_summary[col] = [
                obs.loc[obs["sample_key"].eq(s), col].dropna().astype(str).mode().iloc[0]
                if obs.loc[obs["sample_key"].eq(s), col].dropna().shape[0]
                else ""
                for s in sample_levels
            ]

    module_matrix = {}
    for module, spec in present_modules.items():
        idx = spec["indices"]
        if not idx:
            sample_summary[module] = np.nan
            continue
        vals = np.asarray(X[:, idx].mean(axis=1)).ravel()
        sums = np.bincount(sample_codes, weights=vals, minlength=len(sample_levels))
        counts = np.bincount(sample_codes, minlength=len(sample_levels))
        sample_summary[module] = sums / counts
        module_matrix[module] = vals

    out_csv = OUT_DIR / f"{args.label}_sample_module_scores.csv"
    sample_summary.to_csv(out_csv, index=False)

    if celltype_col and celltype_col in obs.columns:
        ct_rows = []
        obs["_celltype_for_summary"] = obs[celltype_col].astype(str)
        group = obs.groupby(["sample_key", "_celltype_for_summary"], observed=True)
        for (sample, celltype), idx_labels in group.groups.items():
            iloc = obs.loc[idx_labels, "_row_index"].to_numpy(dtype=int)
            row = {"sample_key": sample, "celltype": celltype, "n_cells": len(iloc)}
            for module, vals in module_matrix.items():
                row[module] = float(np.mean(vals[iloc])) if len(iloc) else np.nan
            ct_rows.append(row)
        ct_df = pd.DataFrame(ct_rows)
        out_ct = OUT_DIR / f"{args.label}_sample_celltype_module_scores.csv"
        ct_df.to_csv(out_ct, index=False)
    else:
        out_ct = None

    audit = {
        "input": str(h5ad),
        "label": args.label,
        "shape": list(adata.shape),
        "obs_columns": obs.columns.tolist(),
        "sample_column_used": sample_col or "barcode_suffix",
        "celltype_column_used": celltype_col,
        "output_sample_scores": str(out_csv),
        "output_celltype_scores": str(out_ct) if out_ct else None,
        "modules": {
            k: {
                "requested_genes": v["requested_genes"],
                "present_genes": v["present_genes"],
                "n_present": len(v["present_genes"]),
            }
            for k, v in present_modules.items()
        },
    }
    audit_path = OUT_DIR / f"{args.label}_module_score_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
