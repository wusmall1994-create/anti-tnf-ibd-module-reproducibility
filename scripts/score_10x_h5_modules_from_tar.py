from __future__ import annotations

import argparse
import csv
import json
import shutil
import tarfile
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TAR = Path("D:/TAURUS_data/GSE282122_filtered_processed_data.tar.gz")
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


def read_10x_h5_score(path: Path) -> tuple[dict[str, float | int | str], dict[str, dict[str, object]]]:
    with h5py.File(path, "r") as h:
        g = h["matrix"]
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        shape = tuple(g["shape"][:])
        genes = np.array([x.decode("utf-8") for x in g["features"]["name"][:]])
        x = sparse.csc_matrix((data, indices, indptr), shape=shape, dtype=np.float64)

    lib = np.asarray(x.sum(axis=0)).ravel()
    valid = lib > 0
    scale = np.zeros_like(lib, dtype=np.float64)
    scale[valid] = 1e4 / lib[valid]

    gene_to_idx = {g.upper(): i for i, g in enumerate(genes.astype(str))}
    row: dict[str, float | int | str] = {
        "n_genes": int(shape[0]),
        "n_cells": int(shape[1]),
        "median_counts_per_cell": float(np.median(lib[valid])) if valid.any() else np.nan,
        "total_counts": float(lib.sum()),
    }
    audit: dict[str, dict[str, object]] = {}
    for module, module_genes in MODULES.items():
        idx = [gene_to_idx[g] for g in module_genes if g in gene_to_idx]
        present = [g for g in module_genes if g in gene_to_idx]
        audit[module] = {
            "requested_genes": module_genes,
            "present_genes": present,
            "n_present": len(present),
        }
        if not idx:
            row[module] = np.nan
            row[module + "_detected_fraction"] = np.nan
            continue
        sub = x[idx, :].tocsc(copy=True)
        sub = sub.multiply(scale)
        # log1p only non-zero values; zeros stay zero and are included in sparse mean.
        sub.data = np.log1p(sub.data)
        row[module] = float(sub.mean())
        detected_cells = np.asarray((x[idx, :] > 0).sum(axis=0)).ravel() > 0
        row[module + "_detected_fraction"] = float(detected_cells.mean())
    return row, audit


def sample_from_member_name(name: str) -> str:
    parts = Path(name).parts
    # filtered_processed_data/CID003376-1/filtered_feature_bc_matrix.h5
    if len(parts) >= 2:
        return parts[-2].replace("-reup", "")
    return ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tar", default=str(DEFAULT_TAR))
    parser.add_argument("--label", default="geo_10x")
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Process complete members from a still-downloading tar.gz until EOF.",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Optional limit for smoke tests.",
    )
    args = parser.parse_args()

    tar_path = Path(args.tar)
    if not tar_path.exists():
        raise SystemExit(f"Missing tar: {tar_path}")

    rows = []
    audits: dict[str, object] = {}
    stopped_partial = False
    with tempfile.TemporaryDirectory(prefix="taurus_10x_") as td:
        tmpdir = Path(td)
        try:
            with tarfile.open(tar_path, "r:gz") as tar:
                for member in tar:
                    if not (member.isfile() and member.name.endswith("filtered_feature_bc_matrix.h5")):
                        continue
                    sample = sample_from_member_name(member.name)
                    tmp_h5 = tmpdir / f"{sample}.h5"
                    src = tar.extractfile(member)
                    if src is None:
                        continue
                    with tmp_h5.open("wb") as out:
                        shutil.copyfileobj(src, out)
                    row, audit = read_10x_h5_score(tmp_h5)
                    row["sample_key"] = sample
                    row["source_member"] = member.name
                    rows.append(row)
                    audits[sample] = audit
                    print(f"scored {sample} ({len(rows)})", flush=True)
                    tmp_h5.unlink(missing_ok=True)
                    if args.max_samples and len(rows) >= args.max_samples:
                        break
        except (EOFError, tarfile.ReadError) as e:
            if not args.allow_partial:
                raise
            stopped_partial = True
            print(f"PARTIAL_TAR_STOP: {type(e).__name__}: {e}", flush=True)

    scores = pd.DataFrame(rows)
    meta = load_sample_metadata()
    outcomes = load_outcomes()
    if not scores.empty:
        scores = scores.merge(meta, on="sample_key", how="left")
        scores = scores.merge(outcomes, on="sample_key", how="left")
        scores["Category"] = scores["Category"].fillna("No_outcome")

    suffix = "_partial" if stopped_partial else ""
    out_csv = OUT_DIR / f"{args.label}_sample_module_scores{suffix}.csv"
    scores.to_csv(out_csv, index=False)
    audit_path = OUT_DIR / f"{args.label}_10x_module_score_audit{suffix}.json"
    audit_doc = {
        "tar": str(tar_path),
        "n_samples_scored": int(scores.shape[0]),
        "stopped_partial": stopped_partial,
        "output": str(out_csv),
        "module_gene_audit_by_sample": audits,
    }
    audit_path.write_text(json.dumps(audit_doc, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({k: v for k, v in audit_doc.items() if k != "module_gene_audit_by_sample"}, indent=2))


if __name__ == "__main__":
    main()
