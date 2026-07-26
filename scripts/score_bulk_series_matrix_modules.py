from __future__ import annotations

import argparse
import gzip
import json
import re
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "bulk_validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)
GPL_DIR = ROOT / "data" / "bulk_geo" / "platforms"
GPL_DIR.mkdir(parents=True, exist_ok=True)


MODULES: dict[str, list[str]] = {
    "myeloid_inflammation": ["IL1B", "TNF", "CXCL8", "SPP1", "OSM", "IL1RN", "CCL3", "CCL4", "NFKBIA"],
    "epithelial_ifn_damage": ["ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "MX1", "CXCL10", "HLA-DRA", "HLA-DRB1"],
    "epithelial_barrier_maturity": ["MUC2", "TFF3", "KRT20", "CA1", "CA2", "BEST4", "OTOP2", "SLC26A3"],
    "stromal_fibroinflammatory": ["COL1A1", "COL3A1", "DCN", "LUM", "CXCL12", "IL11", "OSMR", "CCL19", "CXCL13", "WNT2B", "RSPO3"],
}


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    sample_meta: dict[str, list[str]] = {}
    expr_lines: list[str] = []
    in_table = False
    platforms: list[str] = []
    with open_text(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                expr_lines.append(line)
                continue
            if line.startswith("!Sample_"):
                parts = [x.strip('"') for x in line.split("\t")]
                key = parts[0].replace("!Sample_", "")
                if key in sample_meta:
                    i = 2
                    new_key = f"{key}_{i}"
                    while new_key in sample_meta:
                        i += 1
                        new_key = f"{key}_{i}"
                    key = new_key
                sample_meta[key] = parts[1:]
                if key == "platform_id":
                    platforms = parts[1:]

    if not expr_lines:
        raise ValueError(f"No expression table found in {path}")
    header = [x.strip('"') for x in expr_lines[0].split("\t")]
    rows = [[x.strip('"') for x in line.split("\t")] for line in expr_lines[1:] if line]
    expr = pd.DataFrame(rows, columns=header)
    expr = expr.set_index("ID_REF")
    expr = expr.apply(pd.to_numeric, errors="coerce")

    accessions = sample_meta.get("geo_accession", expr.columns.tolist())
    meta = pd.DataFrame({"geo_accession": accessions})
    for key, vals in sample_meta.items():
        if len(vals) == len(accessions):
            meta[key] = vals
    return expr, meta, sorted(set(platforms))


def platform_prefix(gpl: str) -> str:
    m = re.match(r"GPL(\d+)", gpl.upper())
    if not m:
        raise ValueError(f"Not a GPL accession: {gpl}")
    digits = m.group(1)
    return f"GPL{digits[:-3]}nnn"


def platform_annot_url(gpl: str) -> str:
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{platform_prefix(gpl)}/{gpl.upper()}/annot/"
        f"{gpl.upper()}.annot.gz"
    )


def download_file(url: str, out: Path) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r, out.open("wb") as f:
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)


def extract_gene_symbol(value: str) -> str | None:
    if not isinstance(value, str):
        return None
    value = value.strip()
    if not value or value in {"---", "NA", "nan"}:
        return None
    # Affymetrix annotations often use "GENE /// GENE2"; keep the first.
    value = re.split(r"\s*///\s*|\s*;\s*|\s*,\s*", value)[0].strip()
    # Some platform fields include "NM_x // SYMBOL // description"; find likely symbols.
    if " // " in value:
        parts = [p.strip() for p in value.split(" // ")]
        candidates = [p for p in parts if re.fullmatch(r"[A-Z0-9][A-Z0-9.-]{1,15}", p)]
        if candidates:
            value = candidates[0]
    if not re.search(r"[A-Za-z]", value):
        return None
    return value.upper()


def load_platform_mapping(gpl: str, download: bool) -> dict[str, str]:
    annot = GPL_DIR / f"{gpl.upper()}.annot.gz"
    if download and not annot.exists():
        download_file(platform_annot_url(gpl), annot)
    if not annot.exists():
        raise FileNotFoundError(f"Missing platform annotation {annot}. Use --download-platform.")

    # Find table between !platform_table_begin/end.
    lines = []
    in_table = False
    with gzip.open(annot, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                lines.append(line)
    if not lines:
        raise ValueError(f"No platform table found in {annot}")
    header = lines[0].split("\t")
    candidates = [
        "Gene Symbol",
        "Gene symbol",
        "GENE_SYMBOL",
        "Symbol",
        "gene_assignment",
        "Gene assignment",
    ]
    symbol_col = next((c for c in candidates if c in header), None)
    if symbol_col is None:
        # Loose fallback.
        symbol_col = next((c for c in header if "symbol" in c.lower() or "gene_assignment" in c.lower()), None)
    if symbol_col is None:
        raise ValueError(f"No gene symbol-like column found in {annot}; columns={header[:30]}")
    id_col = "ID" if "ID" in header else header[0]

    idx_id = header.index(id_col)
    idx_symbol = header.index(symbol_col)
    mapping = {}
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) <= max(idx_id, idx_symbol):
            continue
        probe = parts[idx_id]
        sym = extract_gene_symbol(parts[idx_symbol])
        if probe and sym:
            mapping[probe] = sym
    return mapping


def collapse_to_gene(expr: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    mapped = expr.loc[expr.index.intersection(mapping.keys())].copy()
    mapped["gene_symbol"] = [mapping[i] for i in mapped.index]
    return mapped.groupby("gene_symbol").mean(numeric_only=True)


def score_modules(gene_expr: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    rows = pd.DataFrame({"geo_accession": gene_expr.columns})
    audit = {}
    for module, genes in MODULES.items():
        present = [g for g in genes if g in gene_expr.index]
        audit[module] = {"requested_genes": genes, "present_genes": present, "n_present": len(present)}
        rows[module] = gene_expr.loc[present].mean(axis=0).values if present else np.nan
    return rows, audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--series-matrix", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--platform", default=None, help="GPL accession; guessed if omitted and unique.")
    parser.add_argument("--download-platform", action="store_true")
    args = parser.parse_args()

    matrix = Path(args.series_matrix)
    expr, meta, platforms = parse_series_matrix(matrix)
    platform = args.platform or (platforms[0] if len(platforms) == 1 else None)
    if not platform:
        raise SystemExit(f"Could not infer one platform from {platforms}; provide --platform GPLxxxx")

    mapping = load_platform_mapping(platform, args.download_platform)
    gene_expr = collapse_to_gene(expr, mapping)
    scores, audit = score_modules(gene_expr)
    scores = scores.merge(meta, on="geo_accession", how="left")

    out_csv = OUT_DIR / f"{args.label}_bulk_module_scores.csv"
    out_audit = OUT_DIR / f"{args.label}_bulk_module_score_audit.json"
    scores.to_csv(out_csv, index=False)
    out_audit.write_text(
        json.dumps(
            {
                "series_matrix": str(matrix),
                "platform": platform,
                "n_probes": int(expr.shape[0]),
                "n_samples": int(expr.shape[1]),
                "n_mapped_genes": int(gene_expr.shape[0]),
                "module_audit": audit,
                "output": str(out_csv),
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    print(out_audit.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
