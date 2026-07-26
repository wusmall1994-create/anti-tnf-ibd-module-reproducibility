from __future__ import annotations

import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path("D:/TAURUS_data")
UMAP = DATA / "UMAP_combined_objects.txt.gz"
META = ROOT / "results" / "feasibility" / "GSE282122_sample_metadata.csv"
PAIRED = ROOT / "data" / "metadata" / "paired_sample_list.csv"
OUT_DIR = ROOT / "results" / "feasibility"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def barcode_to_sample(cell_id: str) -> str:
    # Example: AAACCCAAGATCCCAT-1-CID003352-2
    parts = cell_id.strip().split("-")
    if len(parts) < 4:
        return ""
    return "-".join(parts[-2:])


def load_meta() -> dict[str, dict[str, str]]:
    with META.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for r in rows:
        title = r.get("title", "")
        out[title] = r
        out[title.replace("-reup", "")] = r
    return out


def load_categories() -> dict[str, str]:
    out = {}
    with PAIRED.open(encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sid = r["sample_id"]
            out[sid] = r["Category"]
            out[sid.replace("-reup", "")] = r["Category"]
    return out


def main() -> None:
    if not UMAP.exists():
        raise SystemExit(f"Missing {UMAP}")

    meta = load_meta()
    cats = load_categories()
    sample_cells = Counter()
    n_cells = 0
    with gzip.open(UMAP, "rt", encoding="utf-8", errors="replace") as f:
        header = next(f)
        for line in f:
            if not line.strip():
                continue
            cell = line.split("\t", 1)[0]
            sid = barcode_to_sample(cell)
            sample_cells[sid] += 1
            n_cells += 1

    rows = []
    for sid, n in sample_cells.most_common():
        m = meta.get(sid, {})
        rows.append(
            {
                "sample_id": sid,
                "n_cells_in_umap": n,
                "category": cats.get(sid, "No_outcome"),
                "disease": m.get("disease", ""),
                "treatment": m.get("treatment", ""),
                "inflammation": m.get("inflammation", ""),
                "site": m.get("site", ""),
                "patient": m.get("patient", ""),
                "batch": m.get("batch", ""),
            }
        )

    with (OUT_DIR / "GSE282122_umap_sample_coverage.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    by_category = defaultdict(lambda: {"samples": 0, "cells": 0, "patients": set()})
    by_disease_treatment = Counter()
    unmatched = 0
    for r in rows:
        if not r["disease"]:
            unmatched += 1
        cat = r["category"]
        by_category[cat]["samples"] += 1
        by_category[cat]["cells"] += int(r["n_cells_in_umap"])
        if r["patient"]:
            by_category[cat]["patients"].add(r["patient"])
        by_disease_treatment[(r["disease"], r["treatment"])] += 1

    summary = {
        "umap_file": str(UMAP),
        "n_cells": n_cells,
        "n_samples": len(sample_cells),
        "unmatched_samples": unmatched,
        "category_summary": {
            k: {
                "samples": v["samples"],
                "cells": v["cells"],
                "patients": len(v["patients"]),
            }
            for k, v in sorted(by_category.items())
        },
        "disease_treatment_sample_counts": {
            f"{k[0]}|{k[1]}": v for k, v in sorted(by_disease_treatment.items())
        },
        "min_cells_per_sample": min(sample_cells.values()),
        "median_cells_per_sample": sorted(sample_cells.values())[len(sample_cells) // 2],
        "max_cells_per_sample": max(sample_cells.values()),
    }
    (OUT_DIR / "GSE282122_umap_sample_coverage_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
