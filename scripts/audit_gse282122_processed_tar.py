from __future__ import annotations

import csv
import gzip
import json
import tarfile
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path("D:/TAURUS_data")
TAR_PATH = DATA / "GSE282122_filtered_processed_data.tar.gz"
OUT_DIR = ROOT / "results" / "feasibility"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def load_sample_metadata() -> dict[str, dict[str, str]]:
    meta_path = OUT_DIR / "GSE282122_sample_metadata.csv"
    with meta_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    # title is like CID003352-2-reup; cell barcodes and Zenodo paired list often omit -reup.
    out: dict[str, dict[str, str]] = {}
    for r in rows:
        title = r.get("title", "")
        out[title] = r
        out[title.replace("-reup", "")] = r
    return out


def load_outcomes() -> dict[str, str]:
    path = ROOT / "data" / "metadata" / "paired_sample_list.csv"
    out: dict[str, str] = {}
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            sid = r["sample_id"]
            out[sid] = r["Category"]
            out[sid.replace("-reup", "")] = r["Category"]
    return out


def main() -> None:
    if not TAR_PATH.exists():
        raise SystemExit(f"Missing {TAR_PATH}")
    if TAR_PATH.stat().st_size < 1_000_000_000:
        raise SystemExit(
            f"{TAR_PATH} is only {TAR_PATH.stat().st_size:,} bytes; download likely incomplete."
        )

    meta = load_sample_metadata()
    outcomes = load_outcomes()
    members_summary = []
    sample_counter = Counter()
    h5_samples = set()
    suffix_counter = Counter()

    with tarfile.open(TAR_PATH, "r:gz") as tar:
        members = tar.getmembers()
        for m in members:
            p = Path(m.name)
            suffix = "".join(p.suffixes[-2:]) if p.suffixes[-2:] else p.suffix
            suffix_counter[suffix] += 1
            # Most 10x processed archives are nested by sample id.
            parts = p.parts
            for part in parts:
                key = part.replace("-reup", "")
                if key in meta:
                    sample_counter[key] += 1
                    if m.isfile() and m.name.endswith("filtered_feature_bc_matrix.h5"):
                        h5_samples.add(key)
                    break
            members_summary.append(
                {
                    "name": m.name,
                    "size": m.size,
                    "suffix": suffix,
                    "is_file": m.isfile(),
                }
            )

    matched_samples = sorted(sample_counter)
    disease_counts = Counter(meta[s].get("disease", "") for s in matched_samples)
    category_counts = Counter(outcomes.get(s, "No_outcome") for s in matched_samples)

    report = {
        "tar_path": str(TAR_PATH),
        "tar_size_bytes": TAR_PATH.stat().st_size,
        "n_members": len(members_summary),
        "suffix_counts": dict(suffix_counter.most_common()),
        "n_matched_samples": len(matched_samples),
        "n_10x_h5_samples": len(h5_samples),
        "disease_counts_among_matched": dict(disease_counts),
        "outcome_counts_among_matched": dict(category_counts),
        "h5_sample_examples": sorted(h5_samples)[:20],
        "example_members": members_summary[:30],
    }

    (OUT_DIR / "GSE282122_processed_tar_audit.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    with (OUT_DIR / "GSE282122_processed_tar_members.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(f, fieldnames=["name", "size", "suffix", "is_file"])
        writer.writeheader()
        writer.writerows(members_summary)

    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
