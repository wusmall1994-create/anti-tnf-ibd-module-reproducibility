from __future__ import annotations

import argparse
import gzip
import json
import re
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data" / "bulk_geo"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def series_prefix(gse: str) -> str:
    # GSE16879 -> GSE16nnn, GSE12251 -> GSE12nnn
    m = re.match(r"GSE(\d+)", gse.upper())
    if not m:
        raise ValueError(f"Not a GSE accession: {gse}")
    digits = m.group(1)
    return f"GSE{digits[:-3]}nnn"


def candidate_url(gse: str) -> str:
    prefix = series_prefix(gse)
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse.upper()}/matrix/"
        f"{gse.upper()}_series_matrix.txt.gz"
    )


def download(url: str, out: Path) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r, out.open("wb") as f:
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)


def audit_matrix(path: Path) -> dict[str, object]:
    n_samples = 0
    title = ""
    characteristics: list[str] = []
    platforms: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!Series_title"):
                title = line.strip().split("\t", 1)[-1].strip('"')
            elif line.startswith("!Sample_geo_accession"):
                n_samples = max(0, len(line.rstrip("\n").split("\t")) - 1)
            elif line.startswith("!Sample_platform_id"):
                platforms = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                characteristics.append(line.rstrip("\n")[:2000])
            elif line.startswith("!series_matrix_table_begin"):
                break
    return {
        "path": str(path),
        "title": title,
        "n_samples": n_samples,
        "platforms": sorted(set(platforms)),
        "n_characteristics_rows": len(characteristics),
        "characteristics_examples": characteristics[:10],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gse", nargs="+", help="GEO series accessions, e.g. GSE16879")
    parser.add_argument("--download", action="store_true", help="Actually download files")
    args = parser.parse_args()

    reports = []
    for gse in args.gse:
        gse = gse.upper()
        url = candidate_url(gse)
        out = OUT_DIR / f"{gse}_series_matrix.txt.gz"
        if args.download and not out.exists():
            print(f"Downloading {gse}: {url}")
            download(url, out)
        report = {"gse": gse, "url": url, "local_path": str(out), "exists": out.exists()}
        if out.exists():
            report.update(audit_matrix(out))
        reports.append(report)

    report_path = OUT_DIR / "bulk_geo_series_matrix_audit.json"
    report_path.write_text(json.dumps(reports, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(reports, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
