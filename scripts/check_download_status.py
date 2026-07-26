from __future__ import annotations

from pathlib import Path


DATA = Path("D:/TAURUS_data")
TARGET = DATA / "GSE282122_filtered_processed_data.tar.gz"
EXPECTED_BYTES = 3_027_066_520  # Content-Length from NCBI FTP HTTP header.
LOG = DATA / "download_filtered_processed.err.log"


def main() -> None:
    if not TARGET.exists():
        print(f"NOT_STARTED\t{TARGET}")
        return
    size = TARGET.stat().st_size
    pct = size / EXPECTED_BYTES * 100
    print(f"FILE\t{TARGET}")
    print(f"BYTES\t{size}")
    print(f"APPROX_PERCENT\t{pct:.2f}")
    if LOG.exists():
        lines = LOG.read_text(errors="ignore").splitlines()
        tail = [x for x in lines[-12:] if x.strip()]
        print("LOG_TAIL")
        print("\n".join(tail))


if __name__ == "__main__":
    main()
