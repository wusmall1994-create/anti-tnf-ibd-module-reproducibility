"""Verify that a public repository/archive URL is syntactically valid and reachable.

This is a post-upload gate for the manuscript workflow. It does not modify
submission files. It writes a small Markdown/JSON report under `results/` so
the author team can document that the public code/data link was checked before
final submission.

Example:
    python scripts/verify_public_repository_url.py --url https://doi.org/10.xxxx/zenodo.xxxxx
"""

from __future__ import annotations

import argparse
import json
import re
import ssl
import sys
import urllib.error
import urllib.request
from datetime import datetime
from pathlib import Path
from urllib.parse import urlparse


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
REPORT_MD = RESULTS / "public_repository_url_verification.md"
REPORT_JSON = RESULTS / "public_repository_url_verification.json"

KNOWN_PUBLIC_HOST_HINTS = {
    "doi.org": "DOI resolver",
    "zenodo.org": "Zenodo",
    "github.com": "GitHub",
    "osf.io": "OSF",
    "figshare.com": "Figshare",
    "dataverse.harvard.edu": "Harvard Dataverse",
}


def validate_url(url: str) -> tuple[bool, list[str], dict[str, str]]:
    issues: list[str] = []
    parsed = urlparse(url)
    if parsed.scheme not in {"http", "https"}:
        issues.append("URL must start with http:// or https://.")
    if any(ch.isspace() for ch in url):
        issues.append("URL must not contain whitespace.")
    if not parsed.netloc:
        issues.append("URL must include a host name.")
    if parsed.scheme == "http":
        issues.append("Prefer https:// for a public repository/archive URL.")

    host = parsed.netloc.lower()
    host_hint = "recognized public repository/archive host" if host in KNOWN_PUBLIC_HOST_HINTS else "unrecognized host; manual review recommended"
    return not issues, issues, {"scheme": parsed.scheme, "host": host, "path": parsed.path, "host_hint": host_hint}


def fetch_url(url: str, timeout: float = 10.0) -> dict:
    headers = {
        "User-Agent": "ibd-anti-tnf-public-url-verifier/1.0",
        "Accept": "text/html,application/json,*/*;q=0.8",
    }
    context = ssl.create_default_context()
    attempts = []
    for method in ["HEAD", "GET"]:
        req = urllib.request.Request(url, headers=headers, method=method)
        if method == "GET":
            req.add_header("Range", "bytes=0-2047")
        try:
            with urllib.request.urlopen(req, timeout=timeout, context=context) as response:
                final_url = response.geturl()
                status = getattr(response, "status", None) or response.getcode()
                content_type = response.headers.get("Content-Type", "")
                return {
                    "network_status": "pass" if 200 <= int(status) < 400 else "fail",
                    "method": method,
                    "http_status": int(status),
                    "final_url": final_url,
                    "content_type": content_type,
                    "attempts": attempts,
                }
        except urllib.error.HTTPError as exc:
            attempts.append({"method": method, "error_type": "HTTPError", "status": exc.code, "reason": str(exc.reason)})
            if method == "HEAD" and exc.code in {403, 405}:
                continue
            return {
                "network_status": "fail",
                "method": method,
                "http_status": exc.code,
                "final_url": url,
                "content_type": "",
                "attempts": attempts,
                "error": str(exc),
            }
        except Exception as exc:  # noqa: BLE001 - report the exact network failure.
            attempts.append({"method": method, "error_type": type(exc).__name__, "error": str(exc)})
            if method == "HEAD":
                continue
            return {
                "network_status": "fail",
                "method": method,
                "http_status": None,
                "final_url": url,
                "content_type": "",
                "attempts": attempts,
                "error": str(exc),
            }
    return {"network_status": "fail", "attempts": attempts, "final_url": url}


def write_report(data: dict) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    REPORT_JSON.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# Public repository/archive URL verification",
        "",
        f"Generated: {data['generated']}",
        "",
        f"Overall status: **{data['overall_status'].upper()}**",
        "",
        "| Field | Value |",
        "|---|---|",
        f"| URL | `{data['url']}` |",
        f"| Syntax status | {data['syntax_status']} |",
        f"| Host | `{data['parsed'].get('host', '')}` |",
        f"| Host hint | {data['parsed'].get('host_hint', '')} |",
        f"| Network status | {data['network'].get('network_status', 'not_run')} |",
        f"| HTTP status | {data['network'].get('http_status', '')} |",
        f"| Final URL after redirects | `{data['network'].get('final_url', '')}` |",
        f"| Content type | `{data['network'].get('content_type', '')}` |",
        "",
    ]
    if data["issues"]:
        lines.extend(["## Issues", ""])
        lines.extend(f"- {issue}" for issue in data["issues"])
        lines.append("")
    if data["network"].get("attempts"):
        lines.extend(["## Network attempts", ""])
        for attempt in data["network"]["attempts"]:
            lines.append(f"- `{attempt}`")
        lines.append("")
    lines.extend(
        [
            "## Interpretation",
            "",
            "- `PASS` means the URL is syntactically valid and returned an HTTP 2xx/3xx response during this check.",
            "- `WARN` means the URL is syntactically acceptable but the host is not one of the common public archive hosts; manual review is recommended.",
            "- `FAIL` means the URL should not be inserted into final submission files until fixed or manually verified.",
        ]
    )
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify a public repository/archive URL.")
    parser.add_argument("--url", required=True, help="Public repository/archive URL or DOI URL")
    parser.add_argument("--skip-network", action="store_true", help="Only perform syntax/host checks")
    parser.add_argument("--timeout", type=float, default=10.0, help="Per-request network timeout in seconds")
    args = parser.parse_args()

    syntax_ok, issues, parsed = validate_url(args.url)
    network = {"network_status": "skipped", "final_url": args.url} if args.skip_network else fetch_url(args.url, timeout=args.timeout)

    host_known = parsed.get("host") in KNOWN_PUBLIC_HOST_HINTS
    if not syntax_ok:
        overall = "fail"
    elif network.get("network_status") == "fail":
        overall = "fail"
    elif not host_known:
        overall = "warn"
    else:
        overall = "pass"

    data = {
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "url": args.url,
        "syntax_status": "pass" if syntax_ok else "fail",
        "issues": issues,
        "parsed": parsed,
        "network": network,
        "overall_status": overall,
    }
    write_report(data)
    print(json.dumps({"overall_status": overall, "network_status": network.get("network_status"), "report": str(REPORT_MD.relative_to(ROOT))}, ensure_ascii=False))
    if overall == "fail":
        sys.exit(1)


if __name__ == "__main__":
    main()
