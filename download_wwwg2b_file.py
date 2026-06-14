#!/usr/bin/env python3
"""Download a WWWG2B OneDrive-backed file with token refresh and validation."""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional, Sequence
from urllib import request

if sys.platform == "win32":
    try:
        import io

        if not isinstance(sys.stdout, io.TextIOWrapper) or sys.stdout.encoding != "utf-8":
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
        if not isinstance(sys.stderr, io.TextIOWrapper) or sys.stderr.encoding != "utf-8":
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
    except (ValueError, AttributeError):
        pass


API_TEMPLATE = (
    "https://wwwg2b.com/api/dataAvailable/table/get_download_url_form_onedrive"
    "?type=availableTable&fileId={file_id}"
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download a WWWG2B file by refreshing the SharePoint temporary URL first.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--file-id", required=True, help="WWWG2B fileId from the dataAvailable fileTable API")
    parser.add_argument("--output", required=True, help="Destination file path")
    parser.add_argument("--expected-size", type=int, default=None, help="Expected complete size in bytes")
    parser.add_argument("--retries", type=int, default=20, help="Number of token refresh/download attempts")
    parser.add_argument("--retry-delay", type=int, default=30, help="Delay between attempts in seconds")
    parser.add_argument("--connect-timeout", type=int, default=60, help="curl connect timeout in seconds")
    parser.add_argument("--speed-time", type=int, default=600, help="curl low-speed timeout window in seconds")
    parser.add_argument("--speed-limit", type=int, default=1, help="curl low-speed limit in bytes/sec")
    parser.add_argument("--curl", default="curl.exe" if sys.platform == "win32" else "curl", help="curl executable")
    return parser


def fetch_download_url(file_id: str) -> str:
    req = request.Request(
        API_TEMPLATE.format(file_id=file_id),
        headers={
            "User-Agent": "Mozilla/5.0",
            "Accept": "application/json,text/plain,*/*",
            "Referer": "https://wwwg2b.com/dataAvailable",
        },
    )
    with request.urlopen(req, timeout=120) as response:
        payload = json.loads(response.read().decode("utf-8"))
    if payload.get("code") != 0 or not payload.get("data"):
        raise RuntimeError(f"WWWG2B API returned: {payload}")
    return str(payload["data"])


def file_starts_with_gzip_magic(path: Path) -> bool:
    if not path.exists() or path.stat().st_size < 2:
        return False
    with path.open("rb") as f:
        return f.read(2) == b"\x1f\x8b"


def looks_like_html(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    with path.open("rb") as f:
        prefix = f.read(256).lstrip().lower()
    return prefix.startswith(b"<html") or prefix.startswith(b"<!doctype html")


def remove_invalid_partial(path: Path) -> None:
    if looks_like_html(path) or (path.exists() and path.stat().st_size > 0 and not file_starts_with_gzip_magic(path)):
        bad_path = path.with_suffix(path.suffix + ".invalid")
        if bad_path.exists():
            bad_path.unlink()
        path.rename(bad_path)
        print(f"[WARN] Existing output was not gzip; moved to {bad_path}")


def is_complete(path: Path, expected_size: Optional[int]) -> bool:
    if not path.exists():
        return False
    if expected_size is not None and path.stat().st_size != expected_size:
        return False
    return file_starts_with_gzip_magic(path)


def run_curl(args: argparse.Namespace, url: str, output: Path) -> int:
    cmd = [
        args.curl,
        "-L",
        "-C",
        "-",
        "--fail",
        "--retry",
        "5",
        "--retry-delay",
        str(args.retry_delay),
        "--retry-all-errors",
        "--connect-timeout",
        str(args.connect_timeout),
        "--speed-time",
        str(args.speed_time),
        "--speed-limit",
        str(args.speed_limit),
        "-o",
        str(output),
        url,
    ]
    return subprocess.run(cmd).returncode


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    remove_invalid_partial(output)

    for attempt in range(1, args.retries + 1):
        if is_complete(output, args.expected_size):
            print(f"[INFO] Complete file already present: {output} ({output.stat().st_size} bytes)")
            return 0

        before = output.stat().st_size if output.exists() else 0
        print(f"[INFO] Attempt {attempt}/{args.retries}; current size={before} bytes")
        try:
            url = fetch_download_url(args.file_id)
        except Exception as exc:
            print(f"[WARN] Could not refresh WWWG2B URL: {exc}")
            time.sleep(args.retry_delay)
            continue

        rc = run_curl(args, url, output)
        after = output.stat().st_size if output.exists() else 0
        print(f"[INFO] curl exit={rc}; size {before} -> {after}")

        if looks_like_html(output):
            remove_invalid_partial(output)
            time.sleep(args.retry_delay)
            continue
        if is_complete(output, args.expected_size):
            print(f"[INFO] Download complete: {output} ({after} bytes)")
            return 0
        if after > 0 and not file_starts_with_gzip_magic(output):
            remove_invalid_partial(output)
            time.sleep(args.retry_delay)
            continue
        if rc != 0:
            time.sleep(args.retry_delay)
            continue

    print(f"[ERROR] Download incomplete after {args.retries} attempts: {output}", file=sys.stderr)
    if output.exists():
        print(f"[ERROR] Current size: {output.stat().st_size} bytes", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
