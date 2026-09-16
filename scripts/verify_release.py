#!/usr/bin/env python3
"""Verify the integrity and basic hygiene of the public release archive."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "checksums" / "SHA256SUMS.tsv"
UNWANTED_NAMES = {".DS_Store", "__pycache__"}
UNWANTED_SUFFIXES = {".pyc", ".pyo", ".swp", ".tmp", ".bak"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    failures: list[str] = []
    for row in rows:
        path = ROOT / row["file"]
        if not path.is_file():
            failures.append(f"missing: {row['file']}")
            continue
        if str(path.stat().st_size) != row["size_bytes"]:
            failures.append(f"size mismatch: {row['file']}")
        if sha256(path) != row["sha256"]:
            failures.append(f"checksum mismatch: {row['file']}")

    junk = [
        str(path.relative_to(ROOT))
        for path in ROOT.rglob("*")
        if path.name in UNWANTED_NAMES or path.suffix.lower() in UNWANTED_SUFFIXES
    ]
    failures.extend(f"unwanted file: {item}" for item in junk)

    if failures:
        for item in failures:
            print(item)
        raise SystemExit("RELEASE_INTEGRITY=FAIL")

    print(f"FILES_VERIFIED={len(rows)}")
    print("UNWANTED_FILE_COUNT=0")
    print("RELEASE_INTEGRITY=PASS")


if __name__ == "__main__":
    main()

