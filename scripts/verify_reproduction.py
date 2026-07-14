#!/usr/bin/env python3
"""Verify deposited inputs and the primary publication artifacts."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = REPO_ROOT / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from brca_integration.validation import validate_generated_workbook  # noqa: E402


EXPECTED_INPUT_SHA256 = {
    "dataset/SUPP_TABLES_BRCA12_APR_2026.xlsx": (
        "71e7d32f0ec3ab6f69633b6ae90167a7acd2cdf24567751f61e41bcf71f42904"
    ),
    "dataset/ACMG_other_points.xlsx": (
        "12eb9270ac578bb621a64716a6ab3555f5acc25fd53b40a419977306f232dcd9"
    ),
}
EXPECTED_TABLE_METRICS = {
    "Sup Table 18": {
        "rows": 3_247,
        "eve_scored": 3_218,
        "classes": {"B": 83, "LB": 2_205, "VUS": 585, "LP": 364, "P": 10},
    },
    "Sup Table 19": {
        "rows": 6_177,
        "eve_scored": 5_601,
        "classes": {"B": 104, "LB": 3_881, "VUS": 1_873, "LP": 314, "P": 5},
    },
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify checksums, workbook structure, expected counts, and figures."
    )
    parser.add_argument("workbook", type=Path, help="Generated supplementary workbook.")
    parser.add_argument(
        "--figure-prefix",
        type=Path,
        help="Optional Supp Fig 2 prefix; checks .png, .pdf, and .svg files.",
    )
    parser.add_argument(
        "--skip-input-checksums",
        action="store_true",
        help="Skip checksums when intentionally using replacement source workbooks.",
    )
    return parser.parse_args()


def verify_input_checksums() -> None:
    for relative_path, expected in EXPECTED_INPUT_SHA256.items():
        path = REPO_ROOT / relative_path
        if not path.is_file():
            raise SystemExit(f"[FAIL] Missing deposited input: {path}")
        observed = sha256_file(path)
        if observed != expected:
            raise SystemExit(
                f"[FAIL] SHA256 mismatch for {relative_path}\n"
                f"  expected: {expected}\n  observed: {observed}"
            )
        print(f"[OK] {relative_path} SHA256")


def verify_metrics(workbook: Path) -> None:
    for sheet_name, expected in EXPECTED_TABLE_METRICS.items():
        frame = pd.read_excel(workbook, sheet_name=sheet_name, header=1)
        observed_rows = len(frame)
        observed_scored = int(frame["EVE Score"].notna().sum())
        observed_classes = {
            label: int(count)
            for label, count in frame["FINAL ClinVar Classification"].value_counts().items()
        }
        if observed_rows != expected["rows"]:
            raise SystemExit(
                f"[FAIL] {sheet_name} rows: {observed_rows}; expected {expected['rows']}"
            )
        if observed_scored != expected["eve_scored"]:
            raise SystemExit(
                f"[FAIL] {sheet_name} EVE scores: {observed_scored}; "
                f"expected {expected['eve_scored']}"
            )
        if observed_classes != expected["classes"]:
            raise SystemExit(
                f"[FAIL] {sheet_name} final classes: {observed_classes}; "
                f"expected {expected['classes']}"
            )
        print(
            f"[OK] {sheet_name}: {observed_rows:,} rows; "
            f"{observed_scored:,} EVE scores; final classes {observed_classes}"
        )


def verify_figures(prefix: Path) -> None:
    for suffix in (".png", ".pdf", ".svg"):
        path = prefix.with_suffix(suffix)
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"[FAIL] Missing or empty figure: {path}")
        print(f"[OK] {path} ({path.stat().st_size:,} bytes)")


def main() -> None:
    args = parse_args()
    if not args.skip_input_checksums:
        verify_input_checksums()
    validate_generated_workbook(args.workbook)
    print(f"[OK] Workbook structure and error scan: {args.workbook}")
    verify_metrics(args.workbook)
    if args.figure_prefix:
        verify_figures(args.figure_prefix)
    print("[OK] Reproduction verification passed")


if __name__ == "__main__":
    main()
