from __future__ import annotations

import argparse
from pathlib import Path

from .config import (
    DEFAULT_EVE_WORKBOOK,
    DEFAULT_INPUT_WORKBOOK,
    DEFAULT_OTHER_POINTS_WORKBOOK,
    DEFAULT_RESULTS_DIR,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the BRCA1/BRCA2 supplementary tables workbook."
    )
    parser.add_argument(
        "--input-workbook",
        default=str(DEFAULT_INPUT_WORKBOOK),
        help="Source workbook containing Sup Tables 1-6.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_RESULTS_DIR),
        help="Directory for timestamped SUPP_TABLES_BRCA12_<timestamp>.xlsx output.",
    )
    parser.add_argument(
        "--output-workbook",
        default=None,
        help="Explicit output workbook path. Overrides --output-dir and --timestamp naming.",
    )
    parser.add_argument(
        "--timestamp",
        default=None,
        help="Timestamp token for the default output filename. Defaults to current local time.",
    )
    parser.add_argument(
        "--eve-workbook",
        default=str(DEFAULT_EVE_WORKBOOK),
        help="Normalized EVE predictor workbook used by Sup Tables 18/19.",
    )
    parser.add_argument(
        "--other-points-workbook",
        default=str(DEFAULT_OTHER_POINTS_WORKBOOK),
        help="Workbook containing additional ACMG point inputs.",
    )
    parser.add_argument(
        "--figure-prefix",
        default=None,
        help="Prefix for generated Supp Fig 2 files. Defaults to figures/<workbook_stem>/supp_fig2.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    from .pipeline import build_supplementary_workbook

    outputs = build_supplementary_workbook(
        input_workbook=Path(args.input_workbook),
        output_workbook=Path(args.output_workbook) if args.output_workbook else None,
        output_dir=Path(args.output_dir),
        timestamp=args.timestamp,
        eve_workbook=Path(args.eve_workbook),
        other_points_workbook=Path(args.other_points_workbook),
        figure_prefix=Path(args.figure_prefix) if args.figure_prefix else None,
    )

    print(f"Wrote supplementary workbook: {outputs.workbook}")
    print(f"Wrote Supp Fig 2 files with prefix: {outputs.figure_prefix}")


if __name__ == "__main__":
    main()
