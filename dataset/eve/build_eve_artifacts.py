#!/usr/bin/env python3
"""
Download EVE BRCA1/BRCA2 variant predictions and build the normalized workbook.

This script recreates the ignored large files used by the repository:

  dataset/eve/BRCA1_HUMAN.EVE.variants.zip
  dataset/eve/BRCA2_HUMAN.EVE.variants.zip
  dataset/eve/EVE_BRCA12_scores.xlsx

Run from the repository root:

  python3 dataset/eve/build_eve_artifacts.py

If your local Python TLS stack rejects the EVE certificate chain, rerun with
--insecure-tls after confirming the HTTPS endpoints in the README.
"""

from __future__ import annotations

import argparse
import hashlib
import ssl
import sys
import zipfile
from pathlib import Path
from urllib.error import URLError
from urllib.request import Request, urlopen

try:
    import pandas as pd
except ImportError as exc:  # pragma: no cover - only used for user-facing setup errors
    raise SystemExit(
        "Missing dependency: pandas. Install pandas and openpyxl, then rerun this script."
    ) from exc


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EVE_DIR = Path(__file__).resolve().parent
DEFAULT_SOURCE_WORKBOOK = REPO_ROOT / "dataset" / "SUPP_TABLES_BRCA12_APR_2026.xlsx"
DEFAULT_OUTPUT_WORKBOOK = DEFAULT_EVE_DIR / "EVE_BRCA12_scores.xlsx"

EVE_SOURCE_SITE = "https://evemodel.org/"
EVE_DOWNLOAD_TEMPLATE = "https://evemodel.org/api/proteins/web_pid/{web_pid}/download/?variants=True"

EVE_COLUMNS = [
    "wt_aa",
    "position",
    "mt_aa",
    "EVE_scores_ASM",
    "EVE_classes_75_pct_retained_ASM",
]
MASTER_COLUMNS = ["T1", "T2", "T3", "T4", "T5", "T6", "T7"]

PROTEINS = {
    "BRCA1": {
        "web_pid": "BRCA1_HUMAN",
        "uniprot": "P38398",
        "master_sheet": "Sup Table 1",
        "zip_name": "BRCA1_HUMAN.EVE.variants.zip",
        "member": "BRCA1_HUMAN.csv",
        "sha256": "d99a1f4b383154afdec9cca35e5a27c91184f63e41a01805461a9b8d280a0b39",
    },
    "BRCA2": {
        "web_pid": "BRCA2_HUMAN",
        "uniprot": "P51587",
        "master_sheet": "Sup Table 2",
        "zip_name": "BRCA2_HUMAN.EVE.variants.zip",
        "member": "BRCA2_HUMAN.csv",
        "sha256": "ba15c69673dab8bb6ce73f18407e10bf4106eebdafdf84caaea9db39a01038d8",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Recreate BRCA1/BRCA2 EVE artifacts used by Sup Tables 18/19."
    )
    parser.add_argument(
        "--source-workbook",
        default=str(DEFAULT_SOURCE_WORKBOOK),
        help="Workbook containing Sup Table 1 and Sup Table 2.",
    )
    parser.add_argument(
        "--eve-dir",
        default=str(DEFAULT_EVE_DIR),
        help="Directory for downloaded EVE archives.",
    )
    parser.add_argument(
        "--output-workbook",
        default=str(DEFAULT_OUTPUT_WORKBOOK),
        help="Normalized EVE workbook to write.",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Do not download missing EVE archives; require them to already exist.",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Download EVE archives even when local copies already exist.",
    )
    parser.add_argument(
        "--allow-checksum-mismatch",
        action="store_true",
        help="Warn instead of failing if a downloaded/local EVE archive differs from the recorded SHA256.",
    )
    parser.add_argument(
        "--insecure-tls",
        action="store_true",
        help="Disable TLS certificate verification for EVE downloads.",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_sha256(path: Path, expected: str, allow_mismatch: bool) -> None:
    observed = sha256_file(path)
    if observed == expected:
        print(f"[OK] {path.name} SHA256 verified")
        return
    message = (
        f"{path} SHA256 mismatch.\n"
        f"  expected: {expected}\n"
        f"  observed: {observed}\n"
        "The EVE download may have changed. Use --allow-checksum-mismatch only if you "
        "intend to regenerate results from the newer source file."
    )
    if allow_mismatch:
        print(f"[WARN] {message}", file=sys.stderr)
        return
    raise SystemExit(message)


def download_file(url: str, path: Path, insecure_tls: bool) -> None:
    print(f"[INFO] Downloading {url}")
    request = Request(url, headers={"User-Agent": "BRCA-integration-code/1.0"})
    context = ssl._create_unverified_context() if insecure_tls else None
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    try:
        with urlopen(request, context=context) as response, tmp_path.open("wb") as handle:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
    except URLError as exc:
        raise SystemExit(
            f"Failed to download {url}: {exc}\n"
            "If this is a local certificate-chain issue, retry with --insecure-tls."
        ) from exc
    tmp_path.replace(path)
    print(f"[OK] Wrote {path}")


def ensure_archives(
    eve_dir: Path,
    skip_download: bool,
    force_download: bool,
    insecure_tls: bool,
    allow_checksum_mismatch: bool,
) -> dict[str, Path]:
    archive_paths: dict[str, Path] = {}
    for gene, info in PROTEINS.items():
        path = eve_dir / info["zip_name"]
        url = EVE_DOWNLOAD_TEMPLATE.format(web_pid=info["web_pid"])
        if force_download or not path.exists():
            if skip_download:
                raise SystemExit(f"Missing {path}; rerun without --skip-download to download it.")
            download_file(url, path, insecure_tls)
        else:
            print(f"[OK] Using cached {path}")
        verify_sha256(path, info["sha256"], allow_checksum_mismatch)
        archive_paths[gene] = path
    return archive_paths


def load_eve_archive(zip_path: Path, expected_member: str) -> pd.DataFrame:
    with zipfile.ZipFile(zip_path) as archive:
        members = archive.namelist()
        if expected_member not in members:
            raise SystemExit(
                f"{zip_path} does not contain {expected_member}. Available members: {members}"
            )
        with archive.open(expected_member) as handle:
            eve = pd.read_csv(handle, usecols=EVE_COLUMNS)

    eve = eve.rename(
        columns={
            "wt_aa": "T2",
            "position": "T3",
            "mt_aa": "T4",
            "EVE_scores_ASM": "EVE Score",
            "EVE_classes_75_pct_retained_ASM": "EVE classification",
        }
    )
    eve["T2"] = eve["T2"].astype(str).str.strip()
    eve["T4"] = eve["T4"].astype(str).str.strip()
    eve["T3"] = pd.to_numeric(eve["T3"], errors="coerce")
    if eve.duplicated(["T2", "T3", "T4"]).any():
        duplicates = eve.loc[eve.duplicated(["T2", "T3", "T4"], keep=False), ["T2", "T3", "T4"]]
        raise SystemExit(f"Duplicate EVE variant keys found in {zip_path}:\n{duplicates.head()}")
    return eve


def load_master(source_workbook: Path, sheet_name: str) -> pd.DataFrame:
    master = pd.read_excel(source_workbook, sheet_name=sheet_name, header=1)
    missing = [column for column in MASTER_COLUMNS if column not in master.columns]
    if missing:
        raise SystemExit(f"{source_workbook} sheet {sheet_name} is missing columns: {missing}")
    out = master[MASTER_COLUMNS].copy()
    out["T2"] = out["T2"].astype(str).str.strip()
    out["T4"] = out["T4"].astype(str).str.strip()
    out["T3"] = pd.to_numeric(out["T3"], errors="coerce")
    return out


def build_normalized_workbook(
    source_workbook: Path,
    archive_paths: dict[str, Path],
    output_workbook: Path,
) -> None:
    if not source_workbook.exists():
        raise SystemExit(f"Missing source workbook: {source_workbook}")

    output_workbook.parent.mkdir(parents=True, exist_ok=True)
    coverage_rows = []

    readme = pd.DataFrame(
        [
            ("source_site", EVE_SOURCE_SITE),
            ("download_endpoint_template", EVE_DOWNLOAD_TEMPLATE),
            ("brca1_source_url", EVE_DOWNLOAD_TEMPLATE.format(web_pid=PROTEINS["BRCA1"]["web_pid"])),
            ("brca2_source_url", EVE_DOWNLOAD_TEMPLATE.format(web_pid=PROTEINS["BRCA2"]["web_pid"])),
            ("brca1_model", f'{PROTEINS["BRCA1"]["web_pid"]}; UniProt {PROTEINS["BRCA1"]["uniprot"]}'),
            ("brca2_model", f'{PROTEINS["BRCA2"]["web_pid"]}; UniProt {PROTEINS["BRCA2"]["uniprot"]}'),
            ("score_column", "EVE_scores_ASM"),
            ("class_column", "EVE_classes_75_pct_retained_ASM"),
            ("score_direction", "1 = most pathogenic, 0 = most benign"),
            ("class_mapping", "Benign -> -1; Pathogenic -> +1; Uncertain/no score -> 0"),
        ],
        columns=["field", "value"],
    )

    with pd.ExcelWriter(output_workbook, engine="openpyxl") as writer:
        readme.to_excel(writer, sheet_name="README", index=False)
        for gene, info in PROTEINS.items():
            zip_path = archive_paths[gene]
            eve = load_eve_archive(zip_path, info["member"])
            eve["EVE Source"] = f"{zip_path.name}:{info['member']}"
            master = load_master(source_workbook, info["master_sheet"])
            merged = master.merge(eve, on=["T2", "T3", "T4"], how="left")
            merged.insert(0, "INDEX", range(1, len(merged) + 1))
            merged.to_excel(writer, sheet_name=gene, index=False)
            scored = int(merged["EVE Score"].notna().sum())
            coverage_rows.append((gene, len(merged), scored, len(merged) - scored, scored / len(merged)))

    print(f"[OK] Wrote normalized workbook: {output_workbook}")
    print("[INFO] Coverage after joining EVE to full BRCA master tables:")
    for gene, rows, scored, missing, coverage in coverage_rows:
        print(f"  {gene}: {scored:,}/{rows:,} scored; {missing:,} missing; {coverage:.2%}")


def main() -> None:
    args = parse_args()
    source_workbook = Path(args.source_workbook).resolve()
    eve_dir = Path(args.eve_dir).resolve()
    output_workbook = Path(args.output_workbook).resolve()

    archive_paths = ensure_archives(
        eve_dir=eve_dir,
        skip_download=args.skip_download,
        force_download=args.force_download,
        insecure_tls=args.insecure_tls,
        allow_checksum_mismatch=args.allow_checksum_mismatch,
    )
    build_normalized_workbook(source_workbook, archive_paths, output_workbook)


if __name__ == "__main__":
    main()
