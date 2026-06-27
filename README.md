# BRCA Integration Code

This repository contains the reproducible code used to generate the BRCA1/BRCA2
supplementary tables and supporting figures.

## Repository Layout

```text
dataset/                     Source workbooks and reproducibility notes
scripts/                     Command-line helper scripts
src/brca_integration/        Python package
  pipeline.py                End-to-end supplementary workbook pipeline
  cli.py                     Command-line interface
  tables/                    Supplementary table writers
  figures/                   Figure generation
  reference/                 Curated reference constants
  analyses/                  Secondary analysis scripts
results/                     Generated workbooks (ignored by Git)
figures/                     Generated figure files (ignored by Git)
```

## Environment

Create or activate a Python environment with the packages in `requirements.txt`.
For example:

```sh
python3 -m pip install -r requirements.txt
```

The repository has also been run successfully from Conda base with Python 3.12
when the listed packages are installed.

## Rebuild The Supplementary Tables

From the repository root:

```sh
python3 main.py
```

Equivalent explicit entry point:

```sh
python3 scripts/build_supplementary_tables.py
```

Package module entry point:

```sh
PYTHONPATH=src python3 -m brca_integration
```

If installed as a package, the console script is:

```sh
brca-build-tables
```

By default the output workbook is timestamped:

```text
results/SUPP_TABLES_BRCA12_<YYYYMMDD_HHMMSS>.xlsx
```

Supp Fig 2 files are written under:

```text
figures/<workbook_stem>/supp_fig2.{png,pdf,svg}
```

Useful options:

```sh
python3 main.py --output-dir results
python3 main.py --timestamp 20260626_120000
python3 main.py --output-workbook results/SUPP_TABLES_BRCA12_review.xlsx
```

## EVE Predictor Artifacts

The large EVE archives and generated EVE workbook are intentionally ignored by
Git. Recreate them before running the full integration if
`dataset/eve/EVE_BRCA12_scores.xlsx` is missing:

```sh
python3 scripts/build_eve_artifacts.py
```

See `dataset/eve/README.md` for source URLs, checksum details, and expected
coverage.

## Current Primary Output

The end-to-end pipeline writes Sup Tables 7-13, Supp Tables 14-17, and Sup
Tables 18-19 into a single workbook.
