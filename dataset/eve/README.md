# EVE Reproducibility Artifacts

This directory contains lightweight instructions and code for recreating the
large EVE files used by the BRCA1/BRCA2 integration analysis. The downloaded
archives and generated workbook are intentionally not committed to GitHub.

## Inputs

The script expects the repository input workbook:

```text
dataset/SUPP_TABLES_BRCA12_APR_2026.xlsx
```

It downloads the full-length EVE models:

```text
https://evemodel.org/api/proteins/web_pid/BRCA1_HUMAN/download/?variants=True
https://evemodel.org/api/proteins/web_pid/BRCA2_HUMAN/download/?variants=True
```

Model identifiers:

```text
BRCA1_HUMAN, UniProt P38398
BRCA2_HUMAN, UniProt P51587
```

The BRCA1 `BRCA1_BRCT_HUMAN` and `BRCA1_RING_HUMAN` domain models are not used
because the analyzed BRCA1 variant tables cover the full coding sequence.

## Generate The EVE Workbook

From the repository root, run:

```sh
python3 scripts/build_eve_artifacts.py
```

If the EVE HTTPS certificate chain is rejected by your local Python
installation, rerun after confirming the URLs above:

```sh
python3 scripts/build_eve_artifacts.py --insecure-tls
```

The script creates:

```text
dataset/eve/BRCA1_HUMAN.EVE.variants.zip
dataset/eve/BRCA2_HUMAN.EVE.variants.zip
dataset/eve/EVE_BRCA12_scores.xlsx
```

By default it verifies that the downloaded EVE archives match the SHA256 hashes
used for the reported analysis:

```text
BRCA1_HUMAN.EVE.variants.zip  d99a1f4b383154afdec9cca35e5a27c91184f63e41a01805461a9b8d280a0b39
BRCA2_HUMAN.EVE.variants.zip  ba15c69673dab8bb6ce73f18407e10bf4106eebdafdf84caaea9db39a01038d8
```

If EVE republishes a file and you intentionally want to regenerate with the
newer source, use:

```sh
python3 scripts/build_eve_artifacts.py --force-download --allow-checksum-mismatch
```

## Generate The Integrated Tables

After `dataset/eve/EVE_BRCA12_scores.xlsx` exists, run the repository pipeline:

```sh
python3 main.py
```

This writes a timestamped generated supplementary workbook to:

```text
results/SUPP_TABLES_BRCA12_<timestamp>.xlsx
```

Sup Tables 18 and 19 use the EVE workbook through
`brca_integration.tables.sup_table_18_19`.

## Expected EVE Coverage

The generated EVE workbook should report:

```text
BRCA1: 10,958/11,009 scored; 51 missing; 99.54%
BRCA2: 15,542/20,169 scored; 4,627 missing; 77.06%
```

In the final integrated tables, the expected EVE score coverage is:

```text
Sup Table 18: 3,218/3,247 scored; 29 missing; 99.11%
Sup Table 19: 5,601/6,177 scored; 576 missing; 90.68%
```
