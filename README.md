# BRCA1/BRCA2 functional evidence integration

[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

This repository contains the code and deposited inputs used to generate the
BRCA1/BRCA2 supplementary tables and Supplementary Figure 2 for the associated
human-genetics manuscript. The pipeline integrates multiplexed functional assay
tracks, reference-variant evidence, protein-domain annotations, EVE predictions,
and other ACMG/AMP evidence into one reviewable workbook.

The primary reproducible products are:

- Supplementary Tables 7–19 in a single `.xlsx` workbook (13 worksheets).
- Supplementary Figure 2 as PNG, PDF, and SVG.

## Reproduce the results

The commands below start from a fresh clone. Python 3.14.6 was used for the
final repository audit; Python 3.11 or newer is required.

```sh
git clone https://github.com/paulocilasjr/BRCA_integration-code.git
cd BRCA_integration-code
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-lock.txt
```

Download the two EVE archives and build the normalized predictor workbook:

```sh
python scripts/build_eve_artifacts.py
```

The script verifies the exact EVE archives used in the analysis by SHA-256. If
the EVE server's certificate chain is rejected locally, review the documented
URLs in [`dataset/eve/README.md`](dataset/eve/README.md) and use
`--insecure-tls`; checksums remain mandatory.

Build and verify all publication artifacts:

```sh
python main.py \
  --output-workbook results/SUPP_TABLES_BRCA12_reproduced.xlsx \
  --figure-prefix figures/reproduced/supp_fig2

PYTHONPATH=src python scripts/verify_reproduction.py \
  results/SUPP_TABLES_BRCA12_reproduced.xlsx \
  --figure-prefix figures/reproduced/supp_fig2
```

Or run the equivalent make targets:

```sh
make reproduce
make verify
```

Expect the full run to take about one minute after the dependencies and EVE
archives are available. It requires approximately 100 MB of temporary/download
space. The pipeline writes to temporary paths, validates the workbook, and only
then publishes the requested output, so a failed run does not leave a partial
workbook at the final path.

## Expected verification results

The verifier checks input hashes, required sheets, spreadsheet errors, EVE
coverage, figure presence, and final classification counts.

| Output | Rows | EVE scored | B | LB | VUS | LP | P |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Sup Table 18 (BRCA1) | 3,247 | 3,218 | 83 | 2,205 | 585 | 364 | 10 |
| Sup Table 19 (BRCA2) | 6,177 | 5,601 | 104 | 3,881 | 1,873 | 314 | 5 |

## Repository layout

```text
dataset/                      Deposited source workbooks
  eve/README.md               EVE sources, hashes, and coverage
docs/                         Methods and reproducibility notes
scripts/
  build_eve_artifacts.py      Download/normalize checksum-pinned EVE data
  build_supplementary_tables.py
  verify_reproduction.py      Machine-check the publication outputs
src/brca_integration/
  pipeline.py                 Atomic end-to-end pipeline
  validation.py               Input/output integrity checks
  tables/                     Supplementary table calculations and writers
  figures/                    Supplementary Figure 2 generation
  reference/                  Curated BRCA1/BRCA2 domain constants
tests/                        Fast regression tests
checksums.sha256              SHA-256 manifest for deposited workbooks
requirements-lock.txt         Fully resolved audited runtime
```

Generated `results/`, `figures/`, and large downloaded EVE artifacts are ignored
by Git. The original source workbook contains Supplementary Tables 1–6; this
pipeline regenerates Tables 7–19.

## Inputs, provenance, and manual curation

[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) identifies every active
input, records expected versions and counts, distinguishes legacy workbooks from
active inputs, and lists the small number of explicit manual overrides. The
deposited workbook hashes can also be checked independently:

```sh
shasum -a 256 -c checksums.sha256
```

The EVE files are not redistributed here because they are large and externally
hosted. Their source URLs, model identifiers, exact archive hashes, join keys,
and expected coverage are documented in
[`dataset/eve/README.md`](dataset/eve/README.md).

## Other entry points

All of these invoke the same primary pipeline:

```sh
python main.py
python scripts/build_supplementary_tables.py
PYTHONPATH=src python -m brca_integration
python -m pip install -e .
brca-build-tables
```

Run `python main.py --help` for input/output overrides. Supplying replacement
inputs is supported, but the publication verifier's expected hashes and counts
then no longer apply; use `--skip-input-checksums` only for an intentional
sensitivity or update analysis.

The modules under `src/brca_integration/analyses/` are retained exploratory
manuscript-support analyses and are not part of the primary reproduction path.
Some require historical source files that are not deposited; they must not be
used to regenerate the reported tables.

## Testing

```sh
python -m pip install -e '.[test]'
pytest
```

The full publication check is the end-to-end build followed by
`scripts/verify_reproduction.py`, not the fast unit suite alone.

## Citation and license

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Once the
associated article receives its final DOI, cite the article as the scientific
source and archive the corresponding repository release (for example, on
Zenodo) so the manuscript points to an immutable version.

The software is released under the [MIT License](LICENSE). Source workbooks and
third-party EVE data may be subject to separate terms; the software license does
not relicense those data.

## Scope and interpretation

This repository is research software. The generated classifications support the
reported analysis and are not independently validated for clinical diagnostic
use.
