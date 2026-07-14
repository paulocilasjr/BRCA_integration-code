# Reproducibility record

This record describes the exact primary workflow audited on 2026-07-13. It is
intended to accompany the code/data-availability statement and to make reviewer
checks possible without reconstructing the analysis history.

## Validated environment

- macOS arm64
- Python 3.14.6
- Direct dependencies pinned in `requirements.txt`; the fully resolved runtime
  is pinned in `requirements-lock.txt`
- Headless Matplotlib `Agg` backend
- Source-workbook and EVE archive integrity enforced with SHA-256

The project declares Python 3.11 or newer. For long-term archival reproduction,
record the OS, Python version, Git commit, and `python -m pip freeze` alongside
the generated artifacts.

## Active inputs

| Input | Role | Integrity/provenance |
| --- | --- | --- |
| `dataset/SUPP_TABLES_BRCA12_APR_2026.xlsx` | Master variant, metadata, reference-panel, and comments sheets (Tables 1–6) | Deposited; hash in `checksums.sha256` |
| `dataset/ACMG_other_points.xlsx` | PM2, segregation, and allele-frequency point inputs for Tables 18–19 | Deposited; hash in `checksums.sha256` |
| `BRCA1_HUMAN.EVE.variants.zip` | BRCA1 EVE predictions | EVE full-length model; downloaded and checksum-verified |
| `BRCA2_HUMAN.EVE.variants.zip` | BRCA2 EVE predictions | EVE full-length model; downloaded and checksum-verified |
| `src/brca_integration/reference/domains.py` | Curated BRCA1/BRCA2 protein feature intervals used in Table 16 | Version-controlled constants |

The EVE normalizer joins on wild-type amino acid (`T2`), one-based protein
position (`T3`), and alternate amino acid (`T4`). It retains the source variant
fields (`T1`–`T7`) and records the archive/member name. Expected full-master
coverage is 10,958/11,009 for BRCA1 and 15,542/20,169 for BRCA2.

## Deposited legacy inputs

`dataset/AlphaMissense_Calculations_all.xlsx` and
`dataset/evidence_criteria_v6_BRCA12.xlsx` are retained as analysis-history
artifacts. They are not read by the current primary pipeline. Tables 18–19 use
EVE, not AlphaMissense, to avoid predictor circularity described in
`docs/eve_replacement_implementation.md`.

## Explicit curated overrides

All overrides are visible near the top of
`src/brca_integration/tables/sup_table_18_19.py` rather than being silently
edited in an output workbook:

| Gene | Variant | Override |
| --- | --- | --- |
| BRCA1 | p.R71G | segregation points = 8 |
| BRCA2 | p.D2312V | segregation points = -4 |
| BRCA2 | p.R2336H | segregation points = 4; final class = P |
| BRCA2 | p.Q2829R | segregation points = 2 |

The same module also contains named note overrides for variants with known or
potential splicing effects/co-occurrence. These notes do not alter calculated
points except where the table above explicitly says so.

## Determinism and expected outputs

Stratified cross-validation uses shuffled folds with a fixed random seed of 42;
see `docs/cross_validation_stratified.md`. Workbook ZIP metadata and PDF/SVG
serialization can vary across library/platform versions, so byte-identical
output hashes are not used as the scientific equivalence criterion. The
verification script instead checks invariant sheet structure, row counts, EVE
coverage, and final classifications.

The complete workbook has these sheets, in order:

```text
Sup Table 7, Sup Table 8, Sup Table 9, Sup Table 10, Sup Table 11,
Sup Table 12, Sup Table 13, Supp Table 14, Supp Table 15, SuppTable 16,
Supp Table 17, Sup Table 18, Sup Table 19
```

The naming differences (`Sup`, `Supp`, and `SuppTable`) are retained to match
the manuscript workbook and downstream references.

## Audit commands

```sh
python -m compileall -q main.py scripts src
python -m pytest
python scripts/build_eve_artifacts.py
python main.py \
  --output-workbook results/SUPP_TABLES_BRCA12_reproduced.xlsx \
  --figure-prefix figures/reproduced/supp_fig2
PYTHONPATH=src python scripts/verify_reproduction.py \
  results/SUPP_TABLES_BRCA12_reproduced.xlsx \
  --figure-prefix figures/reproduced/supp_fig2
```

Before journal deposition, create an immutable tagged release and archive it
with a DOI. Update `CITATION.cff` with the final article DOI and full manuscript
author list when those bibliographic details are assigned.
