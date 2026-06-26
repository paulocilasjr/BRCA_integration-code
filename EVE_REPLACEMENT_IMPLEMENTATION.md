# EVE Replacement for AlphaMissense in Sup Tables 18 and 19

This note documents the final in silico predictor replacement used to address the reviewer comment about AlphaMissense circularity. The final method uses EVE rather than AlphaMissense or CPT-1 because EVE is fully unsupervised, population-free, available for both BRCA1 and BRCA2 as the same predictor family, and provides author-defined variant classes.

## Decision

The preferred reviewer response is EVE.

Reasons:

- EVE is described by its authors as a fully unsupervised model trained from evolutionary amino acid sequences.
- EVE scores range from 1, most pathogenic, to 0, most benign.
- The downloaded EVE files include author-provided `EVE_classes_75_pct_retained_ASM` calls.
- Using these author-provided classes avoids fitting new clinical thresholds on our BRCA reference labels.
- Both BRCA1 and BRCA2 have full-protein EVE entries, so the method is consistent across genes.

CPT-1 was evaluated first as a fallback because the reviewer also mentioned CPT-1. It has better integrated-tab score coverage than EVE, but it is not as clean for this specific response because CPT-1 is trained with DMS functional assay data and uses EVE-dependent features when available. Since EVE coverage is usable for both integrated tabs, EVE is the final implementation.

## End-to-End Reproduction

The GitHub repository is:

```text
https://github.com/paulocilasjr/BRCA_integration-code
```

The final EVE archives, normalized EVE workbook, and generated result workbook are not committed because they are large/generated files. They can be recreated from the committed source workbook and the public EVE endpoints.

Repository files required:

```text
dataset/SUPP_TABLES_BRCA12_APR_2026.xlsx
dataset/ACMG_other_points.xlsx
dataset/eve/build_eve_artifacts.py
dataset/eve/README.md
sup_table_18_19.py
main.py
```

Python packages required:

```text
pandas
openpyxl
```

Recreate the EVE artifacts from the repository root:

```sh
python3 dataset/eve/build_eve_artifacts.py
```

This command downloads the EVE per-protein variant archives, verifies the SHA256 hashes recorded below, and writes:

```text
dataset/eve/BRCA1_HUMAN.EVE.variants.zip
dataset/eve/BRCA2_HUMAN.EVE.variants.zip
dataset/eve/EVE_BRCA12_scores.xlsx
```

If local TLS certificate verification fails for the EVE site, confirm the HTTPS URLs below and rerun:

```sh
python3 dataset/eve/build_eve_artifacts.py --insecure-tls
```

Then regenerate the integrated supplementary workbook:

```sh
python3 main.py
```

The pipeline writes:

```text
results/code_generated.xlsx
```

The generated `Sup Table 18` and `Sup Table 19` tabs are the EVE-recalculated integrated classification tables.

## Affected Tabs and Columns

Only the Comment 1 code path is modified.

Affected generated workbook tabs:

- `Sup Table 18`: point system integration for BRCA1 variants
- `Sup Table 19`: point system integration for BRCA2 variants

Column substitutions:

- `Alpha Missense Score` was replaced by `EVE Score`
- `Alpha missense classification` was replaced by `EVE classification`
- `ACMG PP3/BP4 in silico predictor points` is now computed from EVE classes instead of the AlphaMissense-derived column in `dataset/ACMG_other_points.xlsx`

Downstream columns affected by the substitution:

- `FINAL ACMG POINTS`
- `FINAL ClinVar Classification`
- `Notes`, only where the final class changes a reference-panel disregard note

Unchanged non-predictor ACMG evidence columns still loaded from `dataset/ACMG_other_points.xlsx`:

- `ACMG PM2 points`
- `ACMG PP1/BS4 segregation points`
- `ACMG BA1/BS1 allele freq points`

## Source Data

EVE project site:

```text
https://evemodel.org/
```

EVE bulk download page:

```text
https://evemodel.org/download/bulk
```

EVE per-protein download endpoints used here:

```text
https://evemodel.org/api/proteins/web_pid/BRCA1_HUMAN/download/?variants=True
https://evemodel.org/api/proteins/web_pid/BRCA2_HUMAN/download/?variants=True
```

EVE protein identifier list endpoint used to confirm protein IDs:

```text
https://evemodel.org/api/proteins/list/identifiers/
```

Identifier rows found:

```text
BRCA1, BRCA1_HUMAN, UniProt P38398
BRCA1, BRCA1_BRCT_HUMAN, UniProt P38398
BRCA1, BRCA1_RING_HUMAN, UniProt P38398
BRCA2, BRCA2_HUMAN, UniProt P51587
```

The implementation uses the full-protein entries:

- `BRCA1_HUMAN`
- `BRCA2_HUMAN`

The BRCA1 BRCT/RING submodels are not used because the integration tables cover full-length BRCA1 missense variants.

Primary EVE publication:

```text
Frazer et al. Disease variant prediction with deep generative models of evolutionary data. Nature 599, 91-95 (2021).
https://doi.org/10.1038/s41586-021-04043-8
```

Reviewer benchmark paper:

```text
Livesey and Marsh. Variant effect predictor correlation with functional assays is reflective of clinical classification performance. Genome Biology 26, 104 (2025).
https://doi.org/10.1186/s13059-025-03575-w
```

## Local Files Created Or Recreated

The following EVE files are created under `dataset/eve/` by `dataset/eve/build_eve_artifacts.py`. They are ignored by Git and are not exposed in the GitHub repository:

| file | size in bytes | SHA256 |
| --- | ---: | --- |
| `dataset/eve/BRCA1_HUMAN.EVE.variants.zip` | 10,219,048 | `d99a1f4b383154afdec9cca35e5a27c91184f63e41a01805461a9b8d280a0b39` |
| `dataset/eve/BRCA2_HUMAN.EVE.variants.zip` | 16,023,212 | `ba15c69673dab8bb6ce73f18407e10bf4106eebdafdf84caaea9db39a01038d8` |
| `dataset/eve/EVE_BRCA12_scores.xlsx` | 1,804,353 | `c5e9b06fe5d03e5e6cb4603ee83f640f0f5e45d840910f773d0afcf224c0fb59` |

The downloaded zip archives contain:

| archive | member | uncompressed bytes |
| --- | --- | ---: |
| `BRCA1_HUMAN.EVE.variants.zip` | `BRCA1_HUMAN.csv` | 10,218,920 |
| `BRCA2_HUMAN.EVE.variants.zip` | `BRCA2_HUMAN.csv` | 16,023,084 |

## Download And Build Commands

Preferred command:

```sh
python3 dataset/eve/build_eve_artifacts.py
```

Manual download commands, equivalent to the download step in the script:

```sh
mkdir -p dataset/eve

curl -ksS "https://evemodel.org/api/proteins/web_pid/BRCA1_HUMAN/download/?variants=True" \
  -o dataset/eve/BRCA1_HUMAN.EVE.variants.zip

curl -ksS "https://evemodel.org/api/proteins/web_pid/BRCA2_HUMAN/download/?variants=True" \
  -o dataset/eve/BRCA2_HUMAN.EVE.variants.zip
```

`curl -k` was used because the local Python TLS stack rejected the certificate chain in this environment. The saved source URLs are HTTPS EVE endpoints.

## Raw EVE Schema

The source CSV files contain many fields. The columns used by this implementation are:

```text
wt_aa
position
mt_aa
EVE_scores_ASM
EVE_classes_75_pct_retained_ASM
```

The code maps them to:

| EVE source column | normalized column |
| --- | --- |
| `wt_aa` | `T2` |
| `position` | `T3` |
| `mt_aa` | `T4` |
| `EVE_scores_ASM` | `EVE Score` |
| `EVE_classes_75_pct_retained_ASM` | `EVE classification` |

ASM was selected because it is the primary EVE score/class set used by the EVE site and prior benchmark reuse. The alternate BPU columns did not fill any ASM missing-score rows for BRCA1 or BRCA2.

## Normalized Workbook Generation

The normalized workbook is:

```text
dataset/eve/EVE_BRCA12_scores.xlsx
```

It contains:

- `README`
- `BRCA1`
- `BRCA2`

The `BRCA1` and `BRCA2` sheets contain:

```text
INDEX
T1
T2
T3
T4
T5
T6
T7
EVE Score
EVE classification
EVE Source
```

Generation steps:

1. Run `dataset/eve/build_eve_artifacts.py` from the repository root.
2. The script downloads or reuses `BRCA1_HUMAN.EVE.variants.zip` and `BRCA2_HUMAN.EVE.variants.zip`.
3. The script verifies the downloaded archive SHA256 hashes by default.
4. The script reads `dataset/SUPP_TABLES_BRCA12_APR_2026.xlsx`.
5. It reads `Sup Table 1` for BRCA1 and `Sup Table 2` for BRCA2 with `header=1`.
6. It reads each EVE source zip member CSV.
7. It keeps `wt_aa`, `position`, `mt_aa`, `EVE_scores_ASM`, and `EVE_classes_75_pct_retained_ASM`.
8. It renames `wt_aa`, `position`, and `mt_aa` to `T2`, `T3`, and `T4`.
9. It left-joins each master table to EVE by `T2`, `T3`, and `T4`.
10. It preserves master-table row order and adds a 1-based `INDEX`.
11. It adds `EVE Source` as `<zip file>:<csv member>`.
12. It writes `dataset/eve/EVE_BRCA12_scores.xlsx`.

The implementation of these steps is in:

```text
dataset/eve/build_eve_artifacts.py
```

## Coverage

Raw EVE source coverage:

| gene | EVE source rows | score non-null | class non-null |
| --- | ---: | ---: | ---: |
| BRCA1 | 37,260 | 35,245 | 35,245 |
| BRCA2 | 68,360 | 50,046 | 50,046 |

Coverage after joining to the full BRCA master tables:

| gene | master rows | scored rows | missing rows | score coverage |
| --- | ---: | ---: | ---: | ---: |
| BRCA1 | 11,009 | 10,958 | 51 | 99.54% |
| BRCA2 | 20,169 | 15,542 | 4,627 | 77.06% |

Coverage in the integrated output tabs:

| tab | rows | scored rows | missing rows | score coverage |
| --- | ---: | ---: | ---: | ---: |
| `Sup Table 18` | 3,247 | 3,218 | 29 | 99.11% |
| `Sup Table 19` | 6,177 | 5,601 | 576 | 90.68% |

The integrated no-score rows are handled as neutral in silico evidence:

```text
EVE classification = No score
ACMG PP3/BP4 in silico predictor points = 0
```

Integrated no-score rows by reference-panel label:

| tab | no-score rows | missing reference label | benign reference label | pathogenic reference label |
| --- | ---: | ---: | ---: | ---: |
| `Sup Table 18` | 29 | 25 | 4 | 0 |
| `Sup Table 19` | 576 | 556 | 20 | 0 |

No integrated pathogenic reference-panel variants lack EVE scores.

## EVE Class to ACMG Point Mapping

No custom clinical thresholds are used.

The code uses EVE's author-provided `EVE_classes_75_pct_retained_ASM` classes:

| EVE class | `ACMG PP3/BP4 in silico predictor points` |
| --- | ---: |
| `Benign` | -1 |
| `Pathogenic` | 1 |
| `Uncertain` | 0 |
| missing score/class | 0 |

This is implemented in `sup_table_18_19.py` by `_eve_points`.

## Code Changes

### sup_table_18_19.py

Changes:

- Replaced the CPT-1/AlphaMissense predictor path with EVE.
- Added EVE columns:
  - `EVE Score`
  - `EVE classification`
  - `EVE Source`
- Added `_load_eve_predictor`.
- Added `_eve_points`.
- Removed CPT-1 thresholding from the active implementation.
- Kept `ACMG PP3/BP4 in silico predictor points` as the integration column, but now computes it from EVE classes.
- Stopped loading `ACMG PP3/BP4 in silico predictor points` from `dataset/ACMG_other_points.xlsx` to avoid retaining the old AlphaMissense-derived points.
- The command-line option remains `--predictor`; `--alpha` is still accepted only as a deprecated compatibility alias.

### main.py

The Sup Table 18/19 writer now receives:

```text
dataset/eve/EVE_BRCA12_scores.xlsx
```

## Generated Output Checks

The full workbook was regenerated with:

```sh
python3 main.py
```

Output workbook:

```text
results/code_generated.xlsx
```

The generated `Sup Table 18` and `Sup Table 19` no longer contain:

- `Alpha Missense Score`
- `Alpha missense classification`
- CPT-1 columns

They now contain:

- `EVE Score`
- `EVE classification`
- `ACMG PP3/BP4 in silico predictor points`

Final generated EVE counts:

| tab | rows | EVE benign | EVE uncertain | EVE pathogenic | EVE no score | PP3/BP4 -1 | PP3/BP4 0 | PP3/BP4 1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `Sup Table 18` | 3,247 | 650 | 936 | 1,632 | 29 | 650 | 965 | 1,632 |
| `Sup Table 19` | 6,177 | 1,363 | 2,449 | 1,789 | 576 | 1,363 | 3,025 | 1,789 |

Before and after final classification counts:

| tab | source | B | LB | VUS | LP | P |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `Sup Table 18` | previous AlphaMissense integration | 83 | 2,211 | 626 | 318 | 9 |
| `Sup Table 18` | new EVE integration | 83 | 2,207 | 583 | 364 | 10 |
| `Sup Table 19` | previous AlphaMissense integration | 107 | 4,595 | 1,119 | 351 | 5 |
| `Sup Table 19` | new EVE integration | 104 | 3,877 | 1,860 | 331 | 5 |

Row-level comparison against the previous April workbook:

| tab | row set changed? | PP3/BP4 point changes | final point changes | final classification changes |
| --- | --- | ---: | ---: | ---: |
| `Sup Table 18` | no | 1,872 | 1,872 | 70 |
| `Sup Table 19` | no | 3,856 | 3,856 | 874 |

Final classification transitions caused by the EVE substitution:

| tab | transition | count |
| --- | --- | ---: |
| `Sup Table 18` | VUS -> LP | 55 |
| `Sup Table 18` | LP -> VUS | 8 |
| `Sup Table 18` | LB -> VUS | 5 |
| `Sup Table 18` | VUS -> LB | 1 |
| `Sup Table 18` | LP -> P | 1 |
| `Sup Table 19` | LB -> VUS | 727 |
| `Sup Table 19` | LP -> VUS | 79 |
| `Sup Table 19` | VUS -> LP | 59 |
| `Sup Table 19` | VUS -> LB | 6 |
| `Sup Table 19` | B -> LB | 3 |

## Verification Commands

Syntax check:

```sh
python3 -m py_compile dataset/eve/build_eve_artifacts.py sup_table_18_19.py main.py
```

EVE artifact generation check using cached downloads:

```sh
python3 dataset/eve/build_eve_artifacts.py \
  --skip-download \
  --output-workbook /private/tmp/EVE_BRCA12_scores_test.xlsx
```

Isolated Sup Table 18/19 writer check:

```sh
python3 sup_table_18_19.py dataset/SUPP_TABLES_BRCA12_APR_2026.xlsx \
  -o /private/tmp/sup18_19_eve_test.xlsx \
  --predictor /private/tmp/EVE_BRCA12_scores_test.xlsx \
  --other-points dataset/ACMG_other_points.xlsx
```

Full workbook regeneration:

```sh
python3 main.py
```
