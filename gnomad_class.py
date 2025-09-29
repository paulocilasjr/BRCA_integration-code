#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Offline scoring using local gnomAD export files, with robust header handling.

Inputs:
  - SUPP_TABLES_BRCA12_OCT_2025.xlsx
    * For Sup Table 19 and Sup Table 20, headers are on the 2nd row (index 1).
    * Other sheets use the first row as headers (index 0).
  - BRCA1_gnomAD_v4.1.0_ENSG00000012048_*.csv   (tsv/csv export)
  - BRCA2_gnomAD_v4.1.0_ENSG00000139618_*.csv   (tsv/csv export)

Outputs:
  - SUPP_TABLES_BRCA12_OCT_2025_scored.xlsx

Scoring:
  - PM2 is set to 1 if non-founder FAF95 is NaN or 0, else 0.
  - BA1 is mapped from non-founder FAF95:
      > 0.1%  ( > 0.001 )   => -8
      > 0.01% ( > 0.0001 )  => -4
      > 0.002% ( > 0.00002 ) and <= 0.01% => -1
      otherwise => 0
  - Final = PM2 + BA1 + "Alpha Missense - ACMG" + "Capped Final Score"
  - Only the columns PM2, BA1, Final are appended to the target sheets.
"""

import re
import csv
import io
import pandas as pd
import numpy as np
from pathlib import Path

# ========= USER SETTINGS =========
INPUT_XLSX  = "SUPP_TABLES_BRCA12_OCT_2025.xlsx"
OUTPUT_XLSX = "SUPP_TABLES_BRCA12_OCT_2025_scored.xlsx"

SHEETS_TO_PROCESS = {
    "Sup Table 19": {
        "gene": "BRCA1",
        "gnomad_file": "BRCA1_gnomAD_v4.1.0_ENSG00000012048_2025_09_18_14_35_11.csv",
    },
    "Sup Table 20": {
        "gene": "BRCA2",
        "gnomad_file": "BRCA2_gnomAD_v4.1.0_ENSG00000139618_2025_09_18_14_36_35.csv",
    },
}

HEADER_ROW_INDEX = 1  # use row 2 as header ONLY for Sup Table 19/20
PROTEIN_CHANGE_COL = "T5"  # protein change column in sheets 19/20

# BA1 thresholds (FAF95)
THRESH_A = 0.001        # > 0.1%  -> -8
THRESH_B = 0.0001       # > 0.01% -> -4
THRESH_C_LOW  = 0.00002 # > 0.002%
THRESH_C_HIGH = 0.0001  # <= 0.01%

# Non-founder groups in v4 exports.
NON_FOUNDER_GROUPS = {
    "afr", "amr", "eas", "nfe", "sas", "oth", "other", "rem", "remaining", "mid", "middle eastern"
}
FOUNDER_GROUPS = {"fin", "asj", "amish"}

# ========= Protein notation normalization =========
AA3_TO_1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G",
    "His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S",
    "Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Sec":"U","Pyl":"O","Ter":"*"
}
AA3_RE = re.compile(r"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Sec|Pyl|Ter)")

def to_one_letter_protein(hgvsp: str) -> str:
    """
    Convert protein consequence to a consistent one-letter, no 'p.' prefix.
    """
    if not hgvsp or str(hgvsp).strip() == "":
        return ""
    s = str(hgvsp).strip()
    if ":" in s:
        s = s.split(":")[-1]
    if s.lower().startswith("p."):
        s = s[2:]
    def repl(m):
        three = m.group(1)
        return AA3_TO_1.get(three, three)
    s = AA3_RE.sub(repl, s)
    return s.strip()

# ========= Column resolution helpers =========
NBSP = "\u00A0"

def norm_colname(name: str) -> str:
    if name is None:
        return ""
    s = str(name).replace(NBSP, " ")
    s = re.sub(r"\s+", " ", s).strip().lower()
    return s

def resolve_columns(df: pd.DataFrame) -> dict:
    """
    Map flexible aliases -> actual columns present in df.
    Returns dict with keys: 'protein', 'faf_group', 'faf_freq'
    """
    actual = {norm_colname(c): c for c in df.columns}

    aliases = {
        "protein": [
            "protein consequence", "protein_consequence", "hgvsp", "hgvsp (protein)",
            "hgvs protein", "protein"
        ],
        "faf_group": [
            "groupmax faf group", "group max faf group", "groupmax_faf group", "groupmax_faf_group",
            "groupmax faf population", "faf groupmax group", "faf groupmax population"
        ],
        "faf_freq": [
            "groupmax faf frequency", "group max faf frequency", "groupmax_faf frequency",
            "groupmax_faf_frequency", "groupmax faf95 frequency", "groupmax faf95"
        ],
    }

    out = {}
    missing = []
    for key, candidates in aliases.items():
        found = None
        for cand in candidates:
            if cand in actual:
                found = actual[cand]
                break
        if not found:
            missing.append((key, candidates, list(df.columns)))
        else:
            out[key] = found

    if missing:
        msgs = []
        for key, cands, cols in missing:
            msgs.append(
                f"- Could not resolve '{key}'. Tried: {cands}. "
                f"Available columns:\n  {cols}"
            )
        raise KeyError("Header resolution failed:\n" + "\n".join(msgs))

    return out

def resolve_exact_ci(df: pd.DataFrame, expected_name: str) -> str | None:
    """
    Find a column by case-insensitive, whitespace-normalized exact name.
    Returns the actual column name or None.
    """
    want = norm_colname(expected_name)
    for c in df.columns:
        if norm_colname(c) == want:
            return c
    return None

# ========= gnomAD file loading & mapping =========
def sniff_delimiter(sample_bytes: bytes) -> str:
    """Use csv.Sniffer to guess the delimiter; default to tab if unsure."""
    try:
        sample = sample_bytes.decode("utf-8", errors="replace")
        dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";", "|"])
        return dialect.delimiter
    except Exception:
        return "\t"

def read_gnomad_table(path: Path) -> pd.DataFrame:
    """
    Read the gnomAD export, auto-detecting delimiter and normalizing headers.
    """
    raw = path.read_bytes()
    delim = sniff_delimiter(raw)
    df = pd.read_csv(io.BytesIO(raw), sep=delim, dtype=str, na_values=["", "NA", "NaN"], engine="python")

    # Normalize header names (replace NBSP, collapse spaces)
    new_cols = []
    for c in df.columns:
        c2 = c.replace(NBSP, " ")
        c2 = re.sub(r"\s+", " ", c2).strip()
        new_cols.append(c2)
    df.columns = new_cols
    return df

def build_faf_map(gnomad_df: pd.DataFrame) -> dict:
    """
    Build { normalized_one_letter_protein -> max non-founder FAF }.
    Uses 'GroupMax FAF group' and 'GroupMax FAF frequency' ONLY if group is non-founder.
    """
    # Resolve columns robustly
    cols = resolve_columns(gnomad_df)

    # Coerce to numeric FAF
    gdf = gnomad_df.copy()
    gdf["_faf"] = pd.to_numeric(gdf[cols["faf_freq"]], errors="coerce")

    # Normalize groups
    grp = gdf[cols["faf_group"]].fillna("").astype(str).str.strip().str.lower()

    # Keep rows where GroupMax group is a non-founder
    is_nf = grp.isin(NON_FOUNDER_GROUPS)
    df_nf = gdf.loc[is_nf & gdf["_faf"].notna(), [cols["protein"], "_faf"]].copy()

    # Normalize protein to one-letter key
    df_nf["_prot_key"] = df_nf[cols["protein"]].map(to_one_letter_protein).str.upper()

    # Aggregate: choose MAX FAF for identical protein consequence keys
    faf_map = (
        df_nf.groupby("_prot_key")["_faf"]
        .max()
        .to_dict()
    )
    return faf_map

# ========= Excel scoring (adds only PM2, BA1, Final) =========
def score_sheet(df: pd.DataFrame, gene: str, faf_map: dict) -> pd.DataFrame:
    """
    Compute PM2 and BA1 from non-founder FAF95 and then:
      Final = PM2 + BA1 + "Alpha Missense - ACMG" + "Capped Final Score"

    Only appends PM2, BA1, Final to the DataFrame (in this order).
    """
    # Ensure T5 exists (protein change)
    if PROTEIN_CHANGE_COL not in df.columns:
        raise KeyError(f"Expected protein change column '{PROTEIN_CHANGE_COL}' not found in this sheet.")

    # Locate the two required pre-existing columns
    alpha_col = resolve_exact_ci(df, "Alpha Missense - ACMG")
    capped_col = resolve_exact_ci(df, "Capped Final Score")
    if alpha_col is None or capped_col is None:
        missing = []
        if alpha_col is None:
            missing.append("'Alpha Missense - ACMG'")
        if capped_col is None:
            missing.append("'Capped Final Score'")
        raise KeyError(f"Missing required column(s): {', '.join(missing)}")

    # Normalize T5 to one-letter keys (no 'p.')
    prot_keys = df[PROTEIN_CHANGE_COL].fillna("").astype(str).map(to_one_letter_protein).str.upper()

    # Lookup FAF from map (do not persist as a column)
    faf_vals = prot_keys.map(lambda k: np.nan if not k else faf_map.get(k, np.nan))
    faf_series = pd.to_numeric(faf_vals, errors="coerce")

    # PM2: NaN or 0 -> 1, else 0
    pm2 = np.where(faf_series.isna() | (faf_series == 0), 1, 0)

    # BA1 thresholds
    faf_filled = faf_series.fillna(0)
    cond_a = faf_filled > THRESH_A
    cond_b = (~cond_a) & (faf_filled > THRESH_B)
    cond_c = (~cond_a) & (~cond_b) & (faf_filled > THRESH_C_LOW) & (faf_filled <= THRESH_C_HIGH)
    ba1 = np.select([cond_a, cond_b, cond_c], [-8, -4, -1], default=0)

    # Existing numeric columns
    alpha_vals  = pd.to_numeric(df[alpha_col], errors="coerce").fillna(0)
    capped_vals = pd.to_numeric(df[capped_col], errors="coerce").fillna(0)

    final_vals = pm2 + ba1 + alpha_vals + capped_vals

    # Append ONLY the requested columns, in order
    df["PM2"] = pm2
    df["BA1"] = ba1
    df["Final"] = final_vals

    return df

# ========= MAIN =========
def main():
    # 1) Load and build FAF maps from the gnomAD files
    faf_maps = {}
    for sheet_name, meta in SHEETS_TO_PROCESS.items():
        gpath = Path(meta["gnomad_file"])
        if not gpath.exists():
            raise FileNotFoundError(f"Missing gnomAD file for {meta['gene']}: {gpath}")
        gdf = read_gnomad_table(gpath)
        faf_maps[meta["gene"]] = build_faf_map(gdf)

    # 2) Open Excel, process the sheets
    in_path = Path(INPUT_XLSX)
    if not in_path.exists():
        raise FileNotFoundError(f"Input workbook not found: {INPUT_XLSX}")

    xls = pd.ExcelFile(in_path, engine="openpyxl")
    out_sheets = {}

    for sheet in xls.sheet_names:
        # Only Sup Table 19/20 use header row index 1; others use 0
        header_to_use = HEADER_ROW_INDEX if sheet in SHEETS_TO_PROCESS else 0
        df = pd.read_excel(in_path, sheet_name=sheet, header=header_to_use, engine="openpyxl")

        # Score ONLY the target sheets (19/20)
        if sheet in SHEETS_TO_PROCESS:
            gene = SHEETS_TO_PROCESS[sheet]["gene"]
            try:
                df = score_sheet(df, gene=gene, faf_map=faf_maps[gene])
            except KeyError as e:
                # If required columns are missing, keep the sheet unmodified but warn
                print(f"⚠️ Skipping scoring for '{sheet}' ({gene}): {e}")

        out_sheets[sheet] = df

    # 3) Save output
    with pd.ExcelWriter(OUTPUT_XLSX, engine="xlsxwriter") as writer:
        for name, df in out_sheets.items():
            df.to_excel(writer, sheet_name=name, index=False)

    print(f"✅ Done. Wrote: {OUTPUT_XLSX}")

if __name__ == "__main__":
    main()

