import pandas as pd

from sup_table_7 import summarize_tables, write_sup_table_7
from sup_table_8 import summarize_tables as summarize_tables_8, write_sup_table_8
from sup_table_9_10 import build_track_classification_map, write_sup_table_9_10
from supp_fig2 import generate_supp_fig2
from sup_table_11 import write_sup_table_11
from sup_table_12_13 import write_sup_table_12_13
from sup_table_14_15 import write_sup_table_14_15
from sup_table_18_19 import write_sup_tables_18_19
from sup_table_16 import (
    BRCA1_FEATURES,
    BRCA2_FEATURES,
    build_assignment_df,
    build_feature_table,
    load_features,
    write_sup_table_16,
)
from sup_table_17 import write_sup_table_17

FILE_PATH = "dataset/SUPP_TABLES_BRCA12_DEC_2025_WORKING_fixed.xlsx"
OUTPUT_PATH = "results/code_generated.xlsx"

SHEETS = {
    "Sup Table 1": "BRCA1_table",
    "Sup Table 2": "BRCA2_table",
    "Sup Table 3": "BRCA1_metadata",
    "Sup Table 4": "BRCA2_metadata",
    "Sup Table 5": "BRCA1_Reference_panel",
    "Sup Table 6": "BRCA2_Reference_panel",
}

def _norm_col(text: object) -> str:
    return " ".join(str(text).strip().lower().split())


def _find_ref_variant_col(df: pd.DataFrame) -> str:
    candidates = ["T5", "Variant", "Variant ID", "variant", "variant id", "T5 (Variant)"]
    mapping = {_norm_col(c): c for c in df.columns}
    for cand in candidates:
        key = _norm_col(cand)
        if key in mapping:
            return mapping[key]
    return str(df.columns[0])


def _find_include_col(df: pd.DataFrame) -> str | None:
    target = _norm_col("Include in Ref Panel")
    mapping = {_norm_col(c): c for c in df.columns}
    if target in mapping:
        return mapping[target]
    # fallback: any column containing both "include" and "ref panel"
    for col in df.columns:
        norm = _norm_col(col)
        if "include" in norm and "ref panel" in norm:
            return col
    return None


def _excluded_variants_from_ref_panel(df: pd.DataFrame) -> set[str]:
    include_col = _find_include_col(df)
    if not include_col:
        return set()
    var_col = _find_ref_variant_col(df)
    include_vals = df[include_col].astype(str).str.strip().str.lower()
    excluded = df.loc[include_vals.eq("no"), var_col].dropna().astype(str).str.strip()
    return set(excluded.tolist())


def _mask_t6_for_excluded_variants(variant_df: pd.DataFrame, excluded: set[str]) -> pd.DataFrame:
    if not excluded or "T5" not in variant_df.columns or "T6" not in variant_df.columns:
        return variant_df
    df = variant_df.copy()
    t5_vals = df["T5"].astype(str).str.strip()
    mask = t5_vals.isin(excluded)
    if mask.any():
        df.loc[mask, "T6"] = pd.NA
    return df



def load_tables(path: str) -> dict:
    tables = {}
    for sheet_name, var_name in SHEETS.items():
        tables[var_name] = pd.read_excel(path, sheet_name=sheet_name, header=1)
    # Apply ref-panel include/exclude mapping to T6 in Sup Table 1/2.
    if "BRCA1_Reference_panel" in tables and "BRCA1_table" in tables:
        excluded_b1 = _excluded_variants_from_ref_panel(tables["BRCA1_Reference_panel"])
        tables["BRCA1_table"] = _mask_t6_for_excluded_variants(tables["BRCA1_table"], excluded_b1)
    if "BRCA2_Reference_panel" in tables and "BRCA2_table" in tables:
        excluded_b2 = _excluded_variants_from_ref_panel(tables["BRCA2_Reference_panel"])
        tables["BRCA2_table"] = _mask_t6_for_excluded_variants(tables["BRCA2_table"], excluded_b2)
    return tables


def load_functional_table(path: str, preferred_sheet: str) -> pd.DataFrame:
    header_candidates = (1, 0, 2)
    sheet_candidates = [preferred_sheet]
    try:
        xls = pd.ExcelFile(path)
        sheet_candidates.extend([s for s in xls.sheet_names if s not in sheet_candidates])
    except Exception:
        pass

    def norm(text: object) -> str:
        return " ".join(str(text).replace("\u00a0", " ").strip().lower().split())

    def has_required_columns(df: pd.DataFrame) -> bool:
        cols = [norm(c) for c in df.columns]
        has_t5 = any(c == "t5" or c.startswith("t5") for c in cols)
        has_func = any("functional" in c and "category" in c for c in cols)
        return has_t5 and has_func

    for sheet in sheet_candidates:
        for header in header_candidates:
            try:
                df = pd.read_excel(path, sheet_name=sheet, header=header)
            except Exception:
                continue
            if has_required_columns(df):
                return df

        # Fallback: scan for header row containing Functional category (and ideally T5)
        try:
            raw = pd.read_excel(path, sheet_name=sheet, header=None, nrows=100)
        except Exception:
            continue
        header_row = None
        for idx, row in raw.iterrows():
            cells = [norm(c) for c in row.tolist()]
            has_func = any("functional" in c and "category" in c for c in cells)
            has_t5 = any(c == "t5" or c.startswith("t5") for c in cells)
            if has_func and (has_t5 or header_row is None):
                header_row = idx
                if has_t5:
                    break
        if header_row is not None:
            try:
                return pd.read_excel(path, sheet_name=sheet, header=header_row)
            except Exception:
                continue
    raise ValueError(f"Functional category table not found (preferred: {preferred_sheet}).")


if __name__ == "__main__":
    tables = load_tables(FILE_PATH)
    BRCA1_table = tables["BRCA1_table"]
    BRCA2_table = tables["BRCA2_table"]
    BRCA1_metadata = tables["BRCA1_metadata"]
    BRCA2_metadata = tables["BRCA2_metadata"]
    BRCA1_Reference_panel = tables["BRCA1_Reference_panel"]
    BRCA2_Reference_panel = tables["BRCA2_Reference_panel"]

    summary_7 = summarize_tables(BRCA1_table, BRCA2_table, BRCA1_metadata, BRCA2_metadata)
    write_sup_table_7(summary_7, OUTPUT_PATH)

    summary_8 = summarize_tables_8(BRCA1_table, BRCA2_table, BRCA1_metadata, BRCA2_metadata)
    write_sup_table_8(summary_8, OUTPUT_PATH)

    write_sup_table_9_10(
        BRCA1_table,
        BRCA1_metadata,
        BRCA2_table,
        BRCA2_metadata,
        OUTPUT_PATH,
    )
    generate_supp_fig2(input_path=OUTPUT_PATH, output_prefix="figures/supp_fig2")

    write_sup_table_11(
        BRCA1_table,
        BRCA2_table,
        OUTPUT_PATH,
        BRCA1_metadata,
        BRCA2_metadata,
    )

    brca1_class_map = build_track_classification_map(BRCA1_table, BRCA1_metadata, 11009)
    brca2_class_map = build_track_classification_map(BRCA2_table, BRCA2_metadata, 20169)
    sup12_df, sup13_df = write_sup_table_12_13(
        BRCA1_table,
        BRCA2_table,
        brca1_class_map,
        brca2_class_map,
        OUTPUT_PATH,
        BRCA1_metadata,
        BRCA2_metadata,
    )

    write_sup_table_14_15(
        BRCA1_table,
        BRCA2_table,
        BRCA1_Reference_panel,
        BRCA2_Reference_panel,
        sup12_df,
        sup13_df,
        "results/code_generated.xlsx",
    )

    brca1_features = load_features(None, BRCA1_FEATURES, "BRCA1")
    brca2_features = load_features(None, BRCA2_FEATURES, "BRCA2")
    brca1_assign = build_assignment_df(sup12_df)
    brca2_assign = build_assignment_df(sup13_df)
    brca1_tbl = build_feature_table(brca1_assign, brca1_features)
    brca2_tbl = build_feature_table(brca2_assign, brca2_features)
    write_sup_table_16(brca1_tbl, brca2_tbl, OUTPUT_PATH)

    write_sup_table_17(
        FILE_PATH,
        OUTPUT_PATH,
        sup12_sheet="Sup Table 12",
        sup13_sheet="Sup Table 13",
    )

    write_sup_tables_18_19(
        sup12_df,
        sup13_df,
        "dataset/AlphaMissense_Calculations_all.xlsx",
        "dataset/ACMG_other_points.xlsx",
        OUTPUT_PATH,
    )

    print("Sup Table 7, Sup Table 8, Sup Table 9, Sup Table 10, Sup Table 11, Sup Table 12, Sup Table 13, Sup Table 14, Sup Table 15, and Sup Table 16 written to results/code_generated.xlsx")
