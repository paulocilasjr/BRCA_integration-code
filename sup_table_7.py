from pathlib import Path
from typing import Dict, Tuple, Union

import pandas as pd
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font

DATA_START_COL = "T8"
REFERENCE_COL = "T6"
DOCUMENTED_COL = "T7"
CUTOFFS = [1, 2, 10, 15, 20, 25, 30, 35, 40]
VUS_CUTOFFS = [1, 2, 5, 10, 15, 20]


def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = out.columns.astype(str).str.strip()
    return out


def _get_data_block(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = _normalize_columns(df)
    if DATA_START_COL not in df.columns:
        raise KeyError(f"Missing data start column: {DATA_START_COL}")
    if REFERENCE_COL not in df.columns:
        raise KeyError(f"Missing reference column: {REFERENCE_COL}")
    if DOCUMENTED_COL not in df.columns:
        raise KeyError(f"Missing documented column: {DOCUMENTED_COL}")
    start_idx = df.columns.get_loc(DATA_START_COL)
    data_df = df.iloc[:, start_idx:]
    return df, data_df


def _clean_reference_col(series: pd.Series) -> pd.Series:
    return series.replace(r"^\s*$", pd.NA, regex=True)


def _row_data_counts(df: pd.DataFrame) -> pd.Series:
    _, data_df = _get_data_block(df)
    return data_df.notna().sum(axis=1)


def total_data_points(df: pd.DataFrame) -> int:
    return int(_row_data_counts(df).sum())


def total_variants_tested(df: pd.DataFrame) -> int:
    return int(_row_data_counts(df).gt(0).sum())


def total_vus_variants_tested(df: pd.DataFrame) -> int:
    df = _normalize_columns(df)
    tested = _row_data_counts(df).gt(0)
    t6 = _clean_reference_col(df[REFERENCE_COL])
    documented = pd.to_numeric(df[DOCUMENTED_COL], errors="coerce").fillna(0).eq(1)
    return int((tested & t6.isna() & documented).sum())


def total_reference_variants_tested(df: pd.DataFrame) -> int:
    df = _normalize_columns(df)
    tested = _row_data_counts(df).gt(0)
    t6 = _clean_reference_col(df[REFERENCE_COL])
    return int((tested & t6.notna()).sum())


def summarize_table(df: pd.DataFrame) -> Dict[str, Union[int, Dict[int, int]]]:
    df = _normalize_columns(df)
    data_counts = _row_data_counts(df)
    total_rows = int(len(df))
    t7 = pd.to_numeric(df[DOCUMENTED_COL], errors="coerce").fillna(0)
    reported = t7.eq(1)
    t6 = _clean_reference_col(df[REFERENCE_COL])
    is_reference = t6.notna()
    not_tested = int(data_counts.eq(0).sum())
    tested_ge = {cutoff: int(data_counts.ge(cutoff).sum()) for cutoff in CUTOFFS}
    reported_total = int(reported.sum())
    reported_tested = int((reported & data_counts.gt(0)).sum())
    reported_vus_tested = int(
        (reported & data_counts.gt(0) & t6.isna()).sum()
    )
    reference_total = int(is_reference.sum())
    reference_not_tested = int((is_reference & data_counts.eq(0)).sum())
    reference_tested_ge = {
        cutoff: int((is_reference & data_counts.ge(cutoff)).sum()) for cutoff in CUTOFFS
    }
    vus_total = int(t6.isna().sum())
    vus_not_tested = int((t6.isna() & data_counts.eq(0)).sum())
    vus_tested_ge = {
        cutoff: int((t6.isna() & data_counts.ge(cutoff)).sum()) for cutoff in VUS_CUTOFFS
    }

    return {
        "total_rows": total_rows,
        "not_tested": not_tested,
        "tested_ge": tested_ge,
        "total_data_points": int(data_counts.sum()),
        "total_variants_tested": int(data_counts.gt(0).sum()),
        "total_vus_variants_tested": total_vus_variants_tested(df),
        "total_reference_variants_tested": total_reference_variants_tested(df),
        "reported_total": reported_total,
        "reported_tested": reported_tested,
        "reported_vus_tested": reported_vus_tested,
        "reference_total": reference_total,
        "reference_not_tested": reference_not_tested,
        "reference_tested_ge": reference_tested_ge,
        "vus_total": vus_total,
        "vus_not_tested": vus_not_tested,
        "vus_tested_ge": vus_tested_ge,
    }


def summarize_tables(
    brca1_table: pd.DataFrame, brca2_table: pd.DataFrame
) -> Dict[str, Dict[str, Union[int, Dict[int, int]]]]:
    return {
        "BRCA1": summarize_table(brca1_table),
        "BRCA2": summarize_table(brca2_table),
    }


def write_sup_table_7(
    summary: Dict[str, Dict[str, Union[int, Dict[int, int]]]],
    output_path: str,
) -> None:
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if output_file.exists():
        wb = load_workbook(output_file)
        if "Sup Table 7" in wb.sheetnames:
            wb.remove(wb["Sup Table 7"])
        ws = wb.create_sheet("Sup Table 7")
    else:
        wb = Workbook()
        ws = wb.active
        ws.title = "Sup Table 7"

    ws["A1"] = (
        "Supplementary Table 7: Summary statistics for tested BRCA1 and BRCA2 variants "
        "(number of individual instances of an assay performed for each variant)"
    )
    ws["A1"].font = Font(bold=True)
    ws["A1"].alignment = Alignment(horizontal="left", vertical="center", wrap_text=False)

    ws["A4"] = "Total number of data points (individual functional tests)"
    ws["A5"] = "Total number of variants tested (ref + VUS)"
    ws["A6"] = "Total number of reported VUS tested"
    ws["A7"] = "Total number of reference variants tested"

    ws["B3"] = "BRCA1"
    ws["E3"] = "BRCA2"

    ws["B4"] = summary["BRCA1"]["total_data_points"]
    ws["B5"] = summary["BRCA1"]["total_variants_tested"]
    ws["B6"] = summary["BRCA1"]["total_vus_variants_tested"]
    ws["B7"] = summary["BRCA1"]["total_reference_variants_tested"]

    ws["E4"] = summary["BRCA2"]["total_data_points"]
    ws["E5"] = summary["BRCA2"]["total_variants_tested"]
    ws["E6"] = summary["BRCA2"]["total_vus_variants_tested"]
    ws["E7"] = summary["BRCA2"]["total_reference_variants_tested"]

    ws["A9"] = "Possible missense variants"
    ws["B9"] = "N"
    ws["C9"] = "%"
    ws["E9"] = "N"
    ws["F9"] = "%"

    ws["A10"] = "Total number of unique missense variants resulting from single nucleotide change"
    ws["A11"] = "Not yet test"
    ws["B10"] = summary["BRCA1"]["total_rows"]
    ws["E10"] = summary["BRCA2"]["total_rows"]
    ws["C10"] = 100
    ws["F10"] = 100

    def percent(count: int, total: int) -> float:
        return 0.0 if total == 0 else round((count / total) * 100, 2)

    brca1_total_rows = summary["BRCA1"]["total_rows"]
    brca2_total_rows = summary["BRCA2"]["total_rows"]

    ws["B11"] = summary["BRCA1"]["not_tested"]
    ws["E11"] = summary["BRCA2"]["not_tested"]
    ws["C11"] = percent(summary["BRCA1"]["not_tested"], brca1_total_rows)
    ws["F11"] = percent(summary["BRCA2"]["not_tested"], brca2_total_rows)

    start_row = 12
    for offset, cutoff in enumerate(CUTOFFS):
        row = start_row + offset
        ws[f"A{row}"] = f"Tested >= {cutoff} times"
        brca1_val = summary["BRCA1"]["tested_ge"][cutoff]
        brca2_val = summary["BRCA2"]["tested_ge"][cutoff]
        ws[f"B{row}"] = brca1_val
        ws[f"E{row}"] = brca2_val
        ws[f"C{row}"] = percent(brca1_val, brca1_total_rows)
        ws[f"F{row}"] = percent(brca2_val, brca2_total_rows)

    ws["A22"] = "Reported missense variants (BRCA Exchange)"
    ws["B22"] = "N"
    ws["C22"] = "%"
    ws["E22"] = "N"
    ws["F22"] = "%"

    ws["A23"] = "Total number of reported missense variants"
    ws["A24"] = "Reported missense variants tested"
    ws["A25"] = "Reported missense VUS tested"

    ws["B23"] = summary["BRCA1"]["reported_total"]
    ws["E23"] = summary["BRCA2"]["reported_total"]
    ws["C23"] = 100
    ws["F23"] = 100

    ws["B24"] = summary["BRCA1"]["reported_tested"]
    ws["E24"] = summary["BRCA2"]["reported_tested"]
    ws["C24"] = percent(summary["BRCA1"]["reported_tested"], summary["BRCA1"]["reported_total"])
    ws["F24"] = percent(summary["BRCA2"]["reported_tested"], summary["BRCA2"]["reported_total"])

    ws["B25"] = summary["BRCA1"]["reported_vus_tested"]
    ws["E25"] = summary["BRCA2"]["reported_vus_tested"]
    ws["C25"] = percent(summary["BRCA1"]["reported_vus_tested"], summary["BRCA1"]["reported_total"])
    ws["F25"] = percent(summary["BRCA2"]["reported_vus_tested"], summary["BRCA2"]["reported_total"])

    ws["A27"] = "Reference Variants"
    ws["B27"] = "N"
    ws["C27"] = "%"
    ws["E27"] = "N"
    ws["F27"] = "%"

    ws["A28"] = "Total number of reference missense variants in reference panel"
    ws["A29"] = "Not yet tested"

    ws["B28"] = summary["BRCA1"]["reference_total"]
    ws["E28"] = summary["BRCA2"]["reference_total"]
    ws["C28"] = 100
    ws["F28"] = 100

    ws["B29"] = summary["BRCA1"]["reference_not_tested"]
    ws["E29"] = summary["BRCA2"]["reference_not_tested"]
    ws["C29"] = percent(
        summary["BRCA1"]["reference_not_tested"], summary["BRCA1"]["reference_total"]
    )
    ws["F29"] = percent(
        summary["BRCA2"]["reference_not_tested"], summary["BRCA2"]["reference_total"]
    )

    ref_start_row = 30
    for offset, cutoff in enumerate(CUTOFFS):
        row = ref_start_row + offset
        ws[f"A{row}"] = f"Tested >= {cutoff} times"
        brca1_val = summary["BRCA1"]["reference_tested_ge"][cutoff]
        brca2_val = summary["BRCA2"]["reference_tested_ge"][cutoff]
        ws[f"B{row}"] = brca1_val
        ws[f"E{row}"] = brca2_val
        ws[f"C{row}"] = percent(brca1_val, summary["BRCA1"]["reference_total"])
        ws[f"F{row}"] = percent(brca2_val, summary["BRCA2"]["reference_total"])

    ws["A40"] = "VUS only (excluding reference variants)"
    ws["B40"] = "N"
    ws["C40"] = "%"
    ws["E40"] = "N"
    ws["F40"] = "%"

    ws["A41"] = "Total number of possible VUS"
    ws["A42"] = "Not yet test"

    ws["B41"] = summary["BRCA1"]["vus_total"]
    ws["E41"] = summary["BRCA2"]["vus_total"]
    ws["C41"] = 100
    ws["F41"] = 100

    ws["B42"] = summary["BRCA1"]["vus_not_tested"]
    ws["E42"] = summary["BRCA2"]["vus_not_tested"]
    ws["C42"] = percent(summary["BRCA1"]["vus_not_tested"], summary["BRCA1"]["vus_total"])
    ws["F42"] = percent(summary["BRCA2"]["vus_not_tested"], summary["BRCA2"]["vus_total"])

    vus_start_row = 43
    for offset, cutoff in enumerate(VUS_CUTOFFS):
        row = vus_start_row + offset
        ws[f"A{row}"] = f"VUS tested >= {cutoff} times"
        brca1_val = summary["BRCA1"]["vus_tested_ge"][cutoff]
        brca2_val = summary["BRCA2"]["vus_tested_ge"][cutoff]
        ws[f"B{row}"] = brca1_val
        ws[f"E{row}"] = brca2_val
        ws[f"C{row}"] = percent(brca1_val, summary["BRCA1"]["vus_total"])
        ws[f"F{row}"] = percent(brca2_val, summary["BRCA2"]["vus_total"])

    ws["B3"].font = Font(bold=True)
    ws["E3"].font = Font(bold=True)

    for row in (9, 22, 27, 40):
        for col in ("A", "B", "C", "D", "E", "F"):
            cell = ws[f"{col}{row}"]
            if cell.value is not None:
                cell.font = Font(bold=True)

    for row in ws.iter_rows():
        for cell in row:
            if cell.value is None:
                continue
            cell.font = Font(name="Arial", size=10, bold=bool(cell.font.bold))
            if cell.column_letter == "A":
                if cell.row == 1:
                    cell.alignment = Alignment(horizontal="left", vertical="center", wrap_text=False)
                else:
                    cell.alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
            else:
                cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

    for col in ws.columns:
        max_len = 0
        col_letter = col[0].column_letter
        for cell in col:
            if cell.value is None:
                continue
            if col_letter == "A" and cell.row == 1:
                continue
            max_len = max(max_len, len(str(cell.value)))
        if max_len > 0:
            ws.column_dimensions[col_letter].width = max(12, min(max_len + 2, 120))

    wb.save(output_file)
