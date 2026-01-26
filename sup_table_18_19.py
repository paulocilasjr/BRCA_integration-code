#!/usr/bin/env python3
"""
Build Supplementary Tables 18 and 19 from Sup Tables 12/13 and AlphaMissense/ACMG inputs.

Adds columns:
  - ACMG BS3/PS3 functional points computation
  - Number of assays
  - ACMG BS3/PS3 Final Functional Points
  - ACMG BS3/PS3 Capped Final Functional Points
  - Alpha Missense Score
  - Alpha missense classification
  - ACMG PP3/BP4 in silico predictor points
  - ACMG PM2 points
  - ACMG PP1/BS4 segregation points
  - ACMG BA1/BS1 allele freq points
  - FINAL ACMG POINTS
  - FINAL ClinVar Classification
"""

from __future__ import annotations

import argparse
import ast
import re
from pathlib import Path
from typing import Iterable, List, Tuple

import pandas as pd
from openpyxl import load_workbook


def _norm_col(c: str) -> str:
    return re.sub(r"\s+", " ", str(c).strip().lower())


def _find_col(df: pd.DataFrame, candidates: Iterable[str]) -> str:
    mapping = {_norm_col(c): c for c in df.columns}
    for cand in candidates:
        key = _norm_col(cand)
        if key in mapping:
            return mapping[key]
    for cand in candidates:
        key = _norm_col(cand)
        for k_norm, k_orig in mapping.items():
            if key in k_norm:
                return k_orig
    raise KeyError(f"Could not find any of columns: {list(candidates)}")


def resolve_sheet_name(xlsx_path: str, preferred: str, fallbacks: List[str]) -> str:
    wb = load_workbook(xlsx_path, read_only=True, data_only=True)
    names = wb.sheetnames
    if preferred in names:
        return preferred
    lower_map = {n.lower(): n for n in names}
    if preferred.lower() in lower_map:
        return lower_map[preferred.lower()]
    for pat in fallbacks:
        pat_l = pat.lower()
        for n in names:
            if pat_l in n.lower():
                return n
    raise ValueError(f"Sheet not found. preferred={preferred!r}. Available: {names}")


def parse_vote_to_score(vote: str) -> int:
    vote = vote.strip()
    if re.match(r"BS3_supporting\b", vote):
        return -1
    if re.match(r"BS3_moderate\b", vote):
        return -2
    if re.match(r"BS3\b", vote):
        return -4
    if re.match(r"PS3_supporting\b", vote):
        return 1
    if re.match(r"PS3_moderate\b", vote):
        return 2
    if re.match(r"PS3\b", vote):
        return 4
    if re.match(r"hypomorph\b", vote, flags=re.IGNORECASE):
        return 0
    return 0


def safe_parse_votes(votes_raw) -> List[str]:
    if pd.isna(votes_raw) or votes_raw == "":
        return []
    try:
        if isinstance(votes_raw, str):
            return ast.literal_eval(votes_raw)
        return votes_raw if isinstance(votes_raw, list) else []
    except Exception:
        return []


def compute_acmg_score(votes: List[str]) -> List[str]:
    scores = []
    for v in votes:
        score = parse_vote_to_score(v)
        match = re.search(r"T\d+", v)
        score_str = f"{score} ({match.group()})" if match else str(score)
        scores.append(score_str)
    return scores


def compute_final_score(score_list: List[str]) -> int:
    try:
        return sum(int(str(s).split()[0]) for s in score_list)
    except Exception:
        return 0


def cap_score(score: int) -> int:
    if score < -4:
        return -4
    if score > 4:
        return 4
    return score


def classify_final_points(score: float) -> str:
    if score >= 10:
        return "P"
    if 6 <= score <= 9:
        return "LP"
    if -1 <= score <= 5:
        return "VUS"
    if -6 <= score <= -2:
        return "LB"
    if score <= -7:
        return "B"
    return ""


def _load_alpha_missense(path: str, sheet: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=sheet, header=None)
    df = df.rename(columns={0: "T2", 1: "T3", 2: "T4", 3: "Alpha Missense Score", 4: "Alpha missense classification"})
    df = df.dropna(subset=["T2", "T3", "T4"])
    df["T2"] = df["T2"].astype(str).str.strip()
    df["T4"] = df["T4"].astype(str).str.strip()
    df["T3"] = pd.to_numeric(df["T3"], errors="coerce")
    return df[["T2", "T3", "T4", "Alpha Missense Score", "Alpha missense classification"]]


def _load_acmg_other_points(path: str, sheet: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=sheet)
    col_t2 = _find_col(df, ["T2"])
    col_t3 = _find_col(df, ["T3"])
    col_t4 = _find_col(df, ["T4"])
    col_pp3 = _find_col(df, ["ACMG PP3/BP4 in silico predictor points"])
    col_pm2 = _find_col(df, ["ACMG PM2 points"])
    col_pp1 = _find_col(df, ["ACMG PP1/BS4 segregation points"])
    col_ba1 = _find_col(df, ["ACMG BA1/BS1 allele freq points"])
    out = df[[col_t2, col_t3, col_t4, col_pp3, col_pm2, col_pp1, col_ba1]].copy()
    out = out.rename(
        columns={
            col_t2: "T2",
            col_t3: "T3",
            col_t4: "T4",
            col_pp3: "ACMG PP3/BP4 in silico predictor points",
            col_pm2: "ACMG PM2 points",
            col_pp1: "ACMG PP1/BS4 segregation points",
            col_ba1: "ACMG BA1/BS1 allele freq points",
        }
    )
    out["T2"] = out["T2"].astype(str).str.strip()
    out["T4"] = out["T4"].astype(str).str.strip()
    out["T3"] = pd.to_numeric(out["T3"], errors="coerce")
    return out


def _attach_new_columns(
    base_df: pd.DataFrame,
    alpha_df: pd.DataFrame,
    other_points_df: pd.DataFrame,
) -> pd.DataFrame:
    df = base_df.copy()
    df.columns = df.columns.astype(str).str.strip()

    col_votes = next((c for c in df.columns if "all_votes" in _norm_col(c)), None)
    if not col_votes:
        raise ValueError("Could not find an 'All_votes' column in Sup Table 12/13.")
    col_t2 = _find_col(df, ["T2"])
    col_t3 = _find_col(df, ["T3"])
    col_t4 = _find_col(df, ["T4"])

    df["Parsed_votes"] = df[col_votes].apply(safe_parse_votes)
    df["ACMG BS3/PS3 functional points computation"] = df["Parsed_votes"].apply(compute_acmg_score)
    df["Number of assays"] = df["Parsed_votes"].apply(len)
    df["ACMG BS3/PS3 Final Functional Points"] = df["ACMG BS3/PS3 functional points computation"].apply(compute_final_score)
    df["ACMG BS3/PS3 Capped Final Functional Points"] = df["ACMG BS3/PS3 Final Functional Points"].apply(cap_score)

    df = df.drop(columns=["Parsed_votes"])

    # Merge AlphaMissense
    df["T2_key"] = df[col_t2].astype(str).str.strip()
    df["T4_key"] = df[col_t4].astype(str).str.strip()
    df["T3_key"] = pd.to_numeric(df[col_t3], errors="coerce")
    df = df.merge(
        alpha_df,
        how="left",
        left_on=["T2_key", "T3_key", "T4_key"],
        right_on=["T2", "T3", "T4"],
        suffixes=("", "_alpha"),
    )

    # Merge ACMG other points
    df = df.merge(
        other_points_df,
        how="left",
        left_on=["T2_key", "T3_key", "T4_key"],
        right_on=["T2", "T3", "T4"],
        suffixes=("", "_other"),
    )

    drop_cols = ["T2_key", "T3_key", "T4_key", "T2_alpha", "T3_alpha", "T4_alpha", "T2_other", "T3_other", "T4_other"]
    for col in drop_cols:
        if col in df.columns:
            df = df.drop(columns=[col])

    # Final ACMG points + classification
    points_cols = [
        "ACMG PP3/BP4 in silico predictor points",
        "ACMG PM2 points",
        "ACMG PP1/BS4 segregation points",
        "ACMG BA1/BS1 allele freq points",
        "ACMG BS3/PS3 Capped Final Functional Points",
    ]
    for c in points_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    df["FINAL ACMG POINTS"] = df[points_cols].sum(axis=1)
    df["FINAL ClinVar Classification"] = df["FINAL ACMG POINTS"].apply(classify_final_points)

    # Order: keep original columns, then append new ones in required order
    append_cols = [
        "ACMG BS3/PS3 functional points computation",
        "Number of assays",
        "ACMG BS3/PS3 Final Functional Points",
        "ACMG BS3/PS3 Capped Final Functional Points",
        "Alpha Missense Score",
        "Alpha missense classification",
        "ACMG PP3/BP4 in silico predictor points",
        "ACMG PM2 points",
        "ACMG PP1/BS4 segregation points",
        "ACMG BA1/BS1 allele freq points",
        "FINAL ACMG POINTS",
        "FINAL ClinVar Classification",
    ]
    existing = [c for c in append_cols if c in df.columns]
    base_cols = [c for c in base_df.columns if c in df.columns]
    df = df[base_cols + existing]
    return df


def _write_df(output_path: str, sheet_name: str, df: pd.DataFrame) -> None:
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.exists():
        wb = load_workbook(output_file)
        if sheet_name in wb.sheetnames:
            wb.remove(wb[sheet_name])
        wb.save(output_file)
        with pd.ExcelWriter(output_file, engine="openpyxl", mode="a") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)


def write_sup_tables_18_19(
    sup12_df: pd.DataFrame,
    sup13_df: pd.DataFrame,
    alpha_path: str,
    other_points_path: str,
    output_path: str,
) -> None:
    alpha_b1 = _load_alpha_missense(alpha_path, "BRCA1")
    alpha_b2 = _load_alpha_missense(alpha_path, "BRCA2")
    other_b1 = _load_acmg_other_points(other_points_path, "BRCA1")
    other_b2 = _load_acmg_other_points(other_points_path, "BRCA2")

    out18 = _attach_new_columns(sup12_df, alpha_b1, other_b1)
    out19 = _attach_new_columns(sup13_df, alpha_b2, other_b2)

    _write_df(output_path, "Sup Table 18", out18)
    _write_df(output_path, "Sup Table 19", out19)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("workbook", help="Input workbook containing Sup Table 12/13")
    ap.add_argument("-o", "--out", default="SupTable18_19.xlsx", help="Output .xlsx")
    ap.add_argument("--sup12", default="Sup Table 12")
    ap.add_argument("--sup13", default="Sup Table 13")
    ap.add_argument(
        "--alpha",
        default="dataset/AlphaMissense_Calculations_all.xlsx",
        help="AlphaMissense workbook with BRCA1/BRCA2 tabs",
    )
    ap.add_argument(
        "--other-points",
        default="dataset/ACMG_other_points.xlsx",
        help="Workbook with additional ACMG point columns (BRCA1/BRCA2 tabs)",
    )
    args = ap.parse_args()

    sup12_sheet = resolve_sheet_name(args.workbook, args.sup12, ["sup table 12", "table 12"])
    sup13_sheet = resolve_sheet_name(args.workbook, args.sup13, ["sup table 13", "table 13"])

    sup12_df = pd.read_excel(args.workbook, sheet_name=sup12_sheet, header=1)
    sup13_df = pd.read_excel(args.workbook, sheet_name=sup13_sheet, header=1)

    write_sup_tables_18_19(sup12_df, sup13_df, args.alpha, args.other_points, args.out)
    print(f"[OK] Wrote Sup Table 18/19 to {args.out}")


if __name__ == "__main__":
    main()
