#!/usr/bin/env python3
import sys
import re
import pandas as pd

def extract_code(var):
    """
    From a string like "p.A123T" return "AT" (original aa + mutant aa).
    """
    m = re.match(r"p\.([A-Z])\d+([A-Z])", var)
    return m.group(1) + m.group(2) if m else None

def make_df(vus_df, variant_set, suffix):
    """
    Build a counts+frequency DataFrame for the given subset of variants.
    """
    sub = vus_df[vus_df["T5"].isin(variant_set)].copy()
    sub["TYPE"] = sub["T5"].map(extract_code)
    counts = sub["TYPE"].value_counts().sort_index().rename(f"counts_{suffix}")
    freqs = (counts / len(sub)).round(3).rename(f"f_{suffix}")
    return pd.concat([counts, freqs], axis=1)

def add_enrichment(df, ref_freq, suffix):
    """
    Given a df with f_<suffix>, divide by f_all to get ER_<suffix>.
    """
    df[f"ER_{suffix}"] = (
        df[f"f_{suffix}"] /
        df.index.to_series().map(lambda t: ref_freq.get(t, 1e-6))
    ).round(3)
    return df

def main(path):
    # 1) Load Sup Table 13 (skip the first junk row, real header is row 1)
    sup13 = pd.read_excel(path, sheet_name="Sup Table 14", header=1)
    sup13 = sup13[["T5", "Functional category"]].dropna(subset=["T5", "Functional category"])

    # build variant sets by functional category
    patho_vars  = set(sup13.loc[sup13["Functional category"]
                                 .str.contains("Pathogenic", case=False), "T5"])
    hypo_vars   = set(sup13.loc[sup13["Functional category"]
                                 .str.contains("Hypomorph", case=False), "T5"])
    benign_vars = (set(sup13.loc[sup13["Functional category"]
                                 .str.contains("Benign", case=False), "T5"])
                   - hypo_vars)
    all_vars    = set(sup13["T5"])

    # 2) Load Sup Table 1, skip first junk row, pick only T5 and T6
    sup1 = pd.read_excel(path, sheet_name="Sup Table 2", header=1)
    sup1 = sup1[["T5", "T6"]]

    # keep only variants in our classification table…
    sup1 = sup1[sup1["T5"].isin(all_vars)]
    # …and drop any that have a non‐empty T6
    vus  = sup1[sup1["T6"].isna()]

    # 3) All VUS
    df_all = make_df(vus, all_vars, "all")

    # 4) Pathogenic only
    df_path = make_df(vus, patho_vars, "path")

    # 5) Benign only
    df_benign = make_df(vus, benign_vars, "benign")

    # 6) Hypomorph only
    df_hypo = make_df(vus, hypo_vars, "hypo")

    # compute enrichment ratios relative to f_all
    ref_freq = df_all["f_all"].to_dict()
    df_path    = add_enrichment(df_path,    ref_freq, "path")
    df_benign  = add_enrichment(df_benign,  ref_freq, "benign")
    df_hypo    = add_enrichment(df_hypo,    ref_freq, "hypo")

    # merge everything into one wide table
    result = (
        df_all
        .merge(df_path,    left_index=True, right_index=True, how="outer")
        .merge(df_benign,  left_index=True, right_index=True, how="outer")
        .merge(df_hypo,    left_index=True, right_index=True, how="outer")
        .reset_index()
        .rename(columns={"index":"TYPE"})
        .fillna(0)
    )

    # write out and print
    out_file = "frequency_enrichment.xlsx"
    result.to_excel(out_file, index=False)
    print(f"Wrote enrichment table to {out_file}\n")
    print(result)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python frequency_enrichment.py <path_to_workbook.xlsx>")
        sys.exit(1)
    main(sys.argv[1])
