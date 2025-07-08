import sys
import pandas as pd

def load_variant_positions(xlsx, sheet_name):
    # Sup Table 1 (BRCA1) or Sup Table 2 (BRCA2): columns T3 (amino acid index) and T5 (variant)
    df = pd.read_excel(xlsx, sheet_name=sheet_name, header=1, usecols=["T3", "T5"])
    df.columns = ["pos", "Variant"]
    df = df.dropna(subset=["pos", "Variant"]).copy()
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype(int)
    return df


def load_evidence(xlsx, sheet_name):
    # Sup Table 13 (BRCA1) or Sup Table 14 (BRCA2)
    df = pd.read_excel(
        xlsx,
        sheet_name=sheet_name,
        header=1,
        usecols=["T5", "Integrated ACMG evidence criteria", "Functional category", "Hypomorph observation"]
    )
    df.columns = ["Variant", "ACMG", "Class", "HypoObs"]
    df = df.dropna(subset=["Variant"])  # must have variant
    # Normalize
    df["Variant"] = df["Variant"].astype(str).str.strip()
    df["ACMG"] = df["ACMG"].astype(str).str.strip()
    df["Class"] = df["Class"].astype(str).str.strip()
    df["HypoObs"] = df["HypoObs"].astype(str).str.strip().str.upper()
    return df


def load_regions(xlsx, sheet_name, usecols):
    # Sup Table 17: region definitions; header row at index=1
    df = pd.read_excel(xlsx, sheet_name=sheet_name, header=1, usecols=usecols)
    # Rename columns to Feature, Start aa, End aa
    df.columns = ["Feature", "Start aa", "End aa"]
    df = df.dropna(subset=["Feature", "Start aa", "End aa"])  # ensure valid
    df["Start aa"] = pd.to_numeric(df["Start aa"], errors="coerce").astype(int)
    df["End aa"] = pd.to_numeric(df["End aa"], errors="coerce").astype(int)
    return df


def summarize_by_region(var_pos, evidence, regions):
    # Merge evidence with positions
    ev = pd.merge(evidence, var_pos, on="Variant", how="inner")
    records = []

    # categories to count
    acmg_list = [
        "BS3_supporting", "BS3_moderate", "BS3",
        "VUS",
        "PS3_supporting", "PS3_moderate", "PS3"
    ]

    for _, row in regions.iterrows():
        feat = row["Feature"]
        start = row["Start aa"]
        end = row["End aa"]
        sub = ev[(ev["pos"] >= start) & (ev["pos"] <= end)]
        total = len(sub)
        counts = {cat: sub["ACMG"].eq(cat).sum() for cat in acmg_list}
        hyp_obs = sub["HypoObs"].eq("Y").sum()
        path_obs = sub["Class"].eq("Pathogenic").sum()
        hypo_freq = hyp_obs / total if total > 0 else 0
        path_freq = path_obs / total if total > 0 else 0
        rec = {
            "Feature": feat,
            "Start aa": start,
            "End aa": end,
            **counts,
            "Total": total,
            "Hypomorph observations": hyp_obs,
            "Hypo observations frequency": round(hypo_freq, 3),
            "Pathogenic frequency": round(path_freq, 3)
        }
        records.append(rec)

    return pd.DataFrame.from_records(records)


def main(xlsx):
    # Load all four data pieces
    # BRCA1
    br1_pos = load_variant_positions(xlsx, "Sup Table 1")
    br1_ev  = load_evidence(xlsx, "Sup Table 13")
    br1_reg = load_regions(xlsx, "Sup Table 17", usecols="A:C")
    df1 = summarize_by_region(br1_pos, br1_ev, br1_reg)

    # BRCA2
    br2_pos = load_variant_positions(xlsx, "Sup Table 2")
    br2_ev  = load_evidence(xlsx, "Sup Table 14")
    # BRCA2 regions are in columns R:T of Sup Table 17
    br2_reg = load_regions(xlsx, "Sup Table 17", usecols="R:T")
    df2 = summarize_by_region(br2_pos, br2_ev, br2_reg)

    # Write out to Excel
    out = "by_region_summary.xlsx"
    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        df1.to_excel(writer, sheet_name="BRCA1", index=False)
        df2.to_excel(writer, sheet_name="BRCA2", index=False)

    print(f"Wrote region summary to {out}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <SUPP_TABLES.xlsx>")
        sys.exit(1)
    main(sys.argv[1])
