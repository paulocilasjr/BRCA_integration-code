import pandas as pd
import math
import re

# === Hierarchies ===
CLASSIFICATION_HIERARCHY = {
    "BS3": 5, "BS3_moderate": 4, "BS3_supporting": 3,
    "PS3": 5, "PS3_moderate": 4, "PS3_supporting": 3,
    "Hypomorph": 2,
    "Indeterminate": 1
}
PS3_HIERARCHY = CLASSIFICATION_HIERARCHY.copy()
BS3_HIERARCHY = CLASSIFICATION_HIERARCHY.copy()

# === Utility Functions ===
def get_strongest_classification(array, categorie=None):
    strongest = "Indeterminate"
    max_strength = 1
    hierarchy = {
        None: CLASSIFICATION_HIERARCHY,
        "PS3": PS3_HIERARCHY,
        "BS3": BS3_HIERARCHY,
    }[categorie]
    for item in array:
        if item in hierarchy and hierarchy[item] > max_strength:
            max_strength = hierarchy[item]
            strongest = item
    return strongest


def is_valid_ratio(var1, var2, class1, class2):
    if var1 == 0 or var2 == 0:
        return var1, var2, "Indeterminate"
    if var1 / var2 >= 3:
        return var1, var2, class1
    if var2 / var1 >= 3:
        return var2, var1, class2
    return var1, var2, "Indeterminate"

# === Discordance Logic ===
def check_discordance(array):
    bs3_count = ps3_count = hypomorph_count = 0
    cleaned = []
    for item in array:
        plain = re.sub(r"\(T\d+\)", "", item)
        cleaned.append(plain)
        if "BS3" in plain: bs3_count += 1
        if "PS3" in plain: ps3_count += 1
        if "hypomorph" in plain.lower(): hypomorph_count += 1

    bs3_present = bs3_count > 0
    ps3_present = ps3_count > 0
    hypomorph_present = hypomorph_count > 0
    strongest_call = get_strongest_classification(cleaned)
    hypo_obs = "Y" if hypomorph_present else "N"

    # PS3 vs BS3
    if bs3_present and ps3_present:
        if bs3_count > ps3_count:
            var1, var2, c1, c2 = bs3_count, ps3_count, "Benign", "Pathogenic"
            strongest = get_strongest_classification(cleaned, "BS3")
        else:
            var1, var2, c1, c2 = ps3_count, bs3_count, "Pathogenic", "Benign"
            strongest = get_strongest_classification(cleaned, "PS3")
        _, _, prepon = is_valid_ratio(var1, var2, c1, c2)
        if prepon == "Indeterminate":
            return "Discordant", prepon, "VUS", "Indeterminate", hypo_obs
        return "Discordant", prepon, strongest, prepon, hypo_obs

    # BS3 vs Hypomorph
    if bs3_present and hypomorph_present:
        if bs3_count > hypomorph_count:
            var1, var2, c1, c2 = bs3_count, hypomorph_count, "Benign", "Hypomorph"
        else:
            var1, var2, c1, c2 = hypomorph_count, bs3_count, "Hypomorph", "Benign"
        _, _, prepon = is_valid_ratio(var1, var2, c1, c2)
        if prepon == "Indeterminate":
            return "Discordant", prepon, "VUS", "Indeterminate", hypo_obs
        return "Discordant", prepon, strongest_call, prepon, hypo_obs

    # PS3 vs Hypomorph
    if ps3_present and hypomorph_present:
        if ps3_count > hypomorph_count:
            var1, var2, c1, c2 = ps3_count, hypomorph_count, "Pathogenic", "Hypomorph"
        else:
            var1, var2, c1, c2 = hypomorph_count, ps3_count, "Hypomorph", "Pathogenic"
        _, _, prepon = is_valid_ratio(var1, var2, c1, c2)
        if prepon == "Indeterminate":
            return "Discordant", prepon, "VUS", "Indeterminate", hypo_obs
        return "Discordant", prepon, strongest_call, prepon, hypo_obs

    # Concordant/hypomorph-only
    if bs3_present: return "Concordant", "-", strongest_call, "Benign", hypo_obs
    if ps3_present: return "Concordant", "-", strongest_call, "Pathogenic", hypo_obs
    if hypomorph_present: return "Concordant", "-", "VUS", "Hypomorph", hypo_obs

    return "Indeterminate", "Indeterminate", "VUS", "Indeterminate", "N"

# === Build Mapping ===
def BuildDict(tab_data):
    d = {}
    for _, row in tab_data.iterrows():
        key = str(row.iloc[0]).strip()
        d.setdefault(key, {})
        val_b = row.iloc[1]
        val_c = row.iloc[2]
        if isinstance(val_c, str) and val_c in CLASSIFICATION_HIERARCHY: d[key]["0"] = val_c
        if isinstance(val_b, str) and val_b in CLASSIFICATION_HIERARCHY: d[key]["2"] = val_b
    return d

# === Merge and Export ===
def merge_with_metadata(out_dict, metadata):
    rows = []
    for idx, arr in out_dict.items():
        conc, pre, fc, notes, hypo = check_discordance(arr)
        rows.append((idx, arr, conc, pre, fc, notes, hypo))
    df_out = pd.DataFrame(rows, columns=[
        "index", "All_votes(track)", "Concordance", "Preponderance of evidence",
        "Final code", "Notes", "Hypomorph observation"
    ])
    df_out.set_index("index", inplace=True)
    df_out.index = df_out.index.astype(metadata.index.dtype)
    return metadata.merge(df_out, left_index=True, right_index=True, how="left")

# === Main Results ===
def GetResults(df, df2, metadata, name):
    mapping = BuildDict(df2)
    out = {}
    for col in df.columns:
        for idx, val in df[col].items():
            raw = str(val).strip()
            if not raw: continue
            try:
                assay = int(raw)
            except ValueError:
                continue
            key = idx
            if assay == 1:
                cls = "hypomorph"
            else:
                cls = mapping.get(str(assay), {}).get("2") or mapping.get(str(assay), {}).get("0")
                if not cls: continue
            out.setdefault(key, []).append(f"{cls}({col})")
    merged = merge_with_metadata(out, metadata)
    filt = merged[merged["All_votes(track)"].notna() & merged["All_votes(track)"].astype(bool)]
    sel = pd.concat([filt.iloc[:, :7], filt.iloc[:, -6:]], axis=1)
    sel.to_csv(f"merged_output_12_{name}.csv")

file_path = "SUPP_TABLES_BRCA12_JUNE_2025_V2.xlsx"
sheet_name_BRCA1_table = "Sup Table 1"
sheet_name_BRCA2_table = "Sup Table 2"
sheet_name_BRCA1_class = "Sup Table 10"
sheet_name_BRCA2_class = "Sup Table 11"

brca1_df_full = pd.read_excel(file_path, sheet_name=sheet_name_BRCA1_table, engine="openpyxl", header=1)
brca2_df_full = pd.read_excel(file_path, sheet_name=sheet_name_BRCA2_table, engine="openpyxl", header=1)
brca1_class = pd.read_excel(file_path, sheet_name=sheet_name_BRCA1_class, engine="openpyxl", header=2)
brca2_class = pd.read_excel(file_path, sheet_name=sheet_name_BRCA2_class, engine="openpyxl", header=2)

def CleanDf(df, type=None):
    if type == "table":
        metadata_df = df.iloc[:, :]
        start_col_index = df.columns.get_loc('T8')
        df = df.iloc[:, start_col_index:]
        return df, metadata_df
    else:
        df_class = df.iloc[:, [0, -2, -1]]
        return df_class

brca1_df, brca1_metadata_df = CleanDf(brca1_df_full, "table")
brca2_df, brca2_metadata_df = CleanDf(brca2_df_full, "table")
brca1_class = CleanDf(brca1_class)
brca2_class = CleanDf(brca2_class)

GetResults(brca1_df, brca1_class, brca1_metadata_df, "BRCA1")
GetResults(brca2_df, brca2_class, brca2_metadata_df, "BRCA2")
