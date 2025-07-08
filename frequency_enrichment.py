import pandas as pd
from collections import Counter

# --- Load the Excel file ---
excel_file = "SUPP_TABLES_BRCA12_JULY_2025_V2.xlsx"  # Replace with actual path
sup1 = pd.read_excel(excel_file, sheet_name="Sup Table 1")
sup13 = pd.read_excel(excel_file, sheet_name="Sup Table 13")

# --- Step 1: Get variant classes from Sup Table 13 ---
variant_class_map = sup13[["T5", "Functional category"]].dropna()
variant_class_map.columns = ["variant", "class"]
variant_class_map["class"] = variant_class_map["class"].str.strip().str.capitalize()

# --- Step 2: Filter Sup Table 1 ---
# Match variants from Sup13 and T6 must be empty
vus_variants = sup1[sup1["T5"].isin(variant_class_map["variant"]) & sup1["T6"].isna()]
vus_variants = vus_variants[["T5"]]  # Keep variant list only

# Merge class labels
vus_classified = pd.merge(vus_variants, variant_class_map, left_on="T5", right_on="variant", how="left")

# Add TYPE column from amino acid substitution (e.g., "p.M1I" -> "MI")
vus_classified["AA_SUB"] = sup1.set_index("T5").loc[vus_classified["T5"], "T4"].values
vus_classified = vus_classified[vus_classified["AA_SUB"].str.startswith("p.")]
vus_classified["TYPE"] = vus_classified["AA_SUB"].str.extract(r"p\.([A-Z])[0-9]+([A-Z])").agg(''.join, axis=1)

# Function to calculate frequency table
def calc_freq(df, total_count):
    counts = Counter(df["TYPE"])
    freq_df = pd.DataFrame(counts.items(), columns=["TYPE", "counts"])
    freq_df["f"] = freq_df["counts"] / total_count
    return freq_df.sort_values("TYPE").reset_index(drop=True)

# All VUS
vus_total = len(vus_classified)
vus_all = calc_freq(vus_classified, vus_total)

# By class
def classify_subset(class_name):
    subset = vus_classified[vus_classified["class"] == class_name]
    return calc_freq(subset, len(subset)), len(subset)

path_df, path_n = classify_subset("Pathogenic")
benign_df, benign_n = classify_subset("Benign")
hypo_df, hypo_n = classify_subset("Hypomorph")

# Merge and compute enrichment ratios
final = vus_all.copy()
final.columns = ["TYPE", "counts_all", "f_all"]

# Merge class-specific
for label, df, n, col in [
    ("Path", path_df, path_n, "ER (Path/all VUS)"),
    ("Benign", benign_df, benign_n, "ER (Benign/all VUS)"),
    ("Hyp", hypo_df, hypo_n, "ER (Hyp/all VUS)")
]:
    df = df.rename(columns={"counts": f"counts_{label}", "f": f"f_{label}"})
    final = final.merge(df, on="TYPE", how="left")
    final[col] = final[f"f_{label}"] / final["f_all"]

# Fill NaNs with 0 for counts and f columns
final.fillna({"counts_Path": 0, "f_Path": 0, "counts_Benign": 0, "f_Benign": 0,
              "counts_Hyp": 0, "f_Hyp": 0}, inplace=True)

# Optional: Round
final = final.round(3)

# Save to Excel
final.to_excel("BRCA1_AA_Substitution_Frequency_and_ER.xlsx", index=False)
