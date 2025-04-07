
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold

# Load the dataset
file_path = "./dataset/SUPP_TABLES_BRCA12_JAN_2025_V6.xlsx"
BRCA1_table = "Sup Table 1"
BRCA2_table = "Sup Table 2"
OUTPUT_NAME = "classification_results_10fold_V2.xlsx"

print("Reading Excel file...")
df_brca1 = pd.read_excel(file_path, sheet_name=BRCA1_table, header=1)
df_brca2 = pd.read_excel(file_path, sheet_name=BRCA2_table, header=1)
print("Excel file loaded successfully.")

# Function to calculate Wilson score interval
def wilson_score_interval(p, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    denominator = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denominator
    margin = z * np.sqrt((p * (1 - p) / n) + (z**2 / (4 * n**2))) / denominator
    lower_bound = max(0, center - margin)
    upper_bound = min(1, center + margin)
    return lower_bound, upper_bound

# Function to classify odds
def classify_odds(odds):
    if 0.0001 < odds < 0.0029:
        return "BS3_very_strong"
    elif odds < 0.053:
        return "BS3"
    elif odds < 0.23:
        return "BS3_moderate"
    elif odds < 0.48:
        return "BS3_supporting"
    elif odds > 350:
        return "PS3_very_strong"
    elif odds > 18.7:
        return "PS3"
    elif odds > 4.3:
        return "PS3_moderate"
    elif odds > 2.1:
        return "PS3_supporting"
    else:
        return "indeterminate"

# Function to compute final classification
def get_final_classification(classification_str):
    if pd.isna(classification_str) or classification_str.strip() == "":
        return np.nan
    classifications = [item.split(" (")[0] for item in classification_str.split("; ")]
    if all(c.startswith("BS3") for c in classifications):
        return 0
    elif all(c.startswith("PS3") for c in classifications):
        return 2
    elif all(c in ["indeterminate", "hypomorph"] for c in classifications):
        return 1
    else:
        bs3_count = sum(1 for c in classifications if c.startswith("BS3"))
        ps3_count = sum(1 for c in classifications if c.startswith("PS3"))
        if ps3_count >= 3 * bs3_count:
            return 2
        elif bs3_count >= 3 * ps3_count:
            return 0
        else:
            return 1

# Function to process each dataset
def process_table(df, table_name, writer):
    df = df[df["T6"].notna()].copy()
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    classification_results = df[["T1", "T2", "T3", "T4", "T5", "T6"]].copy()
    start_col_index = df.columns.get_loc("T8")

    for fold, (train_index, test_index) in enumerate(kf.split(df), 1):
        print(f"Processing fold {fold} for {table_name}...")
        train_df = df.iloc[train_index].copy()
        test_df = df.iloc[test_index].copy()
        specificity_sensitivity = {}

        for column in df.columns[start_col_index:]:
            tp = ((train_df[column] == 2) & (train_df["T6"].isin([4, 5, "4;5"]))).sum()
            tn = ((train_df[column] == 0) & (train_df["T6"].isin([1, 2, "1;2"]))).sum()
            fp = ((train_df[column] == 2) & ~train_df["T6"].isin([1, 2, "1;2"])).sum()
            fn = ((train_df[column] == 0) & ~train_df["T6"].isin([4, 5, "4;5"])).sum()

            sensitivity = tp / (tp + fn) if (tp + fn) > 0 else np.nan
            specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan

            P1_denominator = tp + tn + fp + fn
            P1 = (tp + fn) / P1_denominator if P1_denominator > 0 else np.nan
            P2_path = (tp + fp) / (tp + fp + 0.5)
            P2_benign = 0.5 / (tn + fn + 0.5)
            oddspath_path = (P2_path * (1 - P1)) / ((1 - P2_path) * P1) if (1 - P2_path) * P1 != 0 else np.nan
            oddspath_benign = (P2_benign * (1 - P1)) / ((1 - P2_benign) * P1) if (1 - P2_benign) * P1 != 0 else np.nan

            acmg_benign = classify_odds(oddspath_benign)
            acmg_path = classify_odds(oddspath_path)

            specificity_sensitivity[column] = {
                "acmg_benign": acmg_benign,
                "acmg_path": acmg_path
            }

        # Classify test set
        def classify_variant(row):
            classification = []
            for col in df.columns[start_col_index:]:
                if row[col] == 2:
                    classification.append(f"{specificity_sensitivity[col]['acmg_path']} ({col})")
                elif row[col] == 1:
                    classification.append(f"hypomorph ({col})")
                elif row[col] == 0:
                    classification.append(f"{specificity_sensitivity[col]['acmg_benign']} ({col})")
            return "; ".join(classification) if classification else np.nan

        classification_results.loc[test_df.index, f"Fold_{fold}"] = test_df.apply(classify_variant, axis=1)

    # Add final classification
    for fold in range(1, 11):
        classification_results[f"Final_Fold_{fold}"] = classification_results[f"Fold_{fold}"].apply(get_final_classification)

    # Save detailed and final classifications
    detailed_cols = ["T1", "T2", "T3", "T4", "T5", "T6"] + [f"Fold_{fold}" for fold in range(1, 11)]
    classification_results[detailed_cols].to_excel(writer, sheet_name=f"{table_name}_detailed", index=False)
    final_cols = ["T1", "T2", "T3", "T4", "T5", "T6"] + [f"Final_Fold_{fold}" for fold in range(1, 11)]
    classification_results[final_cols].to_excel(writer, sheet_name=f"{table_name}_final_classification", index=False)

    # Compute metrics with 95% CI
    metrics_list = []
    for fold in range(1, 11):
        valid_variants = classification_results[f"Fold_{fold}"].notna() & (classification_results[f"Fold_{fold}"] != "")
        final_class = classification_results.loc[valid_variants, f"Final_Fold_{fold}"]
        t6 = df.loc[valid_variants, "T6"]

        pathogenic = t6.isin([4, 5, "4;5"])
        benign = t6.isin([1, 2, "1;2"])

        tp = ((final_class == 2) & pathogenic).sum()
        tn = ((final_class == 0) & benign).sum()
        fp = ((final_class == 2) & benign).sum()
        fn = ((final_class == 0) & pathogenic).sum()

        no_class = (final_class == 1).sum()

        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else np.nan
        specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan

        # 95% CI for sensitivity
        n_sens = tp + fn
        p_sens = tp / n_sens if n_sens > 0 else 0
        sens_lower, sens_upper = wilson_score_interval(p_sens, n_sens)

        # 95% CI for specificity
        n_spec = tn + fp
        p_spec = tn / n_spec if n_spec > 0 else 0
        spec_lower, spec_upper = wilson_score_interval(p_spec, n_spec)

        metrics_list.append({
            "fold": fold,
            "sensitivity": sensitivity,
            "specificity": specificity,
            "sensitivity_lower_95CI": sens_lower,
            "sensitivity_upper_95CI": sens_upper,
            "specificity_lower_95CI": spec_lower,
            "specificity_upper_95CI": spec_upper,
            "true_positive": tp,
            "false_positive": fp,
            "true_negative": tn,
            "false_negative": fn,
            "no_classification_count": no_class
        })

    metrics_df = pd.DataFrame(metrics_list)
    metrics_df.to_excel(writer, sheet_name=f"{table_name}_k_fold_spec_sens", index=False)

# Process tables and save results
with pd.ExcelWriter(OUTPUT_NAME, engine="xlsxwriter") as writer:
    process_table(df_brca1, "BRCA1", writer)
    process_table(df_brca2, "BRCA2", writer)

print(f"Results saved to {OUTPUT_NAME}")

