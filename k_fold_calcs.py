import pandas as pd
import numpy as np
from sklearn.model_selection import KFold

# Load the dataset
file_path = "./dataset/SUPP_TABLES_BRCA12_JAN_2025_V6.xlsx"
BRCA1_table = "Sup Table 1"
BRCA2_table = "Sup Table 2"

OUTPUT_NAME = "classification_results_3.xlsx"

print("Reading Excel file...")
df_brca1 = pd.read_excel(file_path, sheet_name=BRCA1_table, header=1)
df_brca2 = pd.read_excel(file_path, sheet_name=BRCA2_table, header=1)
print("Excel file loaded successfully.")

# Function to classify odds with updated BS3_very_strong condition
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

# Function to compute final classification based on ACMG classifications
def get_final_classification(classification_str):
    if pd.isna(classification_str):
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
def process_table(df, table_name, writer, sensitivity_specificity_dict):
    print(f"Processing {table_name}...")
    df = df[df["T6"].notna()].copy()

    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    classification_results = df[["T1", "T2", "T3", "T4", "T5", "T6"]].copy()
    start_col_index = df.columns.get_loc("T8")
    fold_variants = []

    for fold, (train_index, test_index) in enumerate(kf.split(df), 1):
        print(f"Starting fold {fold}...")
        specificity_sensitivity = {}

        # Extract variant names from T5 column as lists
        train_variants = df.iloc[train_index]["T5"].astype(str).tolist()
        test_variants = df.iloc[test_index]["T5"].astype(str).tolist()
        fold_variants.append({
            "fold": fold,
            "train_variants": train_variants,
            "test_variants": test_variants
        })

        for column in df.columns[start_col_index:]:
            print(f"Calculating metrics for {column}...")
            train_df = df.iloc[train_index].copy()

            true_positive = ((train_df[column] == 2) & (train_df["T6"].isin([4, 5, "4;5"]))).sum()
            true_negative = ((train_df[column] == 0) & (train_df["T6"].isin([1, 2, "1;2"]))).sum()
            false_positive = ((train_df[column] == 2) & ~train_df["T6"].isin([1, 2, "1;2"])).sum()
            false_negative = ((train_df[column] == 0) & ~train_df["T6"].isin([4, 5, "4;5"])).sum()

            sensitivity = true_positive / (true_positive + false_negative) if (true_positive + false_negative) > 0 else np.nan
            specificity = true_negative / (true_negative + false_positive) if (true_negative + false_positive) > 0 else np.nan

            if column not in sensitivity_specificity_dict:
                sensitivity_specificity_dict[column] = {}
            sensitivity_specificity_dict[column][f"fold_{fold}"] = {
                "sensitivity": sensitivity,
                "specificity": specificity
            }

            P1_denominator = true_positive + true_negative + false_positive + false_negative
            P1 = (true_positive + false_negative) / P1_denominator if P1_denominator > 0 else np.nan
            P2_path = (true_positive + false_positive) / (true_positive + false_positive + 0.5)
            P2_benign = 0.5 / (true_negative + false_negative + 0.5)
            oddspath_path = (P2_path * (1 - P1)) / ((1 - P2_path) * P1) if (1 - P2_path) * P1 != 0 else np.nan
            oddspath_benign = (P2_benign * (1 - P1)) / ((1 - P2_benign) * P1) if (1 - P2_benign) * P1 != 0 else np.nan

            acmg_benign = classify_odds(oddspath_benign)
            acmg_path = classify_odds(oddspath_path)

            specificity_sensitivity[column] = {
                "sensitivity": sensitivity,
                "specificity": specificity,
                "P1": P1,
                "P2_path": P2_path,
                "P2_benign": P2_benign,
                "oddspath_path": oddspath_path,
                "oddspath_benign": oddspath_benign,
                "acmg_benign": acmg_benign,
                "acmg_path": acmg_path
            }

        print("Classifying test set...")
        def classify_variant(row):
            classification = []
            for col in df.columns[start_col_index:]:
                if row[col] == 2:
                    classification.append(f"{specificity_sensitivity[col]['acmg_path']} ({col})")
                elif row[col] == 1:
                    classification.append(f"hypomorph ({col})")
                elif row[col] == 0:
                    classification.append(f"{specificity_sensitivity[col]['acmg_benign']} ({col})")
            return "; ".join(classification)

        test_df = df.iloc[test_index].copy()
        classification_results.loc[test_df.index, f"Fold_{fold}"] = test_df.apply(classify_variant, axis=1)

    # Add final classification columns
    for fold in range(1, 11):
        classification_results[f"Final_Fold_{fold}"] = classification_results[f"Fold_{fold}"].apply(get_final_classification)

    # Save detailed classifications
    detailed_cols = ["T1", "T2", "T3", "T4", "T5", "T6"] + [f"Fold_{fold}" for fold in range(1, 11)]
    classification_results[detailed_cols].to_excel(writer, sheet_name=f"{table_name}_detailed", index=False)

    # Save final classifications
    final_cols = ["T1", "T2", "T3", "T4", "T5", "T6"] + [f"Final_Fold_{fold}" for fold in range(1, 11)]
    final_classification_df = classification_results[final_cols]
    final_classification_df.to_excel(writer, sheet_name=f"{table_name}_final_classification", index=False)

    # Compute and save metrics
    metrics_list = []
    for fold in range(1, 11):
        test_variants = classification_results[f"Final_Fold_{fold}"].notna()
        final_class = classification_results.loc[test_variants, f"Final_Fold_{fold}"]
        t6 = df.loc[test_variants, "T6"]

        pathogenic = t6.isin([4, 5, "4;5"])
        benign = t6.isin([1, 2, "1;2"])

        tp = ((final_class == 2) & pathogenic).sum()
        tn = ((final_class == 0) & benign).sum()
        fp = ((final_class == 2) & ~pathogenic).sum()
        fn = ((final_class == 0) & ~benign).sum()

        fn += ((final_class == 1) & pathogenic).sum()
        fp += ((final_class == 1) & benign).sum()

        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else np.nan
        specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan

        metrics_list.append({
            "fold": fold,
            "specificity": specificity,
            "sensitivity": sensitivity,
            "true_positive": tp,
            "false_positive": fp,
            "true_negative": tn,
            "false_negative": fn
        })

    metrics_df = pd.DataFrame(metrics_list)
    metrics_df.to_excel(writer, sheet_name=f"{table_name}_k_fold_spec_sens", index=False)

    # Save fold variants with T5 names in separate rows
    max_train_length = max(len(fold["train_variants"]) for fold in fold_variants)
    max_test_length = max(len(fold["test_variants"]) for fold in fold_variants)
    max_length = max(max_train_length, max_test_length)

    fold_variants_data = {}
    for fold_data in fold_variants:
        fold = fold_data["fold"]
        train_variants = fold_data["train_variants"] + [""] * (max_length - len(fold_data["train_variants"]))
        test_variants = fold_data["test_variants"] + [""] * (max_length - len(fold_data["test_variants"]))
        fold_variants_data[f"Fold_{fold}_Train"] = train_variants
        fold_variants_data[f"Fold_{fold}_Test"] = test_variants

    fold_variants_df = pd.DataFrame(fold_variants_data)
    fold_variants_df.to_excel(writer, sheet_name=f"{table_name}_k_fold_variants", index=False)

    print(f"Processing {table_name} completed.")

# Writing results to a single Excel file
with pd.ExcelWriter(OUTPUT_NAME, engine="xlsxwriter") as writer:
    sensitivity_specificity_dict_brca1 = {}
    sensitivity_specificity_dict_brca2 = {}

    process_table(df_brca1, "BRCA1", writer, sensitivity_specificity_dict_brca1)
    process_table(df_brca2, "BRCA2", writer, sensitivity_specificity_dict_brca2)

    # Save sensitivity and specificity DataFrames
    sensitivity_specificity_df_brca1 = pd.DataFrame.from_dict(sensitivity_specificity_dict_brca1, orient="index")
    sensitivity_specificity_df_brca2 = pd.DataFrame.from_dict(sensitivity_specificity_dict_brca2, orient="index")
    sensitivity_specificity_df_brca1.to_excel(writer, sheet_name="BRCA1_Sensitivity_Specificity", index=True)
    sensitivity_specificity_df_brca2.to_excel(writer, sheet_name="BRCA2_Sensitivity_Specificity", index=True)

print(f"Processing complete. Results saved to {OUTPUT_NAME}")
