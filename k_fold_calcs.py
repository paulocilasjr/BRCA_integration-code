import pandas as pd
import numpy as np
from sklearn.model_selection import KFold

# Load the dataset
file_path = "./dataset/SUPP_TABLES_BRCA12_JAN_2025_V6.xlsx"
BRCA1_table = "Sup Table 1"
BRCA2_table = "Sup Table 2"

print("Reading Excel file...")
df_brca1 = pd.read_excel(file_path, sheet_name=BRCA1_table, header=1)
df_brca2 = pd.read_excel(file_path, sheet_name=BRCA2_table, header=1)
print("Excel file loaded successfully.")

# Function to classify odds
def classify_odds(odds):
    if odds < 0.0029:
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

# Function to process each dataset
def process_table(df, table_name, writer, sensitivity_specificity_dict):
    print(f"Processing {table_name}...")
    df = df[df["T6"].notna()].copy()

    kf = KFold(n_splits=10, shuffle=True, random_state=42)

    # Create a new DataFrame to store classification results with T1 to T6 columns
    classification_results = df[["T1", "T2", "T3", "T4", "T5", "T6"]].copy()

    start_col_index = df.columns.get_loc("T8")

    for fold, (train_index, test_index) in enumerate(kf.split(df), 1):
        print(f"Starting fold {fold}...")

        # Initialize metrics storage
        specificity_sensitivity = {}

        for column in df.columns[start_col_index:]:
            print(f"Calculating metrics for {column}...")
            train_df = df.iloc[train_index].copy()

            true_positive = ((train_df[column] == 2) & (train_df["T6"].isin([4, 5, "4;5"]))).sum()
            true_negative = ((train_df[column] == 0) & (train_df["T6"].isin([1, 2, "1;2"]))).sum()
            false_positive = ((train_df[column] == 2) & ~train_df["T6"].isin([1, 2, "1;2"])).sum()
            false_negative = ((train_df[column] == 0) & ~train_df["T6"].isin([4, 5, "4;5"])).sum()

            # Compute sensitivity & specificity
            sensitivity = true_positive / (true_positive + false_negative) if (true_positive + false_negative) > 0 else np.nan
            specificity = true_negative / (true_negative + false_positive) if (true_negative + false_positive) > 0 else np.nan

            # Store sensitivity and specificity for this fold and column
            if column not in sensitivity_specificity_dict:
                sensitivity_specificity_dict[column] = {}

            sensitivity_specificity_dict[column][f"fold_{fold}"] = {
                "sensitivity": sensitivity,
                "specificity": specificity
            }

            # Compute odds ratios
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

        # Define classify_variant function inside process_table
        def classify_variant(row):
            classification = []
            for col in df.columns[start_col_index:]:
                if row[col] == 2:
                    classification.append(specificity_sensitivity[col]["acmg_path"])
                elif row[col] == 1:
                    classification.append("hypomorph")
                elif row[col] == 0:
                    classification.append(specificity_sensitivity[col]["acmg_benign"])
            return "; ".join(classification)

        test_df = df.iloc[test_index].copy()
        classification_results.loc[test_df.index, f"Fold_{fold}"] = test_df.apply(classify_variant, axis=1)

    # Save the classification results into the Excel file
    classification_results.to_excel(writer, sheet_name=table_name, index=False)
    print(f"Processing {table_name} completed.")

# Writing results to a single Excel file
with pd.ExcelWriter("classification_results_2.xlsx", engine="xlsxwriter") as writer:
    sensitivity_specificity_dict_brca1 = {}
    sensitivity_specificity_dict_brca2 = {}

    process_table(df_brca1, "BRCA1", writer, sensitivity_specificity_dict_brca1)
    process_table(df_brca2, "BRCA2", writer, sensitivity_specificity_dict_brca2)

    # Create DataFrames for sensitivity and specificity per fold and per track for each dataset
    sensitivity_specificity_df_brca1 = pd.DataFrame.from_dict(sensitivity_specificity_dict_brca1, orient="index")
    sensitivity_specificity_df_brca2 = pd.DataFrame.from_dict(sensitivity_specificity_dict_brca2, orient="index")
    
    # Write sensitivity and specificity to separate sheets for BRCA1 and BRCA2
    sensitivity_specificity_df_brca1.to_excel(writer, sheet_name="BRCA1_Sensitivity_Specificity", index=True)
    sensitivity_specificity_df_brca2.to_excel(writer, sheet_name="BRCA2_Sensitivity_Specificity", index=True)

print("Processing complete. Results saved to classification_results.xlsx")

