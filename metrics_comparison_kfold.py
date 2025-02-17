import pandas as pd
import numpy as np
from sklearn.model_selection import KFold

# Load the dataset
file_path = "./dataset/SUPP_TABLES_BRCA12_JAN_2025_V6.xlsx"
BRCA1_table = "Sup Table 1"
BRCA2_table = "Sup Table 2"
BRCA1_metrics_table = "Sup Table 10"
BRCA2_metrics_table = "Sup Table 11"

print("Reading Excel file...")
df_brca1 = pd.read_excel(file_path, sheet_name=BRCA1_table, header=1)
df_brca2 = pd.read_excel(file_path, sheet_name=BRCA2_table, header=1)
df_brca1_metrics = pd.read_excel(file_path, sheet_name=BRCA1_metrics_table, header=3)
df_brca2_metrics = pd.read_excel(file_path, sheet_name=BRCA2_metrics_table, header=3)
print("Excel file loaded successfully.")

# Define classification ranking
classification_ranking = {
    "BS3_very_strong": 1, "BS3": 2, "BS3_moderate": 3, "BS3_supporting": 4, "indeterminate": 5,
    "PS3_very_strong": 6, "PS3": 7, "PS3_moderate": 8, "PS3_supporting": 9
}

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

# Function to compare classifications
def compare_classifications(calculated, reference):
    calculated_rank = classification_ranking.get(calculated, float('inf'))
    reference_rank = classification_ranking.get(reference, float('inf'))
    if calculated_rank < reference_rank:
        return "upgrade"
    elif calculated_rank > reference_rank:
        return "downgrade"
    else:
        return "same"

# Function to process each dataset
def process_table(df, metrics_df, table_name, writer, comparison_results):
    print(f"Processing {table_name}...")
    df = df[df["T6"].notna()].copy()
    metrics_columns = [metrics_df.columns[0], metrics_df.columns[1], metrics_df.columns[9], metrics_df.columns[12], metrics_df.columns[-2], metrics_df.columns[-1]]
    df_metrics = metrics_df[metrics_columns].copy()
    
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    classification_results = df[["T1", "T2", "T3", "T4", "T5", "T6"]].copy()
    start_col_index = df.columns.get_loc("T8")

    for fold, (train_index, test_index) in enumerate(kf.split(df), 1):
        print(f"Starting fold {fold}...")
        
        for i, column in enumerate(df.columns[start_col_index:]):
            print(f"Calculating metrics for {column}...")
            train_df = df.iloc[train_index].copy()
            
            true_positive = ((train_df[column] == 2) & (train_df["T6"].isin([4, 5, "4;5"]))).sum()
            true_negative = ((train_df[column] == 0) & (train_df["T6"].isin([1, 2, "1;2"]))).sum()
            false_positive = ((train_df[column] == 2) & ~train_df["T6"].isin([1, 2, "1;2"])).sum()
            false_negative = ((train_df[column] == 0) & ~train_df["T6"].isin([4, 5, "4;5"])).sum()

            sensitivity = true_positive / (true_positive + false_negative) if (true_positive + false_negative) > 0 else np.nan
            specificity = true_negative / (true_negative + false_positive) if (true_negative + false_positive) > 0 else np.nan
            
            P1_denominator = true_positive + true_negative + false_positive + false_negative
            P1 = (true_positive + false_negative) / P1_denominator if P1_denominator > 0 else np.nan

            P2_path = (true_positive + false_positive) / (true_positive + false_positive + 0.5)
            P2_benign = 0.5 / (true_negative + false_negative + 0.5)

            oddspath_path = (P2_path * (1 - P1)) / ((1 - P2_path) * P1) if (1 - P2_path) * P1 != 0 else np.nan
            oddspath_benign = (P2_benign * (1 - P1)) / ((1 - P2_benign) * P1) if (1 - P2_benign) * P1 != 0 else np.nan

            acmg_benign = classify_odds(oddspath_benign)
            acmg_path = classify_odds(oddspath_path)

            track = df_metrics.iloc[i, 0]
            ref_acmg_path = df_metrics.iloc[i, -2]
            ref_acmg_benign = df_metrics.iloc[i, -1]
            
            if track not in comparison_results:
                comparison_results[track] = {}
            
            comparison_results[track][f"fold_{fold}"] = {
                "acmg_path": compare_classifications(acmg_path, ref_acmg_path),
                "acmg_benign": compare_classifications(acmg_benign, ref_acmg_benign)
            }

    pd.DataFrame.from_dict(comparison_results, orient="index").to_excel(writer, sheet_name=f"{table_name}_Comparison")
    print(f"Processing {table_name} completed.")

# Writing results to a single Excel file
with pd.ExcelWriter("classification_results.xlsx", engine="xlsxwriter") as writer:
    comparison_results_brca1 = {}
    comparison_results_brca2 = {}

    process_table(df_brca1, df_brca1_metrics, "BRCA1", writer, comparison_results_brca1)
    process_table(df_brca2, df_brca2_metrics, "BRCA2", writer, comparison_results_brca2)

print("Processing complete. Results saved to classification_results.xlsx")
