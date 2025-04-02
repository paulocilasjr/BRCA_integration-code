import pandas as pd
import numpy as np
import sys  # Import sys to exit script if condition is not met
import re

from time import sleep

# Load the dataset
file_path = "./dataset/SUPP_TABLES_BRCA12_JAN_2025_V6_alpha.xlsx"
BRCA1_table = "Sup Table 1"
BRCA2_table = "Sup Table 2"
BRCA1_class = "Sup Table 10"
BRCA2_class = "Sup Table 11"

print("Reading Excel file...")
df_brca1 = pd.read_excel(file_path, sheet_name=BRCA1_table, header=1)
df_brca2 = pd.read_excel(file_path, sheet_name=BRCA2_table, header=1)
brca1_class = pd.read_excel(file_path, sheet_name=BRCA1_class, header=2)
brca2_class = pd.read_excel(file_path, sheet_name=BRCA2_class, header=2)
print("Excel file loaded successfully.")


# Function to classify odds
def classify_odds(odds):
    if pd.isna(odds):  # Check for NaN
        return None
    odds = round(odds, 9)  # Round to 2 decimal places before classification

    # Classify the odds value
    if odds == 0:
        return 0
    thresholds = [0.0029, 0.053, 0.231, 0.48, 2.08, 4.33, 18.7, 350]
    scores = [-8, -4, -2, -1, 0, 1, 2, 4]
    for threshold, score in zip(thresholds, scores):
        if odds <= threshold:
            return score
    return 8  # Default case: odds >= 350


# Function to iterate over rows, classify values, and store in dictionary
def get_scores(df):
    result_dict = {}

    for _, row in df.iterrows():
        key = row.iloc[0]  # First column as key
        value_0 = classify_odds(row.iloc[-3])  # Classify value from column -3
        value_1 = classify_odds(row.iloc[-4])  # Classify value from column -4

        if value_0 is not None and value_0 > 0:
            value_0 = 0
        if value_1 is not None and value_1 < 0:
            value_1 = 0

        result_dict[key] = {0: value_0, 1: value_1}

    return result_dict

scores_dict_brca1 = get_scores(brca1_class)
scores_dict_brca2 = get_scores(brca2_class)


def process_tracks(df, result_dict):
    """
    Iterate over all columns in the given dataframe starting from column "T8".
    For each row, retrieve values based on the result_dict mapping.

    Args:
        df (pd.DataFrame): DataFrame containing track data.
        result_dict (dict): Dictionary mapping track numbers (columns) to classification values.

    Returns:
        dict: Dictionary where each row index maps to an array of retrieved values.
    """
    start_col_index = df.columns.get_loc("T8")  # Find the index of the first track column
    row_results = {}  # Dictionary to store lists of values for each row index
    row_results_tracks = {}

    for index, row in df.iterrows():
        track_values = []  # List to store retrieved values for this row
        track_values_tracks = []

        for col in df.columns[start_col_index:]:  # Iterate over track columns
            if pd.notna(row[col]):  # Only consider non-empty values
                track_key = col  # Column name corresponds to track key in result_dict
                if track_key in result_dict:
                    if row[col] == 0:
                        track_values.append(result_dict[track_key][0])  # Retrieve value_0
                        track_values_tracks.append(f"{result_dict[track_key][0]} ({col})")
                    elif row[col] == 2:
                        track_values.append(result_dict[track_key][1])  # Retrieve value_1
                        track_values_tracks.append(f"{result_dict[track_key][1]} ({col})")
                    else:
                        track_values.append(0)  # Default value for other cases
                        track_values_tracks.append(f"0 ({col})")
        # Store the values in the dictionary under the row index
        if track_values:
            row_results[index] = track_values
            row_results_tracks[index] = track_values_tracks

    return row_results, row_results_tracks

# Process BRCA1 and BRCA2 datasets
brca1_results, brca1_results_tracks = process_tracks(df_brca1, scores_dict_brca1)
brca2_results, brca2_results_tracks = process_tracks(df_brca2, scores_dict_brca2)

def add_acmg_points(df, row_results, row_results_tracks):
    """
    Adds the ACMG Points System and Final Score columns to the DataFrame, keeping only relevant columns.

    Args:
        df (pd.DataFrame): Original DataFrame containing track data.
        row_results (dict): Dictionary with row index as key and computed ACMG values as list.

    Returns:
        pd.DataFrame: Filtered DataFrame with the ACMG Points System and Final Score columns added.
    """
    # Keep only columns up to "T7"
    relevant_columns = df.loc[:, : "T7"].copy()

    # Function to extract the numeric value from the string
    def extract_numeric_value(track_data):
        if isinstance(track_data, int):
            return track_data  # Directly return the value if it's an integer
        # Find all numbers before any parentheses
        match = re.findall(r'-?\d+', track_data)
        # Return the first numeric value found (it should be the ACMG points)
        return int(match[0]) if match else 0

    # Map ACMG Points System to the DataFrame
    relevant_columns["ACMG Points System"] = relevant_columns.index.map(row_results)

    # Compute Final Score (sum of ACMG Points System, extracting numbers from the track data)
    relevant_columns["Final Score"] = relevant_columns["ACMG Points System"].apply(
        lambda x: sum(extract_numeric_value(item) for item in x) if isinstance(x, list) else None
    )

    relevant_columns["ACMG Points System"] = relevant_columns.index.map(row_results_tracks)

    # Remove rows where ACMG Points System is empty (None or NaN)
    relevant_columns = relevant_columns.dropna(subset=["ACMG Points System"])

    return relevant_columns

# Apply function to both datasets
df_brca1_final = add_acmg_points(df_brca1, brca1_results, brca1_results_tracks)
df_brca2_final = add_acmg_points(df_brca2, brca2_results, brca2_results_tracks)

# Save the DataFrames as CSV files
df_brca1_final.to_csv("brca1_acmg_score_v4_alpha.csv", index=False)
df_brca2_final.to_csv("brca2_acmg_score_v4_alpha.csv", index=False)

print("CSV files saved successfully.")
